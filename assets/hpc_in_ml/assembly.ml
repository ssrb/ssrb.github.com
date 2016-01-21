(*  Copyright (c) 2015, Sebastien Sydney Robert Bigot
 All rights reserved. *)
open Core.Std
open Ctypes
open Foreign

(* Used for interop/FFI *)
open Bigarray

(* Reads a triangulation saved in the bamg mesh format *)
let read_mesh_exn mesh = 
	let open In_channel in
	with_file mesh ~f: (fun chan ->
		let module Parser = struct
			type state = V1 | V2 | V3 of int * int | T1 | T2 | T3 of int * int
			let fold call (state, vs, ts) line =
				match state with
				(* Vertices *)
				| V1 -> (if line = "Vertices" then V2 else V1), vs, ts
				| V2 -> let nvertices = Scanf.sscanf line "%d" (fun n -> n) in 
					(V3(0,nvertices), Some (Array2.create float64 c_layout nvertices 2), ts)
				| V3 (i,n) when i < n -> 
					Scanf.sscanf line "%f %f %f" (fun x y _ -> 
						match vs with 
						| Some vs' -> 
							vs'.{i,0} <- x; 
							vs'.{i,1} <- y
						| None -> ()
					);
					(V3(i+1,n), vs, ts)						
				| V3 _ -> (T1, vs, ts)
				(* Triangles *)
				| T1 ->  (if line = "Triangles" then T2 else T1), vs, ts
				| T2 -> let ntriangles = Scanf.sscanf line "%d" (fun n -> n) in	
					(T3(0,ntriangles), vs, Some (Array2.create int32 c_layout ntriangles 3))
				| T3(i,n) when i < n ->
					Scanf.sscanf line "%ld %ld %ld %ld" (fun v1 v2 v3 _ -> 
						match ts with 
						| Some ts' -> Int32.(
							ts'.{i,0} <- pred v1;
							ts'.{i,1} <- pred v2; 
							ts'.{i,2} <- pred v3
						)
						| None -> ()
					);
					(T3(i+1,n), vs, ts)
				| T3 _ -> call.return(state, vs, ts)
			let init = (V1, None, None)
		end
		in match with_return (fun call ->
			chan |> fold_lines ~init:Parser.init ~f:(Parser.fold call)) with
		| _, Some vs, Some ts -> vs, ts
		| _ -> failwith "Error parsing mesh file"
	)
;;

(* Builds the reverse connectivity table of the mesh, ie a vertex to triangle mapping as an array of list *)
let build_reverse_connectivity_table ts nv = 
	let rct = Array.create ~len:nv [] in
	for t = 0 to (Array2.dim1 ts) - 1 do
		let v0, v1, v2 = Int32.(to_int_exn ts.{t,0}, to_int_exn ts.{t,1}, to_int_exn ts.{t,2}) in
		rct.(v0) <- t::rct.(v0); 
		rct.(v1) <- t::rct.(v1); 
		rct.(v2) <- t::rct.(v2)
	done;
	rct
;;

(* Builds a perfect coloring of a cordal graph, the result is a BigArray as it's meant to be pushed to the GPU *)
let build_perfect_coloring ts rct =

	(* Builds the dual graph of a mesh, ie each node of the dual is a triangle of the mesh *)
	let build_triangulation_dual_graph ts rct = 	
		let nt = Array2.dim1 ts in
		let graph = Array.create ~len:nt [] in
		for t = 0 to (nt-1) do
			let v0, v1, v2 = Int32.(to_int_exn ts.{t,0}, to_int_exn ts.{t,1}, to_int_exn ts.{t,2}) in
			graph.(t) <- List.dedup (rct.(v0) @ rct.(v1) @ rct.(v2))
		done;
		graph

	(* Lex-BFS *)
	and build_perfect_elemination_ordering graph =	
		let compare (p1,_) (p2,_) = compare p2 p1 in
		let heap = Heap.Removable.create ~min_size:(Array.length graph) ~cmp:compare () in
		let elems = Array.mapi (fun v _ -> Heap.Removable.add_removable heap (0,v)) graph in	
		let rec aux peo = 
			match Heap.Removable.pop heap with
			| Some (p,v) -> List.iter graph.(v) ~f:(fun v' -> 
					try					
						let (p',_) = Heap.Removable.Elt.value_exn elems.(v') in
						elems.(v') <- Heap.Removable.update heap elems.(v') (succ p', v')
					with
					| _ -> ();
				);
				aux (v::peo)
			| None -> List.rev peo
		in 
		aux []

	(* Performs the greedy coloring of a graph (hopefuly cordal) walking the vertices in the provided order *)
	and do_greedy_coloring graph ordering =
		let n = Array.length graph in 
		let colors = Array1.create int32 c_layout n in
		let open Int32 in
		let perfect_color v =			
			let rec aux color = function
				| c::cs' -> if color < c then color else aux (succ color) cs'
				| [] -> color
			in graph.(v) |> List.map ~f:(fun v -> colors.{v}) |> List.sort ~cmp:compare |> aux zero
		in
		Array1.fill colors (of_int_exn n);
		List.iter ordering ~f:(fun v -> 	
			colors.{v} <- perfect_color v
		);
		colors
	in 

	(* As per Tarjan, we get the perfect coloring *)
	let graph = build_triangulation_dual_graph ts rct in
	let ordering = build_perfect_elemination_ordering graph in
	do_greedy_coloring graph ordering
;;

(* Computes the CSR "row ptr", "col ind" as well as an array TT[ti=0,1,...,nt-1][vi=0,1,2][vj=0,1,2] = index in "col ind" of A(T[ti][vi],T[ti][vj]) 
	Once again BigArrays as it's meant to be pushed to the GPU *)
let build_crs_profile ts rct =
	
	(* Assembles a single row of the linear system *)
	let do_row ~init ~f row =
		rct.(row) |> List.fold ~init ~f:(fun acc ti ->
			let col0, col1, col2 = Int32.(to_int_exn ts.{ti,0}, to_int_exn ts.{ti,1}, to_int_exn ts.{ti,2}) in
			let trow = if col0 = row then 0 else if col1 = row then 1 else 2 in
			let do_col = f ti trow row 
			in acc |> do_col 0 col0 |> do_col 1 col1 |> do_col 2 col2
		)	
	in 

	(* Simulates the assembly in order to count the number of non zero coefficients *)
	let nnz = 
		(* color helps us to track already counted coeffs *)
		let color = Array.create ~len:(Array.length rct) (-1) in
		(* "fold" all the rows *)
		rct |> Array.foldi ~init:0 ~f:(fun row cnt _ ->
			row |> do_row ~init:cnt ~f:(fun ti trow row tcol col cnt ->
				(* Do not count coeffs twice *)
				if color.(col) < row then (
					color.(col) <- row;
					succ cnt
				) else
					cnt
			)
		)

	(* Perform the real assembly this time *)
	and nv = Array.length rct in		

	(* color helps us to track already counted coeffs *)
	let color = Array.create ~len:nv (-1)

	(* edges helps us to remember triangles linked to an edge *)
	and edges = Array.create ~len:nnz [|(0,0); (0,0)|]

	(* The result arrays, colidx will later be converted to a BigArray as BigArray does not provide an in-place sort function *)
	and rowptr = Array1.create int32 c_layout (nv + 1)
	and colidx = Array.create ~len:nnz Int32.zero
	and tt = Array2.create int32 c_layout (Array2.dim1 ts) 9 in

	(* This sorts the nz column idx of a row before updating rowptr, colidx and tt *)
	let commit_row rstart rend = 		
		Array.sort ~pos:rstart ~len:(rend - rstart) ~cmp:compare colidx;
		for pos = rstart to (rend-1) do
			let col = Int32.to_int_exn colidx.(pos) 
			and pos' = Int32.of_int_exn pos in
			let t0, idx0 = edges.(col).(0)
			and t1, idx1 = edges.(col).(1) in
			tt.{t0, idx0} <- pos';
			tt.{t1, idx1} <- pos'
		done
	in

	(* "fold" all the rows once again *)
 	let last = rct |> Array.foldi ~init:0 ~f:(fun row pos _ ->		
		let pos' = row |> do_row ~init:pos ~f:(fun ti trow row tcol col pos ->
			let idx = 3 * trow + tcol in
			if color.(col) < row then (
				(* We found the first edge going from row to col, 
					this could be the only one if it's a boundary edge *)
				color.(col) <- row;
				edges.(col).(0) <- (ti, idx);
				edges.(col).(1) <- (ti, idx);
				colidx.(pos) <- Int32.of_int_exn col;
				succ pos
			) else (
				(* there is a second edge going from row to col *)
				edges.(col).(1) <- (ti, idx);
				pos
			)
		) 
		in

		(* columns idx of the current row starts at position pos and ends at position pos' - 1 *)
		rowptr.{row} <- Int32.of_int_exn pos;
		commit_row pos pos';
		pos'
	)
	in rowptr.{nv} <- Int32.of_int_exn last;

	(rowptr, (Array1.of_array int32 c_layout colidx), tt)
;;

let do_global_assembly_on_gpu ~vs ~ts ~colors ~rowptr ~colidx ~tt = 
	let aux = foreign "do_global_assembly_on_gpu" (
		Ctypes.size_t @-> (* nv*)
    	Ctypes.size_t @-> (* nt *)
    	ptr Ctypes.double @-> (* vs *)
    	ptr Ctypes.int32_t @-> (* ts*)
    	ptr Ctypes.int32_t @-> (* color *)
    	ptr Ctypes.int32_t @-> (* rowptr *)
    	ptr Ctypes.int32_t @-> (* colidx *)
    	ptr Ctypes.int32_t @-> (* tt *)
    	ptr Ctypes.double @-> (* output coeffs*)
    	returning void
    )
    and addr1 ba = bigarray_start array1 ba
    and addr2 ba = bigarray_start array2 ba
    and coeffs = Array1.create float64 c_layout (Array2.dim1 vs) in 
    let open Unsigned.Size_t in
	aux	(of_int (Array2.dim1 vs)) 
    	(of_int (Array2.dim1 ts)) 
    	(addr2 vs) 
    	(addr2 ts) 
    	(addr1 colors) 
    	(addr1 rowptr) 
    	(addr1 colidx) 
    	(addr2 tt)
    	(addr1 coeffs);
    coeffs
;;

let vs, ts = read_mesh_exn "lyttelton.mesh" in
let rct = build_reverse_connectivity_table ts (Array2.dim1 vs) in
let colors = build_perfect_coloring ts rct in
(* let () = for i = 0 to (Array1.dim colors) -1 do Printf.printf "%d\n" colors.{i} done in *)
let rowptr, colidx, tt = build_crs_profile ts rct in
do_global_assembly_on_gpu ~vs ~ts ~colors ~rowptr ~colidx ~tt
