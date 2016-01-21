(* Added by OPAM. *)
(*let () =
  try Topdirs.dir_directory (Sys.getenv "OCAML_TOPLEVEL_PATH")
  with Not_found -> ()
;;

#use "topfind";;
#thread;;
#camlp4o;;
#require "core.top";;
#require "core.syntax";;*)

open Core.Std

(* Blah blah blah *)
let vs,ts = 
	let open In_channel in
	with_file "lyttelton.mesh" ~f:
		(fun chan ->
			let module Parser = struct
				type t = V1 | V2 | V3 of int | T1 | T2 | T3 of int
				let fold call (state, vs, ts) line =
					match state with
					| V1 -> (if line = "Vertices" then V2 else V1), vs, ts
					| V2 -> let nvertices = Scanf.sscanf line "%d" (fun n -> n) in (V3 nvertices), vs, ts
					| V3 n when n > 0 -> let v = Scanf.sscanf line "%f %f %f" (fun x y _ -> (x,y)) in (V3 (pred n)), v::vs, ts
					| V3 _ -> (T1, vs, ts)
					| T1 ->  (if line = "Triangles" then T2 else T1), vs, ts
					| T2 -> let ntriangles = Scanf.sscanf line "%d" (fun n -> n) in	(T3 ntriangles, vs, ts)
					| T3 n when n > 0 -> let t = Scanf.sscanf line "%d %d %d %d" (fun v1 v2 v3 _ -> (pred v1, pred v2, pred v3)) in (T3 (pred n)), vs, t::ts
					| T3 _ -> call.return(state, vs, ts)
				let init = (V1, [],[])
			end
			in let _, vs, ts = with_return (fun call ->
				chan |> fold_lines ~init:Parser.init ~f:(Parser.fold call))
			in List.rev vs, List.rev ts)
;;

(* Blah blah blah *)
let build_reverse_connectivity_table ts nv = 
	let rct = Array.create ~len:nv [] in
	List.iteri ts ~f:(fun t (v1,v2,v3) ->	
		rct.(v1) <- t::rct.(v1); 
		rct.(v2) <- t::rct.(v2); 
		rct.(v3) <- t::rct.(v3));
	rct
;;

(* Blah blah blah *)
let build_triangulation_dual_graph ts nv = 	
	let graph = Array.create ~len:(List.length ts) [] in
	let rct = build_reverse_connectivity_table ts nv in
	List.iteri ts ~f:(fun t (v1,v2,v3) ->
		graph.(t) <- List.dedup (rct.(v1) @ rct.(v2) @ rct.(v3)));
	graph
;;

(* Blah blah blah *)
let build_perfect_elemination_ordering graph =	
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
;;

(* Blah blah blah *)
let do_greedy_coloring graph ordering =
	let n = Array.length graph in 
	let colors = Array.create ~len:n n in
	let perfect_color v =
		let rec aux color = function
			| c::cs' -> if color < c then color else aux (succ color) cs'
			| [] -> color
		in graph.(v) |> List.map ~f:(fun v -> colors.(v)) |> List.sort ~cmp:compare |> aux 0
	in
	List.iter ordering ~f:(fun v -> 	
		colors.(v) <- perfect_color v
	);
	colors
;;

(* Blah blah blah *)
let graph = build_triangulation_dual_graph ts (List.length vs)
let ordering = build_perfect_elemination_ordering graph
let perfect_coloring = do_greedy_coloring graph ordering
let () = Array.iter ~f:(Printf.printf "%d\n") perfect_coloring

