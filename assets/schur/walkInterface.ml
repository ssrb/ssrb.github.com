open Scanf

module Int = 
	struct type t = int 
	let compare = compare 
end

module IntMap = Map.Make(Int)
module IntSet = Set.Make(Int)

let (--) i j = 
	let rec aux n acc = if n < i then acc else aux (n-1) (n::acc)
	in aux j []

let multi_add k v mmap =
	try 
		let vals = IntMap.find k mmap in
		IntMap.add k (IntSet.add v vals) mmap
	with
	| Not_found -> IntMap.add k (IntSet.singleton v) mmap

(* Read the interface file *)
let on_interface = 
	let interfaceFile = open_in Sys.argv.(1) in
	let rec readInterface vint = 
		let receiver v = vint |> IntSet.add v in
		try
			let vint' = Scanf.fscanf interfaceFile "%d\n" receiver in
			readInterface vint'
		with
		| _ -> vint
	in 
	let vint = readInterface IntSet.empty in
 	fun v -> IntSet.mem v vint

(* Build interface vertices <=> interface triangles mappings*)
let get_domain_interface did =	
	let metisGraph = open_in Sys.argv.(2)
	and domains = open_in Sys.argv.(3)
	(* First line is the number of triangles *) 	
	and add_if_on_interface v ti v2t = 
		if on_interface v then 
			multi_add v ti v2t
		else v2t 
	in 
	let receiver (v2t, t2v) ti did' v1 v2 v3 = 
		if did == did' then
			(v2t |> add_if_on_interface v1 ti |> add_if_on_interface v2 ti |> add_if_on_interface v3 ti , 
			if on_interface v1 || on_interface v2 || on_interface v3 then 
				IntMap.add ti (v1, v2, v3) t2v 
			else 
				t2v)
		else 
			(v2t, t2v)
	and nbTriangles = Scanf.fscanf metisGraph "%d\n" (fun x -> x) in
	let rec readTriangles maps ti =
		if ti == nbTriangles then 
			maps 
		else
			((* 1 triangle (3 vertices) per line*)
			let did' = Scanf.fscanf domains "%d\n" (fun x -> x) in
			let maps' = Scanf.fscanf metisGraph "%d %d %d\n" (receiver maps ti did') in
			readTriangles maps' (ti+1))
	in readTriangles (IntMap.empty, IntMap.empty) 0

let graph = 
	let build_domain_interface g did =
		let edges_of_triangle (v1, v2, v3) = [(v1, v2); (v2, v3); (v3, v1)] in
		let is_boundary_edge v2t (v1,v2) =
			(if on_interface v1 && on_interface v2 then
				let ts = IntMap.find v1 v2t
				and ts' = IntMap.find v2 v2t	
				in 1 == IntSet.cardinal (IntSet.inter ts ts')
			else
				false)
		in		
		let v2t,t2v = get_domain_interface did in
		let l = IntMap.fold (fun _ v -> (@) (edges_of_triangle v)) t2v [] in
		let l' = List.filter (is_boundary_edge v2t) l in
		let add_edge_to_graph g (v1,v2)= g |> multi_add v1 v2 |> multi_add v2 v1 in
		List.fold_left add_edge_to_graph g l'
	in 
	1--2 |> List.fold_left build_domain_interface IntMap.empty

type 'a state = 
{
	pred : int IntMap.t;
	acc : 'a
}

let walk_interface g = 
	let dfs g c fn state = 
		let rec dfs_visit t u {pred; acc} =
			let edge v (t, state) =
				if IntMap.mem v state.pred then
					(t, state)
				else dfs_visit t v {state with pred = IntMap.add v u state.pred;}
		 	in
		 	let t = t + 1 in
	 		IntSet.fold edge (IntMap.find u g) (t, {pred; acc = fn acc u})
	 	in
	 	dfs_visit 0 c state
	in
	let starts = IntMap.filter (fun _ v -> 1 == IntSet.cardinal v) g in
	let find_start pred =
		let starts' = IntMap.filter (fun k _ -> not (IntMap.mem k pred)) starts in
		if IntMap.is_empty starts' then
			None
		else
			Some(fst (IntMap.choose starts'))
	in
	let rec aux (starts, pred) g =
		match find_start pred with
		| None -> (starts, pred)
		| Some(s) -> (let _, state = (dfs g s (fun _ u -> Printf.printf "%d\n" u) {pred = pred; acc = ()})
						in 	aux (s::starts, state.pred) g)
	in aux ([], IntMap.empty) g;;

 walk_interface graph
