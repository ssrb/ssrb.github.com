module Int =
	struct type t = int
	let compare = compare
end

let (--) i j =
	let rec aux n acc= if n < i then acc else aux (n-1) (n::acc)
	in aux (j-1) []

module IntMap = Map.Make(Int)
module IntSet = Set.Make(Int)

let scan_domains chan =
	let rec aux acc =
		try
			aux (Scanf.fscanf chan "%d\n" (fun d -> d::acc))
		with
		| End_of_file -> List.rev acc
	in aux []

let scan_triangles chan domains =

	let check_iface v2d did v iface =
		try
			let did' = IntMap.find v v2d in
			if did' != did then IntSet.add v iface else iface
		with
		| Not_found -> iface
	in

	let receive_one_triangle (v2d, iface, dids) v1 v2 v3 =
		match dids with
		| did::dids' ->
			let v2d' = v2d |> IntMap.add v1 did |> IntMap.add v2 did |> IntMap.add v3 did in
			let check_iface = check_iface v2d did in
			let iface = iface |> check_iface v1 |> check_iface v2 |> check_iface v3 in
			(v2d', iface, dids')
		| [] -> (v2d, iface, [])
	in 

	let rec aux acc = 
		try
			aux (Scanf.fscanf chan "%d %d %d\n" (receive_one_triangle acc))
		with
		| End_of_file -> acc
	in
	let _ = Scanf.fscanf chan "%d\n" (fun n -> n) in
	let v2d, iface, _ = aux (IntMap.empty, IntSet.empty, domains) in (v2d, iface)

let () =
	let tchan = open_in Sys.argv.(1)
	and dchan = open_in Sys.argv.(2) in
	let domains = scan_domains dchan in
	let ndomains = 1 + List.fold_left max 0 domains in
	let v2d, iface = scan_triangles tchan domains in
	let interior_chans = Array.init ndomains (fun ic -> open_out ("domain." ^ (string_of_int ic) ^ ".vids"))
	and inter_chan = open_out "interface.vids" in
	let chan vid did = if IntSet.mem vid iface then inter_chan else interior_chans.(did) in
	begin
		IntMap.iter (fun vid did -> Printf.fprintf (chan vid did) "%d\n" vid) v2d;
		0--ndomains |> List.iter (fun ic -> close_out interior_chans.(ic))
	end
