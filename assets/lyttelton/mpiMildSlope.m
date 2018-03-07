function mpiMildSlope()

	pkg load mpitb;

	MPI_Init;

	[mpistatus, mpirank] = MPI_Comm_rank(MPI_COMM_WORLD);
	[mpistatus, mpicommsize] = MPI_Comm_size(MPI_COMM_WORLD);

	if mpirank == 0
		profile on;
	end

	domain = mpirank + 1;

	[vertices, triangles, depth, openBoundary, closedBoundary, interface, vids] = readMesh("lyttelton.mesh", "interface2.vids", domain);

	k0 = 1 / (3 * 256 / (2 * pi * 6378137 * 0.00955042966330666));
	wavedirection = [1 -1];
	waveheight = 2;
	R = 0.9;

	% Assemble *local* stiffness matrix and load vector
	[A, b] = assemble(vertices, triangles, depth, openBoundary, closedBoundary, k0, wavedirection, waveheight, R);

	nbInterior = length(b) - length(interface);

	% LU-factorize the interior of the subdomains, we're going to reuse this everywhere 
	[L, U, p, tmp] = lu(A(1:nbInterior, 1:nbInterior), 'vector');
	q(tmp) = 1:length(tmp);

	% We solve for the trace first: we need to compute the second member of the Schur complement system
	% Local contribution
	bti = computeBTildI(A, b, nbInterior, L, U, p, q);

	% Sum contributions
	bt = AllReduceComplex(bti, MPI_SUM, MPI_COMM_WORLD);

	% IncomingWaveElevation(vertices(nbInterior+1:end, :), k0, wavedirection, waveheight)	
	% Solve for the trace doing distributed gradient descent
	trcSol = gmres(	@(x) parallelMultiplyBySchurComplement(A, nbInterior, L, U, p, q, x), ...
					bt);
	
	% Solve for the local interior
	localISol = (U \ (L \ (b(1:nbInterior) - A(1:nbInterior, nbInterior + 1 : end) * trcSol)(p)))(q);
	%localISol = real(IncomingWaveElevation(vertices(1:nbInterior, :), k0, wavedirection, waveheight));

	% Consolidate solutions
	% Gather sizes of the subdomains
	allNbInterior = zeros(mpicommsize, 1);
	MPI_Gather(nbInterior, allNbInterior, 0, MPI_COMM_WORLD);

	sumNbInterior = sum(allNbInterior);

	disps = cumsum([0; allNbInterior(1:end-1)]);

	% Concatenate local => global mappings of the unknowns
	globalISol = GathervComplex(localISol, allNbInterior, disps, 0 , MPI_COMM_WORLD);

	% Concatenate solutions
	allVids = zeros(sumNbInterior, 1);
	vids = vids(1:nbInterior);
	MPI_Gatherv(vids, allVids, allNbInterior, disps, 0, MPI_COMM_WORLD);

	if mpirank == 0
		% Reorder the global solution
	 	solution = zeros(sumNbInterior + length(interface), 1);
	 	solution(allVids) = globalISol;
	 	solution(interface) = trcSol;

	 	dlmwrite("lytteltonsol.txt", solution);

	 	solution = real(solution);

	 	writeSolution("lyttelton.sol", solution);
	end

	MPI_Finalize;

	if mpirank == 0
		profile off;
		prof = profile ("info");
		profshow (prof);
	end

endfunction

function [A, b] = assemble(vertices, triangles, depth, openBoundary, closedBoundary, k0, wavedirection, waveheight, R)

	[openBoundaryElements, closedBoundaryElements] = createEdgeToElementMappings(vertices, triangles, openBoundary, closedBoundary);

	wavedirection /= norm(wavedirection);

	nvertices = length(vertices);
	ntriangles = length(triangles);

	q1 = vertices(triangles(:, 1), :);
	q2 = vertices(triangles(:, 2), :);
	q3 = vertices(triangles(:, 3), :);

	u = q2 - q3;
	v = q3 - q1;
	w = q1 - q2;

	a(:,1,:) = v';
	a(:,2,:) = w';
	areas = 0.5 * cellfun(@det, num2cell(a,[1,2]))(:);

	% 4.2.1 Second integral
	n0s = n0(depth, 1:ntriangles, k0);
	stiffness =  - n0s ./ (4 * k0 * k0 * areas); 

	val(:, 1) = stiffness .* sum(u.*u, 2);
	val(:, 2) = stiffness .* sum(u.*v, 2);
	val(:, 3) = stiffness .* sum(u.*w, 2);
	val(:, 5) = stiffness .* sum(v.*v, 2);
	val(:, 6) = stiffness .* sum(v.*w, 2);
	val(:, 9) = stiffness .* sum(w.*w, 2);
	val(:, [ 4 , 7 , 8]) = val(:, [ 2 , 3 , 6 ]);

	% 4.2.1 First integral
	massDiag = n0s .* areas / 6;
	mass = massDiag / 2;

	val(:, 1) += massDiag;
	val(:, 2) += mass;
	val(:, 3) += mass;
	val(:, 4) += mass;
	val(:, 5) += massDiag;
	val(:, 6) += mass;
	val(:, 7) += mass;
	val(:, 8) += mass;
	val(:, 9) += massDiag;
	
	is = triangles(:, [1 1 1 2 2 2 3 3 3]);
	js = triangles(:, [1 2 3 1 2 3 1 2 3]);

	A = sparse(is(:), js(:), val(:), nvertices, nvertices);
	b = zeros(nvertices, 1);

	closedPrefactor = -i * (1 - R) / (k0 * k0 * (1 + R));
	for edgei=1:length(closedBoundary)
		edge = closedBoundary(edgei, :);

		ii = edge(1);
		jj = edge(2);

		n =  n0s(closedBoundaryElements(edgei));
		
		len = norm(vertices(ii,:) - vertices(jj,:));

		A(ii, jj) += closedPrefactor * n * (k0 * len / 6 + 1 / (2 * k0 * len));
		A(jj, ii) += closedPrefactor * n * (k0 * len / 6 + 1 / (2 * k0 * len));
		A(ii, ii) += closedPrefactor * n * (k0 * len / 3 - 1 / (2 * k0 * len));
		A(jj, jj) += closedPrefactor * n * (k0 * len / 3 - 1 / (2 * k0 * len));
	end

	openPrefactor = -i / (k0 * k0);
	for edgei=1:length(openBoundary)
		edge = openBoundary(edgei, :);

		ii = edge(1);
		jj = edge(2);

		vi = vertices(ii, :);
		vj = vertices(jj, :);

		len = norm(vi - vj);

		n =  n0s(openBoundaryElements(edgei));

		A(ii, jj) += openPrefactor * n * (k0 * len / 6 + 1 / (2 * k0 * len));
		A(jj, ii) += openPrefactor * n * (k0 * len / 6 + 1 / (2 * k0 * len));
		A(ii, ii) += openPrefactor * n * (k0 * len / 3 - 1 / (2 * k0 * len));
		A(jj, jj) += openPrefactor * n * (k0 * len / 3 - 1 / (2 * k0 * len));

		wavei = IncomingWaveElevation(vi, k0, wavedirection, waveheight);
		wavej = IncomingWaveElevation(vj, k0, wavedirection, waveheight);

		b(ii) += -openPrefactor * n * (wavej - wavei) / (2 * k0 * len);
		b(jj) += -openPrefactor * n * (wavei - wavej) / (2 * k0 * len);

		x = vj - vi;
		x = [x(2) -x(1)];
		x /= norm(x);

		b(ii) += -openPrefactor * n * k0 * (x * wavedirection' - 1) * (len / 6) * (2 * wavei + wavej);
		b(jj) += -openPrefactor * n * k0 * (x * wavedirection' - 1) * (len / 6) * (wavei + 2 * wavej);
	end

endfunction

function [elevation] = IncomingWaveElevation(pos, k0, wavedirection, waveheight)
	elevation = waveheight * exp(-i * k0 * pos * wavedirection');
endfunction

function [bti] = computeBTildI(A, b, nbInterior, L, U, p, q)
	bti = b(nbInterior + 1 : end) ...
		- A(nbInterior + 1 : end, 1:nbInterior) * (U \ (L \ b(1:nbInterior)(p)))(q);
endfunction

function [res] = parallelMultiplyBySchurComplement(A, nbInterior, L, U, p, q, x)
	local = A(nbInterior + 1 : end, nbInterior + 1 : end) * x ...
		- A(nbInterior + 1 : end, 1:nbInterior) * (U \ (L \ (A(1:nbInterior, nbInterior + 1 : end) * x)(p)))(q);
	res = AllReduceComplex(local, MPI_SUM, MPI_COMM_WORLD);
endfunction

function [openBoundaryElements, closedBoundaryElements] = createEdgeToElementMappings(vertices, triangles, openBoundary, closedBoundary)

	rct = cell(length(vertices), 1);
	for ti=1:length(triangles)
		for si=1:3
			v = triangles(ti, si);
			rct{v} = [rct{v} ti];
		end
	end

	openBoundaryElements = [];
	closedBoundaryElements = [];

	if length(openBoundary)
		openBoundaryElements = cellfun(@intersect, rct(openBoundary(:, 1)), rct(openBoundary(:, 2)));
	end

	if length(closedBoundary)
		closedBoundaryElements = cellfun(@intersect, rct(closedBoundary(:, 1)), rct(closedBoundary(:, 2)));
	end

endfunction

function [vertices, triangles, depth, openBoundary, closedBoundary, interface, vids] = readMesh(fileName, interfaceName, domainIdx)

	interface = readInterface(interfaceName);

	fid = fopen (fileName, "r");
	depthfid = fopen ("depth/out", "r");

	% Elements
	fgoto(fid, "Triangles");
	[ntriangles] = fscanf(fid, "%d", "C");
	triangles = [];
	depth = [];
	for tid=1:ntriangles
		[v1, v2, v3, domain] = fscanf(fid, "%d %d %d %d\n", "C");
		d = fscanf(depthfid, "%f\n", "C");
		if domain == domainIdx
			triangles = [triangles; v1 v2 v3];
			depth = [depth; d];
		end
	end
	fseek (fid, 0, SEEK_SET);

	sortedVertexIds = unique(triangles(:));

	% Put the interface unknowns at the end, using the numbering we computed in the last post
	vids = [setdiff(sortedVertexIds, interface); interface];

	globalToLocal(vids) = 1:length(vids);
	triangles(:) = globalToLocal(triangles(:));

	% Vertices
	fgoto(fid, "Vertices");
	[nvertices] = fscanf(fid, "%d", "C");
	fgets(fid);

	vertices = zeros(length(vids), 2);

	currentLine = 1;
	for vid=1:length(sortedVertexIds)
		globalId = sortedVertexIds(vid);
		fskipLines(fid, globalId - currentLine);
		currentLine = globalId;
		localId = globalToLocal(globalId);
		[vertices(localId, 1) vertices(localId, 2)] = fscanf(fid, "%f %f %d\n", "C");
		++currentLine;
	end

	% Edges
	fseek (fid, 0, SEEK_SET);

	openBoundary = [];
	closedBoundary = [];
	
	fgoto(fid, "Edges");
	[nedges] = fscanf(fid, "%d", "C");
	for edgei=1:nedges
		[v1, v2, boundary] = fscanf(fid, "%d %d %d\n", "C");
		lv1 = globalToLocal(v1);
		lv2 = globalToLocal(v2);
		if lv1 && lv2
			if boundary == 2
				openBoundary = [openBoundary; lv1 lv2];
			elseif boundary == 1
				closedBoundary = [closedBoundary; lv1 lv2];
			end
		end
	end

	fclose(fid);
endfunction

function [interface] = readInterface(fileName)
	fid = fopen (fileName, "r");
	interface = [];
	while (l = fgetl(fid)) != -1
		[vid] = sscanf(l, "%d\n", "C");
		interface = [interface; vid];
	end
	fclose(fid);
end

function [] = fgoto(fid, tag)
	while !strcmp(fgetl(fid),tag)
	end
endfunction

function [] = fskipLines(fid, nlines)
	for i=1:nlines
		fgets(fid);
	end
endfunction

function [] = writeSolution(fileName, solution)
	fid = fopen (fileName, "w");
	fprintf(fid, "MeshVersionFormatted 1\n\nDimension 2\n\nSolAtVertices\n%d\n1 1\n\n", length(solution));
	for i=1:length(solution)
		fprintf(fid, "%e\n", solution(i));
	end
	fclose(fid);
endfunction

function [value] = n0(depth, element, k0)
	kh = k0 * h(depth, element);
	tanhkh = tanh(kh);
	value = 0.5 * (1 + kh .* (1 - tanhkh .* tanhkh) ./ tanhkh);
endfunction

function [d] = h(depth, element)
	% port depth is maintained around 10.5m
	d = depth(element) * 256 / (2 * pi * 6378137 * 0.00955042966330666);
endfunction

function [res] = AllReduceComplex(sendbuf, op, comm)
	tmp = real(sendbuf);
	tmpres = zeros(size(tmp));
	mpistatus = MPI_Allreduce(tmp, tmpres, op, comm);
	res = tmpres;
	tmp = imag(sendbuf);
	mpistatus = MPI_Allreduce(tmp, tmpres, op, comm);
	res += i * tmpres;
endfunction

function [res] = GathervComplex(sendbuf, recvcounts, disps, root, comm);
	tmp = real(sendbuf);
	tmpres = zeros(sum(recvcounts), 1);
	MPI_Gatherv(tmp, tmpres, recvcounts, disps, root, comm);
	res = tmpres;
	tmp = imag(sendbuf);
	MPI_Gatherv(tmp, tmpres, recvcounts, disps, root , comm);
	res += i * tmpres;
endfunction