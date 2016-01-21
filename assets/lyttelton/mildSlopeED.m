function [] = mildSlopeED()

	profile on;
	
	% Google mercator
	k0 = 1 / (3 * 256 / (2 * pi * 6378137 * 0.00955042966330666));

	% Parameters
	wavedirection = [1 -1];
	wavedirection /= norm(wavedirection);
	waveheight = 2;

	R = 0.9;
	g = 9.81;
	cf = 0.5;

	[vertices, triangles, depth, openBoundary, closedBoundary] = readMesh("lyttelton.mesh");
	nvertices = length(vertices);
	ntriangles = length(triangles);

	% Precompute a few things
	n0s = n0(depth, 1:ntriangles, k0);
	hs = h(depth, 1:ntriangles);

	[openBoundaryElements, closedBoundaryElements] = createEdgeToElementMappings(vertices, triangles, openBoundary, closedBoundary);

	q1 = vertices(triangles(:, 1), :);
	q2 = vertices(triangles(:, 2), :);
	q3 = vertices(triangles(:, 3), :);

	u = q2 - q3;
	v = q3 - q1;
	w = q1 - q2;

	a(:,1,:) = v';
	a(:,2,:) = w';
	areas = 0.5 * cellfun(@det, num2cell(a,[1,2]))(:);

	% Picard iteration
	eta = zeros(1, nvertices);
	for i=1:20		
		[A, b] = assemble(	vertices, triangles, 
							openBoundary, closedBoundary, 
							openBoundaryElements, closedBoundaryElements, 
							k0, wavedirection, waveheight, 
							R, cf, 
							u, v, w, areas, 
							n0s, hs,
							eta);
		eta = A \ b;
	end
	
	% Reorder the global solution
 	writeSolution("lyttelton.sol", eta);
 	dlmwrite("lytteltonsol.txt", eta);
 	profile off;
 	prof = profile ("info");
	profshow (prof);

endfunction

function [A, b] = assemble(	vertices, triangles, 
							openBoundary, closedBoundary, 
							openBoundaryElements, closedBoundaryElements, 
							k0, wavedirection, waveheight, 
							R, cf, 
							u, v, w, areas, 
							n0s, hs, 
							eta0)

	nvertices = length(vertices);

	Wf = (8 * cf * mean(abs(eta0(triangles)), 2)) ./ (3 * pi * sinh(k0 * hs).^3);

	% The modified wave number
	ps = i * k0 * sqrt(1 - (i * Wf) ./ n0s); % +/-

	% 4.2.1 Second integral
	stiffness = - n0s ./ (4 * k0 * k0 * areas); 

	val(:, 1) = stiffness .* sum(u.*u, 2);
	val(:, 2) = stiffness .* sum(u.*v, 2);
	val(:, 3) = stiffness .* sum(u.*w, 2);
	val(:, 5) = stiffness .* sum(v.*v, 2);
	val(:, 6) = stiffness .* sum(v.*w, 2);
	val(:, 9) = stiffness .* sum(w.*w, 2);
	val(:, [ 4 , 7 , 8]) = val(:, [ 2 , 3 , 6 ]);

	% 4.2.1 First integral
	massDiag = (n0s - i * Wf) .* areas / 6;
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

		n = n0s(closedBoundaryElements(edgei));
		p = ps(closedBoundaryElements(edgei));

		len = norm(vertices(ii,:) - vertices(jj,:));

		A(ii, jj) += closedPrefactor * n * (p * len / 6 + 1 / (2 * p * len));
		A(jj, ii) += closedPrefactor * n * (p * len / 6 + 1 / (2 * p * len));
		A(ii, ii) += closedPrefactor * n * (p * len / 3 - 1 / (2 * p * len));
		A(jj, jj) += closedPrefactor * n * (p * len / 3 - 1 / (2 * p * len));
	end

	openPrefactor = -i / (k0 * k0);
	for edgei=1:length(openBoundary)
		edge = openBoundary(edgei, :);

		ii = edge(1);
		jj = edge(2);

		vi = vertices(ii, :);
		vj = vertices(jj, :);

		len = norm(vi - vj);

		n = n0s(openBoundaryElements(edgei));
		p = ps(openBoundaryElements(edgei));

		A(ii, jj) += openPrefactor * n * (p * len / 6 + 1 / (2 * p * len));
		A(jj, ii) += openPrefactor * n * (p * len / 6 + 1 / (2 * p * len));
		A(ii, ii) += openPrefactor * n * (p * len / 3 - 1 / (2 * p * len));
		A(jj, jj) += openPrefactor * n * (p * len / 3 - 1 / (2 * p * len));

		wavei = IncomingWaveElevation(vi, k0, wavedirection, waveheight);
		wavej = IncomingWaveElevation(vj, k0, wavedirection, waveheight);

		b(ii) += -openPrefactor * n * (wavej - wavei) / (2 * p * len);
		b(jj) += -openPrefactor * n * (wavei - wavej) / (2 * p * len);

		x = vj - vi;
		x = [x(2) -x(1)];
		x /= norm(x);

		b(ii) += -openPrefactor * n * p * (x * wavedirection' - 1) * (len / 6) * (2 * wavei + wavej);
		b(jj) += -openPrefactor * n * p * (x * wavedirection' - 1) * (len / 6) * (wavei + 2 * wavej);
	end

endfunction

function [elevation] = IncomingWaveElevation(pos, k0, wavedirection, waveheight)
	elevation = waveheight * exp(-i *k0 * pos * wavedirection');
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

function [vertices, triangles, depth, openBoundary, closedBoundary] = readMesh(fileName)

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
		triangles = [triangles; v1 v2 v3];
		depth = [depth; d];
	end
	fseek (fid, 0, SEEK_SET);

	% Vertices
	fgoto(fid, "Vertices");
	[nvertices] = fscanf(fid, "%d", "C");
	fgets(fid);

	vertices = zeros(length(nvertices), 2);

	for vid=1:nvertices
		[vertices(vid, 1) vertices(vid, 2)] = fscanf(fid, "%f %f %d\n", "C");
	end

	% Edges
	fseek (fid, 0, SEEK_SET);

	openBoundary = [];
	closedBoundary = [];
	
	fgoto(fid, "Edges");
	[nedges] = fscanf(fid, "%d", "C");
	for edgei=1:nedges
		[v1, v2, boundary] = fscanf(fid, "%d %d %d\n", "C");
		if boundary == 2
			openBoundary = [openBoundary; v1 v2];
		elseif boundary == 1
			closedBoundary = [closedBoundary; v1 v2];
		end
	end

	fclose(fid);
endfunction

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
	d = depth(element) * 256 / (2 * pi * 6378137 * 0.00955042966330666);
endfunction