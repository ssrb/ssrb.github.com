function schurComplementMethod2()

	[vertices, triangles, border, domains] = readMesh("poisson2D.mesh");
	interface = readInterface("interface2.vids");

	% Assemble local Dirichlet problems 
	% => could be done in parallel
	[A1, b1, vids1] = subAssemble(vertices, triangles, border, domains, interface, 1);
	[A2, b2, vids2] = subAssemble(vertices, triangles, border, domains, interface, 2);

	nbInterior1 = length(b1) - length(interface);
	nbInterior2 = length(b2) - length(interface);

	% LU-factorize the interior of the subdomains, we're going to reuse this everywhere 
	% => could be done in parallel
	[L1, U1, p1, tmp] = lu(A1(1:nbInterior1, 1:nbInterior1), 'vector');
	q1(tmp) = 1:length(tmp);
	[L2, U2, p2, tmp] = lu(A2(1:nbInterior2, 1:nbInterior2), 'vector');
	q2(tmp) = 1:length(tmp);

	% In order to solve for the flux first, we need to compute the corresponding second member
	% => could be done in parallel
	bt1 = computeBTildI(A1, b1, nbInterior1, L1, U1, p1, q1);
	bt2 = computeBTildI(A2, b2, nbInterior2, L2, U2, p2, q2);

	epsilon = 1.e-30;
	maxIter = 600

	% Each pcg contributing to the second member could be done in parallel
	btt = pcg(@(x) multiplyByLocalSchurComplement(A2, nbInterior2, L2, U2, p2, q2, x), bt2, epsilon, maxIter) ... 
			- pcg(@(x) multiplyByLocalSchurComplement(A1, nbInterior1, L1, U1, p1, q1, x), bt1, epsilon, maxIter);

	% Solve for the flux => each nested pcg contribution could be done in parallel
	flux = pcg(@(x) pcg(@(y) multiplyByLocalSchurComplement(A1, nbInterior1, L1, U1, p1, q1, y), x, epsilon, maxIter) ...
				+ pcg(@(y) multiplyByLocalSchurComplement(A2, nbInterior2, L2, U2, p2, q2, y), x, epsilon, maxIter), btt, epsilon, maxIter);

	% Add the computed Neumann data to the load vectors
	b1(nbInterior1 + 1 : end) += flux;
	b2(nbInterior2 + 1 : end) -= flux;
	% and solve the corresponding problems using a direct method
	% => could be done in parallel
	solution = zeros(length(vertices), 1);
	solution(vids1) = A1 \ b1;
	solution(vids2) = A2 \ b2;

	writeSolution("poisson2D.sol", solution);

endfunction

function [A, b, dVertexIds] = subAssemble(vertices, triangles, border, domains, interface, domainIdx)
	% Keep the triangles belonging to the subdomain
	dTriangles = triangles(find(domains == domainIdx),:);
	dVertexIds = unique(dTriangles(:));

	% Put the interface unknowns at the end, using the numbering we computed in the last post
	dVertexIds = [setdiff(dVertexIds, interface); interface];

	% Switch domain triangles to local vertices numbering
	% Should use a hashtable here but their is none in octave ...
	globalToLocal = zeros(length(vertices), 1);
	globalToLocal(dVertexIds) = 1:length(dVertexIds);
	dTriangles(:) = globalToLocal(dTriangles(:));

	% Assemble the Dirichlet problem
	[A, b] = assemble(vertices(dVertexIds, :), dTriangles, border(dVertexIds));
endfunction

function [A, b] = assemble(vertices, triangles, border)
	nvertices = length(vertices);
	ntriangles = length(triangles);

	iis = [];
	jjs = [];
	vs = [];

	b = zeros(nvertices, 1);
	for tid=1:ntriangles
		q = zeros(3,2);
		q(1,:) = vertices(triangles(tid, 1), :) - vertices(triangles(tid, 2), :);
		q(2,:) = vertices(triangles(tid, 2), :) - vertices(triangles(tid, 3), :);
		q(3,:) = vertices(triangles(tid, 3), :) - vertices(triangles(tid, 1), :);
		area = 0.5 * det(q([1,2], :));
		for i=1:3
			ii = triangles(tid,i);
			if !border(ii)
				for j=1:3
					jj = triangles(tid,j);
					if !border(jj)
						hi = q(mod(i, 3) + 1, :);
						hj = q(mod(j, 3) + 1, :);
						
						v = (hi * hj') / (4 * area);

						iis = [iis ii];
						jjs = [jjs jj];
						vs = [vs v];
					end
				end
				b(ii) += -area / 3;
			end
		end
	end
	A = sparse(iis, jjs, vs, nvertices, nvertices, "sum");
endfunction

function [bti] = computeBTildI(A, b, nbInterior, L, U, p, q)
	bti = b(nbInterior + 1 : end) ...
		- A(nbInterior + 1 : end, 1:nbInterior) * (U \ (L \ b(1:nbInterior)(p)))(q);
endfunction

function [res] = multiplyByLocalSchurComplement(A, nbInterior, L, U, p, q, x)
	res = A(nbInterior + 1 : end, nbInterior + 1 : end) * x ...
		- A(nbInterior + 1 : end, 1:nbInterior) * (U \ (L \ (A(1:nbInterior, nbInterior + 1 : end) * x)(p)))(q);
endfunction

function [vertices, triangles, border, domains] = readMesh(fileName)

	fid = fopen (fileName, "r");

	% Vertices
	fgoto(fid, "Vertices");
	
	[nvertices] = fscanf(fid, "%d", "C");
	vertices = zeros(nvertices, 2);
	for vid=1:nvertices
		[vertices(vid, 1), vertices(vid, 2)] = fscanf(fid, "%f %f %d\n", "C");
	end

	% Edges
	border = zeros(nvertices, 1);
	fgoto(fid, "Edges");
	[nedges] = fscanf(fid, "%d", "C");
	for eid=1:nedges
		[v1, v2] = fscanf(fid, "%d %d %d\n", "C");
		border(v1) = 1;
		border(v2) = 1;
	end

	% Elements
	fgoto(fid, "Triangles");
	[ntriangles] = fscanf(fid, "%d", "C");
	triangles = zeros(ntriangles, 3);
	domains = zeros(ntriangles, 1);
	for tid=1:ntriangles
		[triangles(tid, 1), triangles(tid, 2), triangles(tid, 3), domains(tid)] = fscanf(fid, "%d %d %d %d\n", "C");
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

function [] = writeSolution(fileName, solution)
	fid = fopen (fileName, "w");
		fprintf(fid, "MeshVersionFormatted 1\n\nDimension 2\n\nSolAtVertices\n%d\n1 1\n\n", length(solution));
		for i=1:length(solution)
			fprintf(fid, "%e\n", solution(i));
		end
		fclose(fid);
endfunction