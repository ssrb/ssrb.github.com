function mpiSchur()

	pkg load mpitb;

	MPI_Init;

	[mpistatus, mpirank] = MPI_Comm_rank(MPI_COMM_WORLD);
	[mpistatus, mpicommsize] = MPI_Comm_size(MPI_COMM_WORLD);

	if mpirank == 0
		profile on;
	end

	[vertices, triangles, border, interface, vids] = readMesh("lyttelton.mesh", "interface2.vids", mpirank + 1);
	
	% Assemble *local* stiffness matrix and load vector
	[A, b] = assemble(vertices, triangles, border);

	nbInterior = length(b) - length(interface);

	% LU-factorize the interior of the subdomains, we're going to reuse this everywhere 
	[L, U, p, tmp] = lu(A(1:nbInterior, 1:nbInterior), 'vector');
	q(tmp) = 1:length(tmp);

	% We solve for the trace first: we need to compute the second member of the Schur complement system
	% Local contribution
	bti = computeBTildI(A, b, nbInterior, L, U, p, q);
	bt = zeros(size(bti));

	% Sum contributions
	mpistatus = MPI_Allreduce(bti, bt, MPI_SUM, MPI_COMM_WORLD);

	epsilon = 1.e-30;
	maxIter = 600;

	% Solve for the trace doing distributed gradient descent
	trcSol = pcg(	@(x) parallelMultiplyBySchurComplement(A, nbInterior, L, U, p, q, x), ...
					bt, epsilon, maxIter);
	
	% Solve for the local interior
	localISol = (U \ (L \ (b(1:nbInterior) - A(1:nbInterior, nbInterior + 1 : end) * trcSol)(p)))(q);

	% Consolidate solutions
	% Gather sizes of the subdomains
	allNbInterior = zeros(mpicommsize, 1);
	MPI_Gather(nbInterior, allNbInterior, 0, MPI_COMM_WORLD);

	sumNbInterior = sum(allNbInterior);

	disps = cumsum([0; allNbInterior(1:end-1)]);

	% Concatenate local => global mappings of the unknowns
	globalISol = zeros(sumNbInterior, 1);
	MPI_Gatherv(localISol, globalISol, allNbInterior, disps, 0 , MPI_COMM_WORLD);
	
	% Concatenate solutions
	allVids = zeros(sumNbInterior, 1);
	vids = vids(1:nbInterior);
	MPI_Gatherv(vids, allVids, allNbInterior, disps, 0, MPI_COMM_WORLD);

	if mpirank == 0
		% Reorder the global solution
	 	solution = zeros(sumNbInterior + length(interface), 1);
	 	solution(allVids) = globalISol;
	 	solution(interface) = trcSol;
	 	writeSolution("lyttelton.sol", solution);
	end

	MPI_Finalize;

	if mpirank == 0
		profile off;
		prof = profile ("info");
		profshow (prof);
	end

endfunction

function [A, b] = assemble(vertices, triangles, border)
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
	areas4 = 4 * areas;

	val(:, 1) = sum(u.*u, 2) ./ areas4;
	val(:, 2) = sum(u.*v, 2) ./ areas4;
	val(:, 3) = sum(u.*w, 2) ./ areas4;
	val(:, 5) = sum(v.*v, 2) ./ areas4;
	val(:, 6) = sum(v.*w, 2) ./ areas4;
	val(:, 9) = sum(w.*w, 2) ./ areas4;
	val(:, [ 4 , 7 , 8]) = val(:, [ 2 , 3 , 6 ]);

	col = triangles(:, [1 1 1 2 2 2 3 3 3]);
	row = triangles(:, [1 2 3 1 2 3 1 2 3]);

	A = sparse(col(:), row(:), val(:), nvertices, nvertices);
	% Dirichet penalty
	A += spdiags(1e30 * !!border, [0], nvertices, nvertices);

	b = zeros(nvertices, 1);
	for tid=1:ntriangles
		b(triangles(tid,:)) += -areas(tid) / 3;
	end

endfunction

% function [A, b] = assemble(vertices, triangles, border)
% 	nvertices = length(vertices);
% 	ntriangles = length(triangles);

% 	iis = [];
% 	jjs = [];
% 	vs = [];

% 	b = zeros(nvertices, 1);
% 	for tid=1:ntriangles
% 		q = zeros(3,2);
% 		q(1,:) = vertices(triangles(tid, 1), :) - vertices(triangles(tid, 2), :);
% 		q(2,:) = vertices(triangles(tid, 2), :) - vertices(triangles(tid, 3), :);
% 		q(3,:) = vertices(triangles(tid, 3), :) - vertices(triangles(tid, 1), :);
% 		area = 0.5 * det(q([1,2], :));
% 		for i=1:3
% 			ii = triangles(tid,i);
% 			if !border(ii)
% 				for j=1:3
% 					jj = triangles(tid,j);
% 					if !border(jj)
% 						hi = q(mod(i, 3) + 1, :);
% 						hj = q(mod(j, 3) + 1, :);
						
% 						v = (hi * hj') / (4 * area);

% 						iis = [iis ii];
% 						jjs = [jjs jj];
% 						vs = [vs v];
% 					end
% 				end
% 				b(ii) += -area / 3;
% 			end
% 		end
% 	end
% 	A = sparse(iis, jjs, vs, nvertices, nvertices, "sum");
% endfunction

function [bti] = computeBTildI(A, b, nbInterior, L, U, p, q)
	bti = b(nbInterior + 1 : end) ...
		- A(nbInterior + 1 : end, 1:nbInterior) * (U \ (L \ b(1:nbInterior)(p)))(q);
endfunction

function [res] = parallelMultiplyBySchurComplement(A, nbInterior, L, U, p, q, x)
	local = A(nbInterior + 1 : end, nbInterior + 1 : end) * x ...
		- A(nbInterior + 1 : end, 1:nbInterior) * (U \ (L \ (A(1:nbInterior, nbInterior + 1 : end) * x)(p)))(q);
	res = zeros(size(local));
	MPI_Allreduce(local, res, MPI_SUM, MPI_COMM_WORLD);
endfunction

function [vertices, triangles, border, interface, vids] = readMesh(fileName, interfaceName, domainIdx)

	interface = readInterface(interfaceName);

	fid = fopen (fileName, "r");

	% Elements
	fgoto(fid, "Triangles");
	[ntriangles] = fscanf(fid, "%d", "C");
	triangles = [];
	for tid=1:ntriangles
		[v1, v2, v3, domain] = fscanf(fid, "%d %d %d %d\n", "C");
		if domain == domainIdx
			triangles = [triangles; v1 v2 v3];
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
	border = zeros(length(vids), 1);

	currentLine = 1;
	for vid=1:length(sortedVertexIds)
		globalId = sortedVertexIds(vid);
		fskipLines(fid, globalId - currentLine);
		currentLine = globalId;
		localId = globalToLocal(globalId);
		[vertices(localId, 1) vertices(localId, 2) border(localId)] = fscanf(fid, "%f %f %d\n", "C");
		++currentLine;
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

% function x = conjgrad(A, nbInterior, L, U, pp, qq , bb, tol, maxiter, mpirank)
% 	x = 0 * bb;
%     r=bb;
%     p=r;
%     rsold=r'*r;
 
%     for iter=1:maxiter
%         Ap = parallelMultiplyBySchurComplement(A, nbInterior, L, U, pp, qq, p);
% 		if mpirank == 0
% 			dlmwrite(sprintf("schurIter%d.txt", iter), Ap, 'delimiter', '\n');
% 		end
%         alpha=rsold/(p'*Ap);
%         x=x+alpha*p;
%         r=r-alpha*Ap;
%         rsnew=r'*r;
%         if sqrt(rsnew)<tol
%               break;
%         end
%         p=r+rsnew/rsold*p;
%         rsold=rsnew;
%     end

% endfunction