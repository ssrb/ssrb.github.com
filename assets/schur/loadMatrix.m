filename = "poisson2D.A";
fid = fopen (filename, "r");

[n,m,isSymmetric, nz] = fscanf(fid, "%u %u %u %u\n", "C");
if isSymmetric == 1
	nz *= 2;
end

A = spalloc (m, n, nz)
for coeff=1:nz
    [i, j, a] = fscanf(fid, "%u %u %f", "C");
    A(i,j) = a;
    if isSymmetric == 1
    	A(j,i) = a;
    end
end
fclose(fid);

filename = "poisson2D.b";
fid = fopen (filename, "r");

[n] = fscanf(fid, "%u\n", "C");
b = zeros(n,1);
for i=1:n
    b(i) = fscanf(fid, "%f", "C");
end
fclose(fid);

