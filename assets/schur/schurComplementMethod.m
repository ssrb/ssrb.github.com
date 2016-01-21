pkg load general;
load("-ascii", "domain1.vids");
load("-ascii", "domain2.vids");
load("-ascii", "interface2.vids");

p = [domain1; domain2; interface2];
loadMatrix;
A = A(p,p);
F = b(p);
%spy(A)

i1 = length(domain1);
i2 = i1 + length(domain2);

A11 = A(1:i1, 1:i1);
A1T = A(1:i1, i2+1:end);
A22 = A(i1+1:i2, i1+1:i2);
A2T = A(i1+1:i2, i2+1:end);
AT1 = A(i2+1:end, 1:i1);
AT2 = A(i2+1:end, i1+1:i2);
ATT = A(i2+1:end, i2+1:end);
U = zeros(length(A), 1);
clear A;

F1 = F(1:i1);
F2 = F(i1+1:i2);
FT = F(i2+1:end);
clear F;

%spy(ATT - AT1 * inv(A11) * A1T + AT2 * inv(A22) * A2T);

[L11, U11, p11, q11] = lu(A11, 'vector');
q11t(q11) = 1:length(q11);
clear A11;
clear q11;

[L22, U22, p22, q22] = lu(A22, 'vector');
q22t(q22) = 1:length(q22);
clear A22;
clear q22;

% Bad !
%ATT -= AT1 * (U11 \ (L11 \ A1T)) + AT2 * (U22 \ ( L22 \ A2T));
%spy(ATT);

FT -= AT1 * (U11 \ (L11 \ F1(p11)))(q11t) ...
	+ AT2 * (U22 \ (L22 \ F2(p22)))(q22t);

U(i2+1:end) = pcg(@(x) ...
	ATT * x ...
	- AT1 * (U11 \ (L11 \ (A1T * x)(p11)))(q11t) ...
	- AT2 * (U22 \ (L22 \ (A2T * x)(p22)))(q22t), FT, 1.e-12, 500);

clear ATT;
clear FT;

U(1:i1) = (U11 \ (L11 \ (F1 - A1T * U(i2+1:end))(p11)))(q11t);
U(i1+1:i2) = (U22 \ (L22 \ (F2 - A2T * U(i2+1:end))(p22)))(q22t);

pt(p) = 1:length(p);
U = U(pt);

filename = "poisson2D.sol";
fid = fopen (filename, "w");
fprintf(fid, "MeshVersionFormatted 1\n\nDimension 2\n\nSolAtVertices\n3735\n1 1\n\n");
for i=1:3735
	fprintf(fid, "%e\n", U(i));
end
fclose(fid);
