N = 5; nN = N+1;
nz_A = 26; nz_B_E = 3; nz_B_H = 4;
% Phase E
nnz_E = N * (2*nz_A + 2*nz_B_E) + (N+1)*3 + N*7;
% Phase H
nnz_H = N * (2*nz_A + 2*nz_B_H) + (N+1)*4;
fprintf('E: %d, H: %d\n', nnz_E, nnz_H);
