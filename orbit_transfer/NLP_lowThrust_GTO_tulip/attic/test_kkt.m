N = 10;
C = spdiags([-ones(N,1) randn(N,1)*1e5], [0 1], N-1, N); % unstable bidiagonal
Q = speye(N); Q(1,1)=0; Q(N,N)=0;
KKT = [Q C'; C sparse(N-1, N-1)];
rhs = zeros(2*N-1, 1);
rhs(N) = 1; % Some normalization/forcing
sol = KKT \ rhs;
disp(norm(C*sol(1:N) - rhs(N+1:end)));
