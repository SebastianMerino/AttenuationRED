
% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function u = IRLS_TV(b,A,mu,M,N,tol,mask,minimask)

[u,~] = cgs2(A'*A,A'*b,1e-6,20);
G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_isotropic(u,M,N,minimask);

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dx = kron(speye(N),D);

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dy = kron(D,speye(M));

D = [Dx' Dy']';

ite_irls = 0;
error = 1;

while error > tol

    X = reshape(u,M,N);
    ite_irls = ite_irls + 1;
    Dh = diff(X,[],1);
    Dh = [Dh;zeros(1,N)];
    Dv = diff(X,[],2);
    Dv = [Dv zeros(M,1)];

    P = Dh.^2 + Dv.^2;
    eps = 0.01;
    P = 2*sqrt(P.^2 + eps^2);
    P = P.^(-0.5);
    P = P(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(P,0,omega);
    W = kron(speye(2),omega);

    AtA = A'*A;
    [u] = cgs2( AtA + mu*D'*W*D, A'*b, 1e-6 , 20, u );

    G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_isotropic(u,M,N,minimask);
    error = abs(G(ite_irls+1) - G(ite_irls));

end
end
