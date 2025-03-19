clear

for ii=1:1
    n = 100*2^(ii-1);
    A = 2*diag(ones(n,1))-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
    B = 2*diag(ones(n,1))+A;
    U = rand(n,n);
    [U,~] = qr(U);
    D = diag(1./((1:n).^2));
    V = rand(n,n);
    [V,~] = qr(V);
    C = U*D*V';
    fprintf('n=%d\n',n)
    fprintf('cond(A) =%d\n',cond(A));
    fprintf('cond(B) =%d\n',cond(B));
    fprintf('cond(C) =%d\n\n',cond(C));
    
    x = rand(n,1);
    bA = A * x;
    bB = B * x;
    bC = C * x;

    fprintf('For matrix A\n')    
    nmax = 1000;
    
    fprintf('For Jacobi\n')
    for jj = 1 : 4
        epsit = (10^-3) * (10^(-2*(jj-1)));
        [soljacobi,it] = jacobi_method(A,bA,bA,epsit,nmax);
        fprintf('Eps=%d\n',epsit)
        fprintf('iter=%d\n',it)
        fprintf('residual=%d\n',norm(A*soljacobi-bA)/norm(bA))
        fprintf('error=%d\n\n',norm(soljacobi-x)/norm(x))

    end
    fprintf('For Gauss-Seidel\n\n')
    for jj = 1 : 4
        epsit = (10^-3) * (10^(-2*(jj-1)));
        [soljacobi,it] = gs_method(A,bA,bA,epsit,nmax);
        fprintf('Eps=%d\n',epsit)
        fprintf('iter=%d\n',it)
        fprintf('residual=%d\n',norm(A*soljacobi-bA)/norm(bA))
        fprintf('error=%d\n\n',norm(soljacobi-x)/norm(x))

    end
    
    
    

    fprintf('For matrix B\n')
    
    fprintf('For Jacobi\n')
    for jj = 1 : 4
        epsit = (10^-3) * (10^(-2*(jj-1)));
        [sol,it] = jacobi_method(B,bB,bB,epsit,nmax);
        fprintf('Eps=%d\n',epsit)
        fprintf('iter=%d\n',it)
        fprintf('residual=%d\n',norm(B*sol-bB)/norm(bB))
        fprintf('error=%d\n\n',norm(sol-x)/norm(x))

    end
    fprintf('For Gauss-Seidel\n\n')
    for jj = 1 : 4
        epsit = (10^-3) * (10^(-2*(jj-1)));
        [sol,it] = gs_method(B,bB,bB,epsit,nmax);
        fprintf('Eps=%d\n',epsit)
        fprintf('iter=%d\n',it)
        fprintf('residual=%d\n',norm(B*sol-bB)/norm(bB))
        fprintf('error=%d\n\n',norm(sol-x)/norm(x))

    end

    fprintf('For matrix C\n')
    
    fprintf('For Jacobi\n')
    for jj = 1 : 4
        epsit = (10^-3) * (10^(-2*(jj-1)));
        [sol,it] = jacobi_method(C,bC,bC,epsit,nmax);
        fprintf('Eps=%d\n',epsit)
        fprintf('iter=%d\n',it)
        fprintf('residual=%d\n',norm(C*sol-bC)/norm(bC))
        fprintf('error=%d\n\n',norm(sol-x)/norm(x))

    end
    fprintf('For Gauss-Seidel\n\n')
    for jj = 1 : 4
        epsit = (10^-3) * (10^(-2*(jj-1)));
        [sol,it] = gs_method(C,bC,bC,epsit,nmax);
        fprintf('Eps=%d\n',epsit)
        fprintf('iter=%d\n',it)
        fprintf('residual=%d\n',norm(C*sol-bC)/norm(bC))
        fprintf('error=%d\n\n',norm(sol-x)/norm(x))

    end    
    
end
    
function [x,it] = jacobi_method(A,b,x0,tol,Nmax)
    Q = diag(diag(A));
    B = A-Q;
    it = 0;
    x = x0;
    res = norm(A*x-b)/norm(b);
    while res > tol
        x = Q \ (b-B*x);
        it = it+1;
        if it == Nmax       
           fprintf('Max numit reached!\n')
           break
        end    
        res = norm(A*x-b)/norm(b);
    
    end

end

function [x,it] = gs_method(A,b,x0,tol,Nmax)
    D = diag(diag(A));
    Q = tril(A);
    B = triu(A)-D;
    it = 0;
    x = x0;
    res = norm(A*x-b)/norm(b);
    while res > tol
        x = Q \ (b-B*x);
        it = it+1;
        if it == Nmax       
           fprintf('Max numit reached!\n')
           break
        end    
        res = norm(A*x-b)/norm(b);
    end
end

