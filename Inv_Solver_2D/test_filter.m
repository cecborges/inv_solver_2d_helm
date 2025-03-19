clear
for ii = 1 : 5
N = 100*ii
fprintf('Nr points =%e\n',N)
h=1.0d0/N;
x= 0.0:h:((N-1.0)*h);
[X,Y] = meshgrid(x);


Ncoeffs = 3;
coeffs = rand(Ncoeffs);

q = q_domain(X,Y,11,Ncoeffs,coeffs);


cf = filter_adj(N,q ,Ncoeffs);

norm(cf(:)-coeffs(:))/norm(coeffs(:))

end

