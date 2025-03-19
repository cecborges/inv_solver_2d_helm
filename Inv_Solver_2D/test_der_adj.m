% test_adjoint
clear;
fprintf('This code checks the LS solution using self convergence.\n')

%frequency
k  = 5.0d0;
lambda = 2.0d0*pi/k;
Nw = ceil(1.0d0/lambda);
if (Nw < 5)
     Nw = 5;
end
fprintf('lambda    =%e\n',lambda)
fprintf('Nw        =%e\n',Nw)

%direction of incidence of the plane wave
d  = [1.0d0;0.0d0];

%setting sensors
sensors = [ 10 10 -10 -10; 10 -10 10 -10];

%setting up domain
N = 30*Nw;
fprintf('Nr points =%e\n',N)
h=1.0d0/N;
x= 0.0:h:((N-1.0)*h);
[X,Y] = meshgrid(x);

% domain info
dflag = 1;
Ncoeffs = 1;
coeffs = zeros(Ncoeffs,Ncoeffs);
coeffs = 0.0d0;
q = q_domain(X,Y,dflag);

% set-up for forward solver
GG = volume_density_setup_2D(k,N);      

% solving forward problem
[Fq, field_domain] = forward_problem_sd(X,Y,GG,k,d,sensors,q);
umeas = Fq;

u_tot = field_domain + exp(1i*k*(d(1)*X+d(2)*Y));

% choose dq = delta*q;
delta = 0.1;

dq = delta * q;

%solving forward problem F(q)
% [Fnew, ~] = forward_problem(X,Y,GG,k,d,sensors,q+dq);

%finding the derivative Jdq
Jdq = frechet_der_sd(dq,X,Y,GG,k,sensors,q,u_tot);

Jadj = frechet_der_adj_sd(umeas,X,Y,GG,k,sensors,q,u_tot);

%  Checking error by doing <Jdq,conj(umeas)> and <conj(Jadj),dq>
l2Ssum = sum(Jdq.*conj(umeas));

Jadj = reshape(Jadj,N,N);
Jadj = conj(Jadj).*dq;

l2Rsum = trap2d(Jadj);

fprintf('L2S=%e\n',l2Ssum)
fprintf('L2R=%e\n',l2Rsum)
fprintf('Error(L2R)=%e\n',abs(l2Ssum-l2Rsum)/abs(l2Rsum))

