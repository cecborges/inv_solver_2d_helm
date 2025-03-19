function Jadj = frechet_der_adj_sd(x,X,Y,GG,k,sensors,q,u)
% This function calculates the adjoint operator. The explanation of the
% calculation i below.
%  To calculate the adjoint operator applied to a vector f, we 
%  need to solve [Laplace+k^2(1+q)]w = -k^2 \xi(f,tgt)
%      = -k^2 \sum_j=1^Ntgt f(x_j) \delta(x-x_j)
%  To find w we should write w = w1 + w2
%  Take w1 = -k^2 \sum_j=1^Ntgt f(x_j) Gkadj(x-xj)
%  where Gkadj(x-xj)=Gk(x-xj)' and
%  [Laplace+k^2] Gkadj(x-xj)=\delta(x-xj)
%  We have [Laplace+k^2(1+q)]w1=
%     = -k^2 \sum_j=1^Ntgt f(x_j)\delta(x-xj) + k^2qw1
%  then [Laplace+k^2(1+q)]w2 = -k^2qw1. Summing we obtain
%   [Laplace+k^2(1+q)](w1+w2)=
%         = -k^2 \sum_j=1^Ntgt f(x_j) \delta(x-xj)
%  Steps:
%  1 ->  Calculate w1. Note that Gkadj is not Gk.
%  2 -> Solve Lippmann-Schwinger for w2
%  3 -> make w = w1 + w2
%
% input: 
%       x : Nt (number of sensors) x 1 vector to apply the adjoint operator
%       X : NxN matrix with X coord of the points of the disc domain
%       Y : NxN matrix with Y coord of the points of the disc domain
%       GG : 2Nx2N matrix of weights for the LS equation
%       k : wavenumber
%       sensors: 2 X Nt  vector with sensors coordinates. First row has x
%       coords and second row y coords
%       q : NxN matrix representing the sound profile of the domain
%       u : NxN vector with total field u^i+u^s
%
% output:
%       Jadjx : N^2x1 vector with the value of the adjoint applied in the
%       vector x evaluated at the points [X(:) Y(:)]
%


%First problem
Ntgt = size(sensors,2);
N = size(q,1);
k2 = k*k;

w1 = 0.0d0;
for ii = 1 : Ntgt
    
    S = sqrt((X-sensors(1,ii)).^2+(Y-sensors(2,ii)).^2);
    w1 = w1 + x(ii) * conj(1i/4*besselh(0,1,k*S));
    
end

w1 = -k2*w1;

%calculate w2
%rhs
rhs2 = -k2*q.*w1;

[w2,~] = forward_solver_2D(N,conj(GG),q,k,rhs2(:));

w= w1 + w2;

Jadj = -w(:).*conj(u(:));

end









