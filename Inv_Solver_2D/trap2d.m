function out = trap2d(f)
% function to calculate the integral of f in the interval [0,1]^2 using the
% trapezoidal rule.
% input :
%       f : NxN function in the [0,1]^2 interval
% output :
%       out : integral of f
%
N = size(f,1);
h = 1.0d0/N;
coefs = zeros(N,N);
c3 = h * h;    
c2 = c3 * 2.0d0 / 4.0d0;
c1 = c3 / 4.0d0;
coefs(:,:) = c3;

% Vertices
coefs(1,1) = c1;
coefs(1,N) = c1;
coefs(N,1) = c1;
coefs(N,N) = c1;

% Edges
coefs(1,2:N-1) = c2;    
coefs(N,2:N-1) = c2;
coefs(2:N-1,1) = c2;
coefs(2:N-1,N) = c2;
wf = f.*coefs;  

out = sum(wf(:));
