function out_field = matvec_helmholtz_2D(X,GG,q,N,k)
% This function calculates the matvec for the integral solution of the 
% Helmholtz equation. We have us = V rho. The integral equation is
% (-I+k^2 q V) rho = rhs.
% input :
%	X : N^2x1 input vector for the matvec of the Helmholtz equation
%	GG : tensor for solving the forward problem
%	q : NxN tensor representing the sound profile of the domain
%	N : number of discretization points
%	k : wavenumber	
% output:
% 	out_field : value of the matvec (-I+k^2 q V)X

%matvec for solution of the system
    rho = reshape(X,N,N);         
    [Field] = volume_density_fast_2D(GG,rho,N);    
    Field_aux = k.^2.*q.*Field - rho;
    out_field =Field_aux(:);    

return