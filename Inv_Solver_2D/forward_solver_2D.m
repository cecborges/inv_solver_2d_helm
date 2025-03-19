function [field,density] = forward_solver_2D(N,GG,q,k,rhs)
%
% solver for the forward problem
% This function solves the forward problem using the integral equation formulation.
% The scattered field is represented as us = V rho,
% where V is the volume potential.
% input :
%	N : number of points in the domain
%	GG : NxNxN tensor with weights for the volume problem
%	q : NxNxN tensor representing the sound profile of the domain
%	k : wavenumber fo the problem
%	rhs : N^3 right-hand-side for the problem
% output :
%	field : scattered field us = V rho
%	density: NxNxN tensor with density rho of the volume potential
%%%This is out
%	fwinfo : variable with information about the solution
%		fwinfo(1) : time for the solution of the integral equation
%		fwinfo(2) : flag for convergence, see doc for bicgstab
%		fwinfo(3) : relative residue
%		fwinfo(4) : number of iterations
%		fwinfo(5) : time to calculate fielg in the domain us = V rho

% finding the density    
    tol = 1.0E-10;
    maxit = 5000;    
    
    [sol,~] = bicgstab(@(x)matvec_helmholtz_2D(x,GG,q,N,k),rhs,tol,maxit);

%finding the field
    density   = reshape(sol,N,N);
    field = volume_density_fast_2D(GG,density,N);    



