function [field_sensors, field_domain] = forward_problem_sd(X,Y,GG,k,d,sensors,q)
% function calculating the field in the domain discretization points
% and in the sensors for incoming direction d.
%
% Input :
%       X : NxN matrix with X coord of the points of the disc domain
%       Y : NxN matrix with Y coord of the points of the disc domain
%       GG : 2Nx2N matrix for the weights in the volume potential
%       k : wavenumber incoming wave
%       d : direction incoming wave
%       sensors : 2xN matrix with coords of the sensors. First row has x
%       coords and Second row has y coords
%       q : NxN matrix with values of the sound profile
%
% Output :
%       field_sensors : Nt x1 scattered field at the sensors points
%       field_domain: NxN matrix with scattered field at discretization
%       points (this is necessary for the derivative and its adjoint)
%

    % set-up rhs, using plane wave function
    k2 = k*k;
    rhs = -k2*q.*exp(1i*k*(d(1)*X+d(2)*Y));
    rhs=rhs(:);

    % solving forward problem 
    N = size(X,1);
    [field_domain,density] = forward_solver_2D(N,GG,q,k,rhs);
    
    % calculating field at sensors
    field_sensors = volume_density_target_2D(X,Y,k,sensors,density);
    
    return

      
      
