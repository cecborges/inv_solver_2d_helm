function Jdq = frechet_der_sd(dq,X,Y,GG,k,sensors,q,u_tot)
% function to calculate the derivative of the forward operator for one
% single direction. This function was dexigend to provide a matrix vector
% multiplication to be used with an iterative method.
% Input :
%       dq : N^2x1 update of the domain. 
%       X : NxN matrix with X coord of the points of the disc domain
%       Y : NxN matrix with Y coord of the points of the disc domain
%       GG : 2Nx2N matrix of weights for the LS equation
%       k : wavenumber
%       sensors: 2 X Nt  vector with sensors coordinates. First row has x
%       coords and second row y coords
%       q : NxN matrix representing the sound profile of the domain
%       u_tot : NxN matrix with the total field (scattered+ incident) in
%       the discretization points in the domain
%
% Output :
%       Jdq : Ntx1 vector application of the derivative operator in the 
%       vector dq.
%

    % set-up rhs, using total field function
    N = size(X,1);
    dq = reshape(dq,N,N);
    k2 = k*k;
    rhs = -k2*dq.*u_tot;
    rhs=rhs(:);

    % solving forward problem     
    [~,density] = forward_solver_2D(N,GG,q,k,rhs);
    
    % calculating field at sensors
    Jdq = volume_density_target_2D(X,Y,k,sensors,density);
    
    return


