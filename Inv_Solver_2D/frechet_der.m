function Jdq = frechet_der(dq,X,Y,GG,k,sensors,q,u_tot)
% function to calculate the derivative of the forward operator. It
% applies the single direction function for each of the given directions.
% This function was desigend to provide a matrix vector
% multiplication to be used with an iterative method.
% Input :
%       dq : N^2x1 update of the domain. 
%       X : NxN matrix with X coord of the points of the disc domain
%       Y : NxN matrix with Y coord of the points of the disc domain
%       GG : 2Nx2N matrix of weights for the LS equation
%       k : wavenumber
%       sensors : Struct with the sensors for each direction. It is of the 
%            form
%            sensors(id).coords : 2xNtgt matrix with coordinates sensors  
%            for direction id 
%            sensors(id).coords(1,:) -> coordinates x 
%            sensors(id).coords(2,:) -> coordinates y
%       q : NxN matrix representing the sound profile of the domain
%       u_tot : Struct with the total field  for each direction. It is of 
%            the form 
%            u_tot(id).field : NxN matrix with the total field 
%            (scattered+ incident) in the discretization points in the
%             domain for direction id
%
% Output :
%       Jdq : Struct with the derivative of the field at the sensors for 
%            each direction using update dq. It has the form
%            Jdq(id).field -> Nsensrors for direction id x 1 derivative at the 
%            sensors points for direction id
%

    %size of domain and number of directions
    Nd = length(sensors);
    
    for id = 1 : Nd
        
        %setting up directions and sensors
        sensors_dir = sensors(id).coords;
        utot_dir = u_tot(id).field;
        % calculating de derivative for each direction
        Jdq(id).field = frechet_der_sd(dq,X,Y,GG,k,sensors_dir,q,utot_dir);

    end
    
return