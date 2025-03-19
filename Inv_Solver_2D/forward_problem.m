function [field_sensors, field_domain] = forward_problem(X,Y,GG,k,dir,sensors,q)
% function calculating the field in the domain discretization points
% and in the sensors for incoming direction d.
%
% Input :
%       X : NxN matrix with X coord of the points of the disc domain
%       Y : NxN matrix with Y coord of the points of the disc domain
%       GG : 2Nx2N matrix for the weights in the volume potential
%       k : wavenumber incoming wave
%       dir : 2xNd matrix with coordinates of the incoming directions
%            dir(1,:) -> coordinates x of direction id
%            dir(2,:) -> coordinates y of direction id
%       sensors : Struct with the sensors for each direction. It is of the 
%            form
%            sensors(id).coords : 2xNtgt matrix with coordinates sensors  
%            for direction id 
%            sensors(id).coords(1,:) -> coordinates x 
%            sensors(id).coords(2,:) -> coordinates y
%       q : NxN matrix with values of the sound profile
%
% Output :
%       field_sensors : Struct with the scattered field at the sensors for 
%            each direction. It has the form
%            field_sensors(id).field -> Nt x1 scattered field at the 
%            sensors points for direction id
%       field_domain: Struct with scattered field at discretization points 
%            (this is necessary for the derivative and its adjoint). It
%            has the form
%            field_domain(id).field -> NxN matrix with scattered field for
%            direction id
%

    %size of domain and number of directions
    Nd = size(dir,2);
    
    for id = 1 : Nd
        
        %setting up directions and sensors
        d= dir(:,id);
        sensors_dir = sensors(id).coords;
        
        % solving forward problem for d
        [field_sensors(id).field,field_domain(id).field] = ...
            forward_problem_sd(X,Y,GG,k,d,sensors_dir,q);
        
    end
    
return
