function [output] = lsqr_op_domain(input,opt,X,Y,GG,k,sensors,q,u_tot)
% This function calculates J h and J^T h and is used for the LSQR method.
% This function can also be used separately to calculate J dq and J^T dq.
% (Maybe I should change the name)
% This function is personalized using a sine basis.
% input :
%	input : vector input h for multiplication by J or J^T.
%	opt : option to be used in the multiplication
%		opt == notransp : do Jh
%		otherwise : J^T h
%	Ncoeffs :  number of coefficients used in the basis for representing 
%		the domain. We have a total of Ncoeffs^2 total coefficients
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
% output :
%	output : vector with J input or J^T input depending on the opt chosen.

if strcmp(opt,'notransp')    
    %input function values in the domain in R^{N^2}
    %output (2Nd*Ntgt,1) in R       
    inp_domain = input;    

    outputaux = frechet_der(inp_domain,X,Y,GG,k,sensors,q,u_tot);
%     size(outputaux)

    Nd = length(sensors);
    output1 = [];
    for id = 1 : Nd

        output1 = [output1; outputaux(id).field];

    end

    output = [real(output1);imag(output1)];


else
    %input (2Nd*Ntgt,1) in R
    %output (N^2,1) in R^{N^2}
    N = size(X,1);
    Ndata = length(input);    
    inputaux = input(1:Ndata/2,1)+1i*input(Ndata/2+1:end,1);    

    outputaux = frechet_der_adj(inputaux,X,Y,GG,k,sensors,q,u_tot);

    Nd = length(sensors);
    outputaux1 = zeros(size(X));    
    for id = 1 : Nd

        outputaux1 = outputaux1 + reshape(outputaux(id).field,N,N);

    end

    outputaux = real(outputaux1);
    

    output = outputaux(:);
    
end