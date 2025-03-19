function field = volume_density_target_2D(X,Y,k,targ,rho)  
% this function calculates the field in the target/sendor points
% We use u(sensor)=V*rho(sensor), where V is the volume potential and rho
% is the density of the field calculated with Lippmann-Schwinger. We
% calculate the field using trapezoidal rule in 2D in the domain.
%
% Input :
%       X : NxN matrix with X coord of the points of the disc domain
%       Y : NxN matrix with Y coord of the points of the disc domain
%       k : wavenumber
%       targ : 2 X Nt  vector with sensors coordinates. First row has x
%       coords and second row y coords
%       rho : density of the field for the volume potential
%
% Output :
%         field : scattered field calculated at the target/sensor points.
%
%

%     %setting up domain
    N = size(X,1);
    h = 1.0d0/N;
%     x = 0 : h : (N-1)*h;
%     [X,Y] = meshgrid(x);
    
    
    %setting up weights for trapezoidal rule
    %probably consider creating a trapezoidal function for those integrals
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
    
    % multiplying density by weights
    rho_aux = rho.*coefs;    
        
    % this would use the FMM2D package(the direct goes fast(must install
    % FMM2D)
    % I am just leaving the code here and doing the direct calculation
%     srcinfo.sources = [X(:)';Y(:)'];
%     rho_aux = rho_aux(:);
%     srcinfo.charges = transpose(rho_aux);
%     
%     % setting up variables
%     zk      = complex(k);
%     pgt     = 1;
%
%     % calculating field
%     UF = h2ddir(zk,srcinfo,targ,pgt); 
%     field = 1i*UF.pottarg/4.0;    
    
    % doing all calculations in MATLAB - thinking about using the above to
    % speed up - the code with h22dir is faster because it is in fortran,
    % but you need to install another package.
    
    Skxy = k * sqrt((targ(1,:)'-X(:)').^2+(targ(2,:)'-Y(:)').^2);

    field = 1i/4*sum(besselh(0,1,Skxy).*(rho_aux(:).'),2);
        
end


    
    
    