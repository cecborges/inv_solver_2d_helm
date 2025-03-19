function Field = volume_density_setup_2D(k,N)
%
% This function sets up the density for the volume solver.
% input :
%	k : wavenumber
%	N : number of discretization points
% output :
%	Field : 2Nx2N matrix with values for the volume potential kernel
%
% setup dont change those. This is for the cube domain.
%

    h = 1.0/N;
    R_max = 1.5;


    GG = fftn_green_setup(k,h,4*N+1,R_max);
    Field = fftshift(ifftn(ifftshift(GG)))/(h^2);
    Field = Field(N+1:3*N,N+1:3*N)*(N^2);
    Field = fftshift(fftn(ifftshift(Field)))*(h^4);
    Field = fftshift(Field);

return
        
function [GG] = fftn_green_setup(k,dx,N,R_max)

    df = 1 / dx / (N-1);
    indices = -(N-1)/2:((N-1)/2-1);
    fx  = indices*df;

    [FX,FY] = meshgrid(fx,fx);
    S = 2*pi*sqrt(FX.^2+FY.^2);
                
    GG = (1+1i*pi/2 * R_max*S.*besselj(1,R_max*S)*besselh(0,1,R_max*k) + ...
         -1i*pi/2 * R_max*k*besselj(0,R_max*S)*besselh(1,1,R_max*k))./...
         ((S-k).*(S+k));        
          
    indices = find(abs(S-k)<=10^-7);    
    %check these calculations
    GG(indices) = 1i*pi*R_max/4/k*(besselh(0,1,R_max*k)*(R_max*k*besselj(0,R_max*k)+...
        2*besselj(1,R_max*k)-R_max*k*besselj(2,R_max*k))/2+R_max*k*...
        besselj(1,R_max*k)*besselh(1,1,R_max*k));                  
            
return