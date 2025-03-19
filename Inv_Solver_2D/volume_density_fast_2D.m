function [Field]=volume_density_fast_2D(GG,rho,N)
% This function is used to compute the volume potential V\rho
% input :
%	GG : 2Nx2N tensor with weights for the volume problem
%	N : number of points in the domain
%	rho : density for calculating the volume potential
% output :
%	Field: NxN matrix with volume potential applied in the density
%
% Note: GG should be of the form gpuArray for the fftn to be fast.
%       In the main code, add the command 
%                      GG = gpuArray(GG);
%       just after calculating it.


%     if canUseGPU
% 
% 	Field = vol_dens_fast_2D_gpu(GG,rho,N);
% 
%     else    

	Field = vol_dens_fast_2D(GG,rho,N);

%     end

end

function Field = vol_dens_fast_2D(GG,rho,N)

        FT_source = fftn(rho,[2*N,2*N]);
        Field = ifftn(FT_source.*GG);
        Field = Field(1:N,1:N);

end

function FieldGPU = vol_dens_fast_2D_gpu(GG,rho,N)
       rhog = gpuArray(rho);
       rhog = fftn(rhog,[2*N,2*N]);
       rhog = ifftn(rhog.*GG);
       FieldGPU = gather(rhog(1:N,1:N));

end