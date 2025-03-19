clear;
fprintf('This code checks the volume potential using a Gaussian.\n')

%frequency
k  = 5.0d0;
lambda = 2.0d0*pi/k;
Nw = ceil(1.0d0/lambda);
if (Nw < 5)
     Nw = 5;
end

% number of points used in the discretization
Nvec = 200;%[5; 10; 15; 20; 25; 30; 35; 40; 45; 50];
Nvec = Nvec*Nw;
fprintf('lambda    =%e\n',lambda)
fprintf('Nw        =%e\n',Nw)

%direction of incidence of the plane wave
d  = [1.0d0;0.0d0];

for it = 1 : length(Nvec) 

    %setting up domain
    N = Nvec(it);
    fprintf('Nr points =%e\n',N)
    h=1.0d0/N;
    x= 0.0:h:((N-1.0)*h);
    [X,Y] = meshgrid(x);

    % set-up for forward solver
    GG = volume_density_setup_2D(k,N);    

    % create Gaussian density and solution
    sigma = 0.05d0;
    sig2 = sigma * sigma;
    S = sqrt((X-0.5).^2+(Y-0.5).^2);

    %densg
    densg = exp(-S.*S/2.0/sig2)/sig2/(2*pi);
    
    % solving V*densg   
    Fieldg = volume_density_fast_2D(GG,densg,N);

    R1 = S(1,1);
    
    %adjust this number to increase precision
    Nt1 = 10^5;
    h = R1 / Nt1;
    t = 0 : h : R1;
    
    Int1 = 0;
    Int1 = Int1 + besselj(0,k*t(1))*exp(-t(1).^2/2/sig2)*t(1);
    Int1 = Int1 + besselj(0,k*t(end))*exp(-t(end).^2/2/sig2)*t(end);
    Int1 = Int1 + 2*sum(besselj(0,k*t(2:end-1)).*exp(-t(2:end-1).^2/2/sig2).*t(2:end-1));
    Int1 = Int1*h/2;

    Rmax = 10;
    Nt2 = 10^5;
    h = (Rmax-R1) / Nt2;
    t = R1 : h : Rmax;
    
    Int2 = 0;
    Int2 = Int2 + besselh(0,k*t(1))*exp(-t(1).^2/2/sig2)*t(1);
    Int2 = Int2 + besselh(0,k*t(end))*exp(-t(end).^2/2/sig2)*t(end);
    Int2 = Int2 + 2*sum(besselh(0,k*t(2:end-1)).*exp(-t(2:end-1).^2/2/sig2).*t(2:end-1));
    Int2 = Int2*h/2;
    
    solg = besselh(0,1,k*R1)/4/sig2*Int1 + besselj(0,k*R1)/4/sig2*Int2;
    solg = solg*1i;
    
    fprintf('The error at the point (0,0) is %d!\n', abs(Fieldg(1)-solg)./abs(solg))

end
