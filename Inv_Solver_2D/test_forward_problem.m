%forward problem test
clear;
fprintf('This code checks the LS solution using self convergence.\n')

%frequency
k  = 5.0d0;
lambda = 2.0d0*pi/k;
Nw = ceil(1.0d0/lambda);
if (Nw < 5)
     Nw = 5;
end

% number of points used in the discretization
Nvec = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50];
Nvec = Nvec*Nw;
fprintf('lambda    =%e\n',lambda)
fprintf('Nw        =%e\n',Nw)

%direction of incidence of the plane wave
d  = [1.0d0;0.0d0];

%setting sensors
sensors = [ 10 10 -10 -10; 10 -10 10 -10];

%iterations to check the self convergence
for it = 1 : length(Nvec)

      %setting up domain
      N = Nvec(it);
      fprintf('Nr points =%e\n',N)
      h=1.0d0/N;
      x= 0.0:h:((N-1.0)*h);
      [X,Y] = meshgrid(x);
      
      % domain info
      dflag = 1;
      Ncoeffs = 1;
      coeffs = zeros(Ncoeffs,Ncoeffs);
      coeffs = 0.0d0;
      q = q_domain(X,Y,dflag);

      % set-up for forward solver
      GG = volume_density_setup_2D(k,N);      

      % solving forward problem
      [field_sensors, field_domain] = forward_problem_sd(X,Y,GG,k,d,sensors,q);
      
      solaux(it) = field_domain(1,1);
      maxsol(it) = max(abs(field_domain(:)));
      
      solsen(it) = field_sensors(1,1);
      maxsen(it) = max(abs(field_sensors(:)));
end
      
for it =1 : 5
       fprintf('Self-convergence for domain values!\n')
       fprintf('For N =%e\n',Nvec(it))
       fprintf('AErr  =%e\n',abs(solaux(6)-solaux(it)))
       fprintf('RErr  =%e\n',abs(solaux(6)-solaux(it))/abs(solaux(6)))
       fprintf('RErr1 =%e\n\n',abs(solaux(6)-solaux(it))/maxsol(6))
       
       fprintf('Self-convergence for sensor values!\n')
       fprintf('For N =%e\n',Nvec(it))
       fprintf('AErr  =%e\n',abs(solsen(6)-solsen(it)))
       fprintf('RErr  =%e\n',abs(solsen(6)-solsen(it))/abs(solsen(6)))
       fprintf('RErr1 =%e\n\n',abs(solsen(6)-solsen(it))/maxsen(6))
       
end
