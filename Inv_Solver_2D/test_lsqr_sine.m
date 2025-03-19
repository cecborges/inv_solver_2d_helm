% test_adjoint
clear;
fprintf('This code checks the LS solution using self convergence.\n')

%frequency
k  = 5.0d0;
lambda = 2.0d0*pi/k;
Nw = ceil(1.0d0/lambda);
if (Nw < 5)
     Nw = 5;
end
fprintf('lambda    =%e\n',lambda)
fprintf('Nw        =%e\n',Nw)


%direction of incidence of the plane wave
Nd = 2;
theta_d = 0:2*pi/Nd:((1-1/Nd)*2*pi);

%sensors
Nt = 3;
Rt =10;
theta_t = 0:2*pi/Nt:((1-1/Nt)*2*pi);

%setting up struct for directions and sensors
dir = [cos(theta_d);sin(theta_d)];
for id = 1: Nd
    sensors(id).coords = Rt*[cos(theta_t);sin(theta_t)];
end

%generating data
N = 120;
fprintf('Nr points =%e\n',N)
h=1.0d0/N;
x= 0.0:h:((N-1.0)*h);
[X,Y] = meshgrid(x);

% domain info
dflag = 1;
q = q_domain(X,Y,dflag);

Ncoeffs = 2;
dqvec_coeffs = zeros(1,Ncoeffs^2);
dqvec_coeffs(1) = 1;
dqvec_coeffs = reshape(dqvec_coeffs,Ncoeffs,Ncoeffs);
dq = q_domain(X,Y,11,Ncoeffs,dqvec_coeffs);

% set-up for forward solver
GG = volume_density_setup_2D(k,N);      

%find field for each direction
[field_sensors, field_domain] = forward_problem(X,Y,GG,k,dir,sensors,q);

umeas = [];

dqvec = dqvec_coeffs(:);

for id = 1 : Nd
    d = dir(:,id);
    umeas = [umeas; field_sensors(id).field];   
    u_tot(id).field = field_domain(id).field + exp(1i*k*(d(1)*X+d(2)*Y));

end

Jdq = lsqr_op_sine(dqvec,'notransp',Ncoeffs,X,Y,GG,k,sensors,q,u_tot);
Ndata = length(Jdq);
Jdq = Jdq(1:Ndata/2,1)+1i*Jdq(Ndata/2+1:end,1);
l2Ssum = real(sum(Jdq.*conj(umeas)));


umeas = [real(umeas); imag(umeas)];
Jadj = lsqr_op_sine(umeas,'transp',Ncoeffs,X,Y,GG,k,sensors,q,u_tot);

Jadjvec = reshape(Jadj,Ncoeffs,Ncoeffs);
dqadj = q_domain(X,Y,11,Ncoeffs,Jadjvec);
l2Rsum = trap2d(dqadj.*dq);

fprintf('Error =%d\n',abs(l2Ssum-l2Rsum)/abs(l2Ssum))











