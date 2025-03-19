%derivative test
%this code tests the derivative for the forward operator.
%check that F(q+dq)-F(q)-Jdq=O(h^2)
% We check by seeing that if h goes from 0.1^j to 0.1^(j+1), the error
% above decreases 0.001.
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
d  = [1.0d0;0.0d0];

%setting sensors
sensors = [ 10 10 -10 -10; 10 -10 10 -10];

%setting up domain
N = 30*Nw;
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
[Fq, field_domain] = forward_problem_sd(X,Y,GG,k,d,sensors,q);

u_tot = field_domain + exp(1i*k*(d(1)*X+d(2)*Y));


%test of the derivative e(h) = F(q+delta*q)-(F(q)+Jdq)

for ii = 1 : 5

    % choose dq = delta*q;
    delta = 0.1^ii;

    dq = delta * q;

    %solving forward problem F(q)
    [Fnew, ~] = forward_problem_sd(X,Y,GG,k,d,sensors,q+dq);
    
    %finding the derivative Jdq
    Jdq = frechet_der_sd(dq,X,Y,GG,k,sensors,q,u_tot);
    
    e(:,ii)=Fnew-(Fq+Jdq);

end

for ii = 1 : 4
    
    a=diff(e(ii,:));
    fprintf('Error at sensors 1 =%d %d %d %d\n',abs(a(1,1)),abs(a(1,2)),...
        abs(a(1,3)),abs(a(1,4)))   

end

