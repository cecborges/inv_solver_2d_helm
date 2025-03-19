% inverse solver test
clear

%set-up data
kvec = 1 : 0.25 : 30;
Nk = length(kvec);
lambda = 2.0d0*pi./kvec;
Nw = ceil(1.0d0./lambda);
Nw(Nw<1) = 1;
Nvec = 30*Nw;%using 30 points per wavelength. Should be enough for this problem


%setting up struct for directions and sensors
for ik = 1: Nk
    
    k = kvec(ik);    
    %we use the same amount of data for each step
    Nd = ceil(4*Nw(ik));
    theta_d = 0:2*pi/Nd:((1-1/Nd)*2*pi);

    %sensors
    Nt = ceil(4*Nw(ik));
    Rt =10;
    theta_t = 0:2*pi/Nt:((1-1/Nt)*2*pi);

    dvec(ik).dir = [cos(theta_d);sin(theta_d)];
    for id = 1: Nd
        sensorsvec(ik).dir(id).coords = Rt*[cos(theta_t);sin(theta_t)];
    end
    
end

%generating data
for ik = 1 : Nk
    
    %setting up domain
    k = kvec(ik);
    N = Nvec(ik);
    fprintf('Nr points =%e\n',N)
    h=1.0d0/N;
    x= 0.0:h:((N-1.0)*h);
    [X,Y] = meshgrid(x);

    % domain info
    dflag = 1;
    q = q_domain(X,Y,dflag);

    % set-up for forward solver
    GG = volume_density_setup_2D(k,N);      
    
    %find field for each direction
    dir = dvec(ik).dir;
    Nd = length(sensorsvec(ik).dir);
    for id = 1 : Nd
        sensors(id).coords = sensorsvec(ik).dir(id).coords;
    end
    [field_sensors, field_domain] = forward_problem(X,Y,GG,k,dir,sensors,q);

    for id = 1 : Nd
        data(ik).dir(id).field = field_sensors(id).field;
    end

end

qtrue = q;

%inverting domain
%using conjudate gradient

%stopping parameters
itmax = 500;
eps_res = 1e-3;
eps_dq  = 1e-5;
eps_bd  = 1e-5;

for ik = 1: Nk
    %not worrying about changing the domain right now...but we should
    %consider bandlimited domains
    %setting up domain
    k = kvec(ik);
    lambda = 2.0d0*pi/k;
    Nw = ceil(1.0d0/lambda);    
    N = 20*Nw;
    h=1.0d0/N;
    x= 0.0:h:((N-1.0)*h);
    [X,Y] = meshgrid(x);
        
    Ncoeffs = floor(3*Nw); % it works with 2Nwv with 1% error and 3Nw with 0.6%
    if Ncoeffs<1
        Ncoeffs =1 ;
    end

    if ik ==1 
        
        qvec = zeros(Ncoeffs);
        q = zeros(size(X));    
    else
    
        qvecaux = qvec;
        qvec = zeros(Ncoeffs);
        qvec(1:Ncoeffs_old,1:Ncoeffs_old) = qvecaux;        
        q = q_domain(X,Y,11,Ncoeffs,qvec);
        
    end
    
    % number of directions
    Nd = length(sensorsvec(ik).dir);

    %preparing data for rhs
    data_fw =[];    
    for id = 1 : Nd
        data_fw = [data_fw; data(ik).dir(id).field];
    end
    
    % set-up for forward solver
    GG = volume_density_setup_2D(k,N);      
    
    %find forward of initial guess
    dir = dvec(ik).dir;    
    sensors = [];
    for id = 1 : Nd
        sensors(id).coords = sensorsvec(ik).dir(id).coords;
    end    

    [field_sensors, field_domain] = forward_problem(X,Y,GG,k,dir,sensors,q);
    Fq = [];
    for id = 1 : Nd       
        d = dir(:,id);
        Fq = [Fq; field_sensors(id).field];
        utot(id).field = field_domain(id).field + exp(1i*k*(d(1)*X+d(2)*Y));
    end

    it = 1;
    rhs = data_fw - Fq;
    norm_rhs = norm(rhs)./norm(Fq);
    norm_dq = eps_dq+1;
    flag_newton = true;


    sol_it(ik).iter(it).norm_rhs = norm_rhs;
    sol_it(ik).iter(it).domain = q;

    while(flag_newton)
             
        Jadj = lsqr_op_sine([real(rhs);imag(rhs)],'transp',Ncoeffs,X,Y,GG,k,sensors,q,utot);
        dqvec = reshape(Jadj,Ncoeffs,Ncoeffs);       
%         This is not an ideal method, but I didn't want to spend on looking 
%         for the best optimization method possible. I pretty much do steepest 
%         descent, and minimize in that direction. Supposed to work if you are 
%         really close to the solution
        objfunc1 = @(input)invobj(input,Nd,Ncoeffs,data_fw,X,Y,GG,k,dir,sensors,qvec,dqvec);
        [t,fval] = fminsearch(objfunc1,0.1);
        
        dqvec = t*dqvec;        
        qvec_old = qvec;
        qvec = qvec + dqvec;    

        
        norm_dq = norm(dqvec(:))/norm(qvec(:));
        sol_it(ik).iter(it).norm_dq = norm_dq;

        it = it + 1;
        q_old = q;
        q = q_domain(X,Y,11,Ncoeffs,qvec);

        [field_sensors, field_domain] = forward_problem(X,Y,GG,k,dir,sensors,q);
        Fq = [];
        for id = 1 : Nd        
            d = dir(:,id);
            Fq = [Fq; field_sensors(id).field];
            utot(id).field = field_domain(id).field + exp(1i*k*(d(1)*X+d(2)*Y));
        end

        norm_rhs_old = norm_rhs;
        rhs_old = rhs;
        rhs = data_fw - Fq;
        
        norm_rhs = norm(rhs)/norm(Fq);
        sol_it(ik).iter(it).domain = q;        

        ier = [];
        if (it>=itmax)
            flag_newton = false;
            ier = [ier 1];
            fprintf('Total number of it reached!\n')
        end

        if (norm_rhs<=eps_res)
            flag_newton = false;
            ier = [ier 2];
            fprintf('RHS small!\n')
        end

        if ((norm_rhs>=norm_rhs_old))
            flag_newton = false;
            qvec = qvec_old;
            sol_it(ik).iter(it).domain = q_old;
    	    sol_it(ik).iter(it).norm_rhs = norm_rhs_old;
            rhs = rhs_old;
            ier = [ier 3];
            fprintf('RHS is increasing!\n')
        end

    	if ((abs(norm_rhs-norm_rhs_old)<eps_bd) &&(norm_rhs<=norm_rhs_old))
            flag_newton = false;
            ier = [ier 4];
            fprintf('RHS bogged down!\n')
    	end

        if (norm_dq<=eps_dq)
            flag_newton = false;
            ier = [ier 5];
            fprintf('dq small!\n')
        end
        
        if ~flag_newton
            sol_it(ik).ier = ier;
            
            fprintf('For k=%d\n',k)
            fprintf('Nr iterations=%d\n',it)
            fprintf('Ncoeffs=%d\n', Ncoeffs)
            fprintf('||Rhs||=%d\n',norm(rhs))            
            fprintf('rhs=%d\n',norm(rhs)/norm(Fq))
            fprintf('ier=%d\n\n',ier')

        end
                
        Ncoeffs_old = Ncoeffs;
        
    end
    
end

%results
q1=q_domain(X,Y,1);
cf = filter_adj(N,q1,Ncoeffs);
q11=q_domain(X,Y,11,Ncoeffs,cf);
q22=q_domain(X,Y,11,Ncoeffs,qvec);
surf(q11)
figure;surf(q22)
fprintf('Error=%d\n',norm(cf(:)-qvec(:))/norm(cf(:)))
