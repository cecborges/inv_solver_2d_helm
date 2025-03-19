function output = invobj(input,Nd,Ncoeffs,data_fw,X,Y,GG,k,dir,sensors,qvec,dqvec)

    qsol = qvec+input*dqvec;
    q = q_domain(X,Y,11,Ncoeffs,qsol);
    
    [field_sensors, field_domain] = forward_problem(X,Y,GG,k,dir,sensors,q);
    Fq = [];
    for id = 1 : Nd       
        d = dir(:,id);
        Fq = [Fq; field_sensors(id).field];
        utot(id).field = field_domain(id).field + exp(1i*k*(d(1)*X+d(2)*Y));
    end

    it = 1;
    output = 0.5*norm(data_fw - Fq)^2;





