function q = q_domain(X,Y,dflag,Ncoeffs,coeffs)
    
% This function calculates the values of the sound profile in the
% domain input:
%	X : X coordinates to evaluate the sound profile
%   Y : Y coordinates to evaluate the sound profile
%	dflag : indicates the type of the domain chosen
%		1 : gaussian 
%		2 : sphere
%		11 : sine series around 1/2 with length 2a
%		    a = 1/2 - maxh, where maxh is the maximum discretization step
%		    used in the forward problem for all frequencies. (Since we are 
%		    using the minimum amount of points of 30, set maxh=1/30)
%   Ncoeffs : number of coefficients
%   coeffs  : Ncoeffs x Ncoeffs matrix of coefficients
%
%   output:
%       q : evaluated at (x,Y)
%    

    if (dflag == 1)

        fprintf('Gaussian!\n')             
        x0 = 0.5d0;   
        y0 = 0.5d0;
        sigma_q = 0.05d0;
        R = sqrt((X-x0).^2 + (Y-y0).^2);
        q = 0.1d0*exp(-R.*R/2.0d0/sigma_q/sigma_q);
        return
        
    elseif (dflag == 2)
        
        fprintf('Sphere!\n')
        x0 = 0.5d0;      
        y0 = 0.5d0;
        radius = 0.5d0;
        R = sqrt((X-x0).^2 + (Y-y0).^2 );
        
        q=zeros(size(X));
        q(R<radius) = 0.5d0;          
        return
                

    elseif (dflag == 11)

        % fprintf('Evaluating a sine series on [0.5-a, 0.5+a]\n')
        % fprintf('sine series around 1/2 with length 2a!\n')
        N = size(X,1);
%         ind = ones(N);
        h = 1.0/N;
        x = 0:h:(N-1)/N;

        % choose here the h for the smallest problem being solved
        % if you choose an maxh smaller than the
        % h for the forward problem, it breaks.
        maxh = 1/10;%1/30;
        a=1/2-maxh;
        x1 = (x-1/2+a)/(2*a);          
%         ind(X<=1/2-a|X>=1/2+a)=0;
%         ind(Y<=1/2-a|Y>=1/2+a)=0;

        %this can be much faster
        q = zeros(size(X));
        for jc = 1:Ncoeffs
            sinj = sin(jc*pi*x1);
            sinj(x<=maxh|x>=1-maxh)=0;

            for ic = 1:Ncoeffs
                sini = sin(ic*pi*x1);
                sini(x<=maxh|x>=1-maxh)=0;

                vals = sini' * sinj;                
                q = q + coeffs(ic,jc)*vals;                    
            end
        end

        % and finally hard threshold the q, need to understand why
        % this is
        % necessary...
%         q = q.*ind;
        
        %%%%
        %%%%For sines begin
% % % %        N=params(3);        coefs=params(4:end);
% % % %        a=pi/2;%size of domain
% % % % 
% % % %        [n1,n2]=size(XX1);
% % % % 
% % % %        vec1=(1:N)';
% % % %     %    nx=kron(vec1,XX1/2+pi/2);
% % % %     %    ny=kron(vec1',XX2/2+pi/2);
% % % %        nx=pi/(2*a)*kron(vec1,XX1+a);
% % % %        ny=pi/(2*a)*kron(vec1',XX2+a);
% % % %        aux1=repmat(nx,1,N);
% % % %        aux2=repmat(ny,N,1);
% % % %        saux1=sin(aux1);
% % % %        saux2=sin(aux2);
% % % %     %    saux1=exp(1i*aux1);
% % % %     %    saux2=exp(1i*aux2);
% % % %     %    saux1=cos(aux1);
% % % %     %    saux2=cos(aux2);
% % % % 
% % % %        paux=saux1.*saux2;
% % % %        caux=reshape(coefs',[N,N]);
% % % %        caux1=kron(caux,ones(n1,n2));
% % % % 
% % % %        paux1=paux.*caux1;
% % % %        b=zeros(n1,n2);
% % % %        for ii=1:N
% % % %           for jj=1:N
% % % %             b=b+paux1(1+(ii-1)*n1:ii*n1,1+(jj-1)*n2:jj*n2);
% % % %           end
% % % %        end
% % % %        b(abs(XX1)>a)=0; b(abs(XX2)>a)=0;
        
        
        
        return     
            
    else
        
        fprintf('Option not defined!\n')
        q = 0.0d0;
        return    
    end
    

end

