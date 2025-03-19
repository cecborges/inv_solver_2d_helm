function cf = filter_adj(N,func,Ncoeffs)
% this function gets the function func in the grid [0,1]^2 with N points in
% each coordinate and projects it in the space of sine functions with
% Ncoeffs^2 coefficients.
% Some thought should be put to use a dst or fft here.
% For the moment use mash=1/10, and use your grid as a multiply of 10, this
% will give everything that you need.

    
    h = 1.0/N;
    x = 0:h:(N-1)/N;
    maxh = 1/10;%1/30;%N should be a multiple of this guy
    a=1/2-maxh;
    x1 = (x-1/2+a)/(2*a);          

    %this needs to be changed to a dst.
    % I didnt want to think about the translational step among other
    % things. If you set maxh =1/10, and  use N as a multiply of 10
    % everything works out perfectly
    %change this to the appropriate dst
    cf = zeros(Ncoeffs,Ncoeffs);
    for jc = 1:Ncoeffs
        sinj = sin(jc*pi*x1);
        sinj(x<=maxh|x>=1-maxh)=0;
        for ic = 1:Ncoeffs
            sini = sin(ic*pi*x1);
            sini(x<=maxh|x>=1-maxh)=0;

            func_b = sini' * sinj;

            func_val = func_b.*func;

            cf(ic,jc) = trap2d(func_val)/a/a;

        end
    end
                
end