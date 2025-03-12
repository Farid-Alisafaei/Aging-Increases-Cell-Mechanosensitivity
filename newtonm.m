function [x,iter] = newtonm( x0 , E_v1 , E_v2 , E_p )

N = 100; % define max. number of iterations
epsilon = 1e-10; % define tolerance
maxval = 10000000000000.0; % define value for divergence
xx = x0; % load initial guess

while (N>0)
    JJ = jacob6x6( xx , E_v1 , E_v2 , E_p ) ;
    if abs(det(JJ))<epsilon
        return
    end
    f = f6( xx , E_v1 , E_v2 , E_p ) ;
    xn = xx - inv(JJ)*f;
    if norm(f)<epsilon
        x=xn;
        iter = 100-N;
        return
    end
    if norm(f)>maxval
        iter = 100-N;
        disp(['iterations = ',num2str(iter)]);
        error('Solution diverges');
        return;
    end
    N = N - 1;
    xx = xn;
end
error('No convergence after 100 iterations.');
abort;




