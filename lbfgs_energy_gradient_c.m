function p = lbfgs_energy_gradient_c(c, z, w, Iref, reg_rho, iter, niter)
    [dimx, dimy] = size(Iref);

    persistent g;
    persistent y; 
    persistent s; 
    persistent x; 
    
    persistent m = 5;
    
    persistent alpha;
    persistent beta;
    persistent rho;

    if iter == 1
        g = zeros(dimx, dimy, m+1);
        y = zeros(dimx, dimy, m);
        s = zeros(dimx, dimy, m);
        x = zeros(dimx, dimy, m+1);
        
        % Initialize arrays for L-BFGS variables
        alpha = zeros(1, m);
        beta  = zeros(1, m);
        rho   = zeros(1, m);
    endif


    % Initialize arrays to store estimate of x, the gradient and the
    % differences.
    x(:, :, m+1) = c; 
    g(:, :, m+1) = q = -1 * energy_gradient_c(c, z, w, Iref, reg_rho);
    
    s(:, :, m) = x(:, :, m+1) - x(:, :, m);
    y(:, :, m) = g(:, :, m+1) - g(:, :, m);
    
    if (iter <= 2)
        x(:, :, 1:m) = x(:, :, 2:m+1);
        g(:, :, 1:m) = g(:, :, 2:m+1);
        
        s(:, :, 1:m-1) = s(:, :, 2:m);
        y(:, :, 1:m-1) = y(:, :, 2:m);
        
        p = -q;
        return;
    endif

    % Get values of rho, s and y
    for i = m:-1:1
        yTs = innerprod(s(:,:,i), y(:,:,i));
        
        if (yTs == 0); rho(i) = 0; else; rho(i) = 1/yTs; end;
    endfor

    % Get values of alpha and update q
    for i = m:-1:1
        sTq = innerprod(s(:,:,i), q);
        
        alpha(i) = rho(i) * sTq;
        q = q - alpha(i) * y(:, :, i);
    endfor

    % Get the value of gamma     
    yTy = innerprod(y(:,:,m));
    if yTy == 0
        gamma = 1;
    else
        sTy = innerprod(s(:,:,m), y(:,:,m));

        gamma = sTy / yTy;
    endif

    % Get initial value of search direction
    p = gamma * q;

    % Get values of beta and update p
    for i = 1:m
        yTp = innerprod(y(:,:,i), p);
        
        beta(i) = rho(i) * yTp;
        p = p + s(:, :, i) * (alpha(i) - beta(i));
    endfor
    % Get correct direction
    p = -p;
    
    % Update registers
    x(:, :, 1:m) = x(:, :, 2:m+1);
    g(:, :, 1:m) = g(:, :, 2:m+1);
    
    s(:, :, 1:m-1) = s(:, :, 2:m);
    y(:, :, 1:m-1) = y(:, :, 2:m);
endfunction

function I = innerprod(A, B)
    if nargin == 2
        A_vec = reshape(A, [1, prod(size(A))]);
        B_vec = reshape(B, [1, prod(size(B))]);
    
        I = A_vec(:)' * B_vec(:);
    elseif nargin == 1
        A_vec = reshape(A, [1, prod(size(A))]);
        
        I = A_vec(:)' * A_vec(:);
    else
        error('Error: invalid number of input arguments in innerprod.')
    end
endfunction
