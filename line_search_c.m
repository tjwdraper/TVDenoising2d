function alpha = line_search_c(c, z, w, p, Iref, rho)
    alpha = 1.0;
    gamma = 0.80;
    c1 = 1e-4;
    
    % Gradient of current and next c
    Dc_curr = image_gradient(c, 'forward');
    Dc_next = image_gradient(c+alpha*p, 'forward');
    
    Ecurr = 1/2 * sum( (Iref(:) - c(:)).^2 ) + rho/2 * sum( (Dc_curr(:) - z(:) + w(:) ).^2 );
    Enext = 1/2 * sum( (Iref(:) - c(:) - alpha*p(:)).^2 ) + rho/2 * sum( (Dc_next(:) - z(:) + w(:) ).^2 ) - rho/2 * sum( w(:).^2 );
    
    while (Enext > Ecurr + c1 * alpha * p(:)' * p(:) )
        alpha = gamma * alpha;
        
        Dc_next = image_gradient(c+alpha*p, 'forward');
        
        Enext = 1/2 * sum( (Iref(:) - c(:) - alpha*p(:)).^2 ) + rho/2 * sum( (Dc_next(:) - z(:) + w(:) ).^2 ) - rho/2 * sum( w(:).^2 );
        
    endwhile
endfunction