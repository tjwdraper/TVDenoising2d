function [E, Essd, Elambda, Erho, Ew] = augmented_lagrangian(c, z, w, Iref, lambda, rho)
  % SSD contribution
  Essd = 1/2 * sum( (Iref(:) - c(:)).^2 );

  % L1 contribution
  dcdx = z(:,:,1);
  dcdy = z(:,:,2);
  Elambda = lambda * sum(sqrt( dcdx(:).^2 + dcdy(:).^2 ));

  % rho contribution
  Dc = image_gradient(c, 'forward');
  Erho = rho/2 * sum( (Dc(:) - z(:) + w(:)).^2 ) ;
  
  % w contribution
  Ew = -rho/2 * sum( w(:).^2 );

  % Total energy
  E = Essd + Elambda + Erho + Ew;

endfunction
