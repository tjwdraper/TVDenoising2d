function dE = energy_gradient_c(c, z, w, Iref, rho)
  % Gradient of c
  Dc = image_gradient(c, 'forward');
  
  % divergence of Dc - z + w
  div_rhs = div(Dc - z + w);

  % Update
  dE = c - Iref + rho * div_rhs;
endfunction
