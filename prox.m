function z = prox(c, w, lambda, rho)
  % Construct gradient of c
  Dc = image_gradient(c, 'forward');

  % Construct v
  v = Dc + w;

  z = soft_thresholding(v,lambda/rho);
endfunction

function x = soft_thresholding(a, k)
  x = zeros(size(a));

  x(a > k) = a(a > k) - k;
  x(a < -k) = a(a < -k) + k;
endfunction