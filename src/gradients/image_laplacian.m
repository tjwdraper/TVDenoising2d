function delI = image_laplacian(I)
  % Approximate second order finite difference approximations
  dxx = zeros(size(I));
  dyy = zeros(size(I));

  % Central finite difference approximations
  dxx(2:end-1, :) = I(3:end, :) - 2*I(2:end-1, :) + I(1:end-2, :);
  dyy(:, 2:end-1) = I(:, 3:end) - 2*I(:, 2:end-1) + I(:, 1:end-2);

  % Boundary
  dxx(1, :) = 2*I(1, :) -5*I(2, :) +4*I(3, :) -I(4, :);
  dyy(:, 1) = 2*I(:, 1) -5*I(:, 2) +4*I(:, 3) -I(:, 4);

  dxx(end, :) = 2*I(end, :) -5*I(end-1, :) + 4*I(end-2, :) - I(end-3, :);
  dyy(:, end) = 2*I(:, end) -5*I(:, end-1) + 4*I(:, end-2) - I(:, end-3);

  % Add the second order finite difference approximations of the partial derivatives.
  delI = dxx + dyy;
endfunction
