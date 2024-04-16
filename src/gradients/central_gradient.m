function dI = central_gradient(I)
  dIdx = zeros(size(I));
  dIdy = zeros(size(I));

  % Central grid
  dIdx(2:end-1, :) = 0.5 * (I(3:end, :) - I(1:end-2, :));
  dIdy(:, 2:end-1) = 0.5 * (I(:, 3:end) - I(:, 1:end-2));

  % Boundary
  dIdx(1, :) = 0.5*(-3*I(1, :) + 4*I(2, :) - I(3, :));
  dIdy(:, 1) = 0.5*(-3*I(:, 1) + 4*I(:, 2) - I(:, 3));

  dIdx(end, :) = 0.5*(3*I(end, :) - 4*I(end-1, :) + I(end-2, :));
  dIdy(:, end) = 0.5*(3*I(:, end) - 4*I(:, end-1) + I(:, end-2));

  dI = cat(3, dIdx, dIdy);
endfunction
