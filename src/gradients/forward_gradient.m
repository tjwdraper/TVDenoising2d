function dI = forward_gradient(I)
    % Allocate
    dIdx = zeros(size(I));
    dIdy = zeros(size(I));

    dIdx(1:end-1, :) = I(2:end, :) - I(1:end-1, :);
    dIdy(:, 1:end-1) = I(:, 2:end) - I(:, 1:end-1);

    % Boundary
    dIdx(end, :) = I(end, :) - I(end-1, :);
    dIdy(:, end) = I(:, end) - I(:, end-1);

    dI = cat(3, dIdx, dIdy);
endfunction
