function dI = backward_gradient(I)
    % Allocate
    dIdx = zeros(size(I));
    dIdy = zeros(size(I));

    dIdx(2:end, :) = I(2:end, :) - I(1:end-1, :);
    dIdy(:, 2:end) = I(:, 2:end) - I(:, 1:end-1);

    % Boundary
    dIdx(1, :) = I(2, :) - I(1, :);
    dIdy(:, 1) = I(:, 2) - I(:, 1);

    dI = cat(3, dIdx, dIdy);

endfunction
