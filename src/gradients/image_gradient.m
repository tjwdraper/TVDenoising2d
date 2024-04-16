function dI = image_gradient(I, opt)
    if nargin == 1
      dI = central_gradient(I);
    else
      if strcmp(opt, 'central')
        dI = central_gradient(I);
      elseif strcmp(opt, 'forward')
        dI = forward_gradient(I);
      elseif strcmp(opt, 'backward')
        dI = backward_gradient(I);
      else
        error('invalid option argument in image_gradient');
      endif
    endif
endfunction
