function divz = div(z)
  if ndims(z) == 3
    dz1 = image_gradient(squeeze(z(:,:,1)), 'backward');
    dz2 = image_gradient(squeeze(z(:,:,2)), 'backward');

    divz = squeeze(dz1(:,:,1) + dz2(:,:,2));

  else
    error('Invalid dimensions of input z to div function.')
  endif
endfunction
