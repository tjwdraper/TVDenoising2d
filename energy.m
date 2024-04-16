function E = energy(image, Iref, lambda)
    Essd = 1/2 * sum( (image(:) - Iref(:)).^2 );
    
    Dimage = image_gradient(image, 'forward');
    dIdx = Dimage(:,:,1);
    dIdy = Dimage(:,:,2);
    
    
    Etv = lambda * sum( sqrt(dIdx(:).^2 + dIdy(:).^2) );
    
    E = Essd + Etv;
endfunction