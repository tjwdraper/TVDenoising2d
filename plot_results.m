function plot_results(Iref, Iref_noise, c)
    figure();
    subplot(331); imagesc(Iref); title('Iref');
    subplot(332); imagesc(Iref_noise); title('Iref noisy');
    subplot(333); imagesc(c); title('Iref TV-regularized');
    subplot(335); imagesc(Iref_noise - Iref); title('Iref noisy - Iref'); caxis([-1/2 1/2]);
    subplot(336); imagesc(c - Iref); title('Iref TV - Iref'); caxis([-1/2 1/2]);
    subplot(337); imagesc(l1_image(Iref)); 
    subplot(338); imagesc(l1_image(Iref_noise));
    subplot(339); imagesc(l1_image(c))
    colormap gray;
endfunction

function imageout = l1_image(image)
    dI = image_gradient(image, 'forward');
    
    imageout = sqrt(dI(:,:,1).^2 + dI(:,:,2).^2);
endfunction