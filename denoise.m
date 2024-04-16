clc;
clear all;
close all;

addpath(genpath(pwd));

pkg load image


Iref = phantom();
Iref = (Iref - min(Iref(:))) / (max(Iref(:)) - min(Iref(:)));

Iref_noise = imnoise(Iref, 'gaussian', 0, 0.01);

%% ADDMM model for denoising
[dimx, dimy] = size(Iref_noise);

c = Iref_noise;
z = image_gradient(c, 'forward');
w = zeros(dimx, dimy, 2);

niter = 100;
lambda = 1.0;
rho = 100.0;

convergence = 1e-4;

E = zeros(1, niter);

for iter = 1:niter
  % Update c with L-BFGS
  dc = lbfgs_energy_gradient_c(c, z, w, Iref_noise, rho, iter, niter);

  % Line search satisfying Wolfe-conditions
  step = line_search_c(c, z, w, dc, Iref_noise, rho);
  c = c + step * dc;

  % Update z with l1-proximity operator (= soft-thresholding with lambda/rho)
  z = prox(c, w, lambda, rho);

  % Update w with standard ADDMM update
  Dc = image_gradient(c, 'forward');
  w = w + Dc - z;  
    
  % Calculate energy E = Essd + lambda*E_TV and check for convergence 
  E(iter) = energy(c, Iref, lambda);
  printf('Iteration %d - E: %.3f\n', iter, E(iter));
  
  if iter > 1
      if (E(iter-1) - E(iter)) < E(1)*convergence
          fprintf("Convergence reached on iteration %d (step size: %.10f).\n", iter, (E(iter) - E(iter-1))/E(1));
          break;
      endif
  endif
endfor

% Clear persistent variables in L-BFGS update
clear lbfgs_energy_gradient_c

% show overlap
plot_results(Iref, Iref_noise, c);

% Calculate and compare the peak-signal-to-noise ratio (PSNR) before and after 
% TV denoising.
printf("PSNR\n");
printf("Before denoising: %.3f\n", psnr(Iref_noise, Iref));
printf("After denoising: %.3f\n", psnr(c, Iref));