clc
clear

%% User choices %%
x_max = 2; % % maximum x in mm
r = 2; % radius of the beam in mm
f = 100; % focal length of the lens
lambda = 1e-3; % wavelemgth in mm
N = 2^10; % number of points
%%%%%

k0 = 2 * pi / lambda; % wavenumber
Lz = 2 * f; % z-limit of visualization in mm
Lx = 2 * x_max; % x-limit of visualization in mm

dz = Lz / (N - 1);
dx = Lx / (N - 1);
dkx = 2 * pi / Lx;

z = 0:dz:Lz;
x = -Lx/2:dx:Lx/2;
kx = -pi/dx:dkx:pi/dx;

tau = exp(1i * k0 * x.^2 * 0.5 / f); % lens phase function

E = exp(-(x / r).^2) .* tau; % Gaussian beam
% E = rect(x / (2 * r)) .* tau; % Uniform beam

I = zeros(N, N); % resulting intensity

for n = 1:N
    E = ifft(fftshift(exp(-1i * sqrt(k0^2 - kx.^2) * dz) .* fftshift(fft(E))));
    I(n, :) = abs(E);    
end

imagesc(z, x, rot90(I));
title(sprintf("Focal length = %.1f mm. Beam radius = %.1f mm", f, r))
xlabel("z, mm")
ylabel("x, mm")
colorbar
% colormap(gray)
