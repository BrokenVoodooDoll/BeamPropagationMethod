clc
clear

%% User choices %%
lambda = 1; % wavelentgh in um
N = 2^9; % number of points
Lz = 3000; % z-limit of visualization in um
Lx = 300; % x-limit of visualization in um
A = 0.8; % amplitude of the beam
r = 30; % radius of the beam in um;
%%%%%

k0 = 2 * pi / lambda; % wavenumber

dz = Lz / (N - 1);
dx = Lx / (N - 1);
dkx = 2 * pi / Lx;

z = 0:dz:Lz;
x = -Lx/2:dx:Lx/2;
kx = -pi/dx:dkx:pi/dx;

E = A * exp(-(x / r).^2); % Gaussian beam
% E = A * rect(x / (2 * r)); % uniform beam

dn = @(E) (10^-3) * abs(E).^2; % change of refractive index vs field

I = zeros(N, N); % resulting intensity
for n = 1:N
    E = exp(-1i * k0 * dn(E) * dz) .* ifft(fftshift(exp(-1i * sqrt(k0^2 - kx.^2) * dz) .* fftshift(fft(E))));
    I(n, :) = abs(E);    
end

imagesc(z, x, rot90(I));
title(sprintf("beam amplitude = %.1f. Beam radius = %.1f mm", A, r))
xlabel("z, mm")
ylabel("x, mm")
colorbar