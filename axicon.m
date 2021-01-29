clc
clear

%% User choices %%
lambda = 1; %% wavelength in um
N = 2^10; % number of points
Lz = 200; % z-limit of visualization in um
Lx = 600; % x-limit of visualization in um
r = 100; % radius of the beam in um
n0 = 1.5; % refractive index of the axicon
a = 140; % apex angle of the axicon in degrees
%%%%%

k0 = 2*pi/lambda; % wavenumber

dz = Lz / (N - 1);
dx = Lx / (N - 1);
dkx = 2 * pi / Lx;

z = 0:dz:Lz;
x = -Lx/2:dx:Lx/2;
kx = -pi/dx:dkx:pi/dx;

%% Axicon parameters %%
g = 90 - a / 2;
b = asind(n0 * cosd(a / 2)) + a / 2 - 90;
p = sind(b);
%%%%%

tau = exp(1i * k0 * p * abs(x)); % phase function of the axicon
E = exp(-x.^2 / r^2) .* tau; % Gaussian beam
% E = rect(x / (2 * r)) .* tau; % uniform beam

I = zeros(N, N); % resulting intensity
for n = 1:N
    E = exp(-1i * k0 * dz) .* ifft(fftshift(exp(-1i * sqrt(k0^2 - kx.^2) * dz) .* fftshift(fft(E))));
    I(n, :) = abs(E).^2;
end

imagesc(z, x, rot90(I));
title(sprintf("Beam radius = %.1f {\\mu}m. Apex angle = %.1f{\\deg}. Refractive index = %.1f", r, a, n0))
xlabel("z, mm")
ylabel("x, mm")
colorbar
%colormap(gray)
