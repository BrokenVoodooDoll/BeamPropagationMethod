clc
clear

%% User choices %%
lambda = 1; % wavelength in um
N = 2^9; % number of points
Lz = 3000; % z-limit of visualization in um
Lx = 1500; % x-limit of visualization in um
r = 240; % radius of the beam in um

phi = 7 % incident angle of the beam in degrees
T = lambda / (2 * sind(phi)); % period of the grating in um
% T = 5 * lambda;

g_start = 500; % z-coordinate of the input surface of the grating
g_end = 1000; % z-coordinate of the output surface of the grating
%%%%%
k0 = 2*pi/lambda; % wavenumber

dz = Lz / (N - 1);
dx = Lx / (N - 1);
dkx = 2 * pi / Lx;

z = 0:dz:Lz;
x = -Lx/2:dx:Lx/2;
kx = -pi/dx:dkx:pi/dx;

tau = exp(-1i * k0 * sind(phi) * x); % phase function of the grating
E = exp(-(x / r).^2) .* tau;
%E = rect(x / (2*r)) .* tau; %uniform beam

g_width = g_end - g_start; % width of the gradient grating
dn = @(x, z) 0.001 * rect((z - g_start - (g_width / 2)) / g_width) * cos(2 * pi / T * x);

I = zeros(N, N); % resulting intensity
for n = 1:N
    E = exp(-1i * k0 * dn(x, dz*n) * dz) .* ifft(fftshift(exp(-1i * sqrt(k0^2 - kx.^2) * dz) .* fftshift(fft(E))));
    I(n, :) = abs(E);    
end   

figure
imagesc(z, x, rot90(I))
colorbar
title(sprintf("Beam radius = %.1f {\\mu}m. Incident angle = %.1f{\\deg} Grating start = %.1f {\\mu}m. Grating end = %.1f {\\mu}m. Grating period = %.3f {\\mu}m", r, phi, g_start, g_end, T))
xlabel("z, mm")
ylabel("x, mm")

line([g_start g_start], [-Lx/2 Lx/2], 'Color','red','LineStyle','--');
line([g_end g_end], [-Lx/2 Lx/2], 'Color','red','LineStyle','--');