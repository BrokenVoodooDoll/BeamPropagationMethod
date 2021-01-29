clc
clear

%% User choices %%
lambda = 1; % wavelength in um
N = 2^10; % number of points
Lz = 2000; % z-limit of visualization in um
Lx = 600; % x-limit of visualization in um
r = 60; % radius of the beam in um

g_start = 500; % z-coordinate of the input surface of the glass
g_end = 1500; % z-coordinate of the output surface of the glass
p = 0.07; % parameter of the glass
%%%%%

k0 = 2 * pi / lambda; % wavenumber

dz = Lz / (N - 1);
dx = Lx / (N - 1);
dkx = 2 * pi / Lx;

z = 0:dz:Lz;
x = -Lx/2:dx:Lx/2;
kx = -pi/dx:dkx:pi/dx;

E = exp(-(x / r).^2); % Gaussian beam
%E = rect(x / (2 * r)); % uniform beam

g_width = g_end - g_start; % width of the gradient glass
dn = @(x, z) p * rect((z - g_start - (g_width / 2)) / g_width) * (x / Lx); % function of the gradient glass

I = zeros(N, N); % resulting intensity
for n = 1:N
    E = exp(-1i * k0 * dn(x, dz*n) * dz) .* ifft(fftshift(exp(-1i * sqrt(k0^2 - kx.^2) * dz) .* fftshift(fft(E))));
    I(n, :) = abs(E);    
end   

figure
imagesc(z, x, rot90(I))
colorbar
title(sprintf("Beam radius = %.1f {\\mu}m. Glass start = %.1f {\\mu}m. Glass end = %.1f {\\mu}m. Medium parameter = %.3f", r, g_start, g_end, p))
xlabel("z, mm")
ylabel("x, mm")

line([g_start g_start], [-Lx/2 Lx/2], 'Color','red','LineStyle','--');
line([g_end g_end], [-Lx/2 Lx/2], 'Color','red','LineStyle','--');

% colormap(gray)