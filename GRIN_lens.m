clc
clear

%% User choices %%
lambda = 1; % wavelength in um
N = 2^10; % number of points
Lz = 600; % z-limit of visualization in um
Lx = 600; % x-limit of visualization in um
r = 50; % radius of the beam in um

l_start = 150; % z-coordinate of the input surface of the lens
l_end = 250; % z-coordinate of the output surface of the lens
p = 0.01; % parameter of the lens
%%%%%

k0 = 2 * pi / lambda; % wavenumber

dz = Lz / (N - 1);
dx = Lx / (N - 1);
dkx = 2 * pi / Lx;

z = 0:dz:Lz;
x = -Lx/2:dx:Lx/2;
kx = -pi/dx:dkx:pi/dx;

E = exp(-(x / r).^2); % Gaussian beam
% E = rect(x / (2 * r)); % uniform beam

l_width = l_end - l_start; % width of the gradient lens
dn = @(x, z) -0.5 * p^2 * x.^2 * rect((z - l_start - (l_width / 2)) / l_width); % function of the GRIN lens

I = zeros(N, N); % resulting intensity
for n = 1:N
    E = exp(-1i * k0 * dn(x, dz * n) * dz) .* ifft(fftshift(exp(-1i * sqrt(k0^2 - kx.^2) * dz) .* fftshift(fft(E))));
    I(n, :) = abs(E);    
end

figure
imagesc(z, x, rot90(I))
colorbar
title(sprintf("Beam radius = %.1f {\\mu}m. Lens start = %.1f {\\mu}m. Lens end = %.1f {\\mu}m. Lens parameter = %.3f", r, l_start, l_end, p))
xlabel("z, mm")
ylabel("x, mm")

line([l_start l_start], [-Lx/2 Lx/2], 'Color','red','LineStyle','--');
line([l_end l_end], [-Lx/2 Lx/2], 'Color','red','LineStyle','--');




% colormap(gray)