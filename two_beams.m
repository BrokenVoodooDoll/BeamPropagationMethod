clc
clear

%% User choices %%
lambda = 1; % wavelength in um
k0 = 2 * pi / lambda; % wavenumber
N = 2^10; % number of points
Lz = 1000; % z-limit of visualization in um
Lx = 1200; % x-limit of visualization in um

x0_1 = 200; % x-position of the 1st beam in um
x0_2 = -200; % x-position of the 2nd beam in um

r_1 = 60; % radius of the 1st beam in um
r_2 = 60; % radius of the 2nd beam in um

t_1 = 10; % tilt of the 1st beam in degrees
t_2 = 10; % tilt of the 2nd beam in degrees
%%%%%

dz = Lz / (N - 1);
dx = Lx / (N - 1);
dkx = 2 * pi / Lx;

z = 0:dz:Lz;
x = -Lx/2:dx:Lx/2;
kx = -pi/dx:dkx:pi/dx;

%% Gaussian beam %%
E1 = exp(-((x - x0_1) / r_1).^2 + 1i * k0 * sind(t_1) * x); 
E2 = exp(-((x - x0_2) / r_2).^2 - 1i * k0 * sind(t_2) * x);
%% Uniform beam %%
% E1 = rect((x - x0_1) / (2 * r_1)) .* exp(1i * k0 * sind(t_1) * x); 
% E2 = rect((x - x0_2) / (2 * r_2)) .* exp(-1i * k0 * sind(t_2) * x);
%%%%%
E = E1 + E2;

I = zeros(N, N); % resulting intensity
for n = 1:N
    E = exp(-1i * k0 * dz) .* ifft(fftshift(exp(-1i * sqrt(k0^2 - kx.^2) * dz) .* fftshift(fft(E))));
    I(n, :) = abs(E);    
end   

imagesc(z, x, rot90(I));
s_1 = "Two crossing beams";
s_2 = sprintf("First beam: {x_0} = %.1f {\\mu}m; radius = %.1f {\\mu}m; tilt = %.1f{\\deg}", x0_1, r_1, t_1);
s_3 = sprintf("Second beam: {x_0} = %.1f {\\mu}m; radius = %.1f {\\mu}m; tilt = %.1f{\\deg}", x0_2, r_2, t_2);
title({s_1; s_2; s_3})
xlabel("z, {\\mu}m")
ylabel("x, {\\mu}m")
colorbar
% colormap(gray)