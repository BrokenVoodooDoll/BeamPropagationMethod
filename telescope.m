clc
clear

%% User choices %%
x_max = 2; % maximum x in mm
f1 = 100; % focal length of the 1st lens of the telescope in mm
f2 = 100; % focal length of the 2nd lens of the telescope in mm
r = 0.5; % radius of the beam in mm
lambda = 1e-3; % wavelength in mm
N = 2^10; % number of points
%%%%%

k0 = 2 * pi / lambda; % wavenumber
Lz = 2 * (f1 + f2); % z-limit of visualization in mm
Lx = 2 * x_max; % x-limit of visualization in mm

dz = Lz / (N - 1);
dx = Lx / (N - 1);
dkx = 2 * pi / Lx;

z = 0:dz:Lz; % z coordinate array
x = -Lx/2:dx:Lx/2; % x coordinate array
kx = -pi/dx:dkx:pi/dx; % kx coordinate array

tau1 = exp(1i * k0 * x.^2 * 0.5 / f1); % 1st lens phase function
tau2 = exp(-1i * k0 * x.^2 * 0.5 / f2); % 2nd lens phase function

E = exp(-(x / r).^2); % Gaussian beam
% E = rect(x / (2 * r)); % Uniform beam

t = 0; % iterator for z
n = 1; % counter
I = zeros(N, N); % resulting intensity 

while t <= f1
    E = exp(-1i * k0 * dz) * ifft(fftshift(exp(-1i * sqrt(k0^2 - kx.^2) * dz) .* fftshift(fft(E))));
    I(n, :) = abs(E); 
    t = t + dz;
    n = n + 1;
end

E = E .* tau1;

while t <= 2 * f1 + f2
    E = exp(-1i * k0 * dz) * ifft(fftshift(exp(-1i * sqrt(k0^2 - kx.^2) * dz) .* fftshift(fft(E))));
    I(n, :) = abs(E); 
    t = t + dz;
    n = n + 1;
end

E = E.*tau1;

while t <= 2 * (f1 + f2)
    E = exp(-1i * k0 * dz) * ifft(fftshift(exp(-1i * sqrt(k0^2 - kx.^2) * dz) .* fftshift(fft(E))));
    I(n, :) = abs(E); 
    t = t + dz;
    n = n + 1;
end

figure
imagesc(z, x, rot90(I))
line([f1 f1], [-x_max x_max], 'Color','red','LineStyle','--');
line([2*f1+f2 2*f1+f2], [-x_max x_max], 'Color','red','LineStyle','--');

colorbar
title(sprintf("Kepler telescope with {f_1 =} %.2f mm and {f_2} = %.2f mm", f1, f2))
xlabel("z, mm")
ylabel("x, mm")
% colormap(gray)
