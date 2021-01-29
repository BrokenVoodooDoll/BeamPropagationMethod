clc
clear

%% User choices %%
lambda = 1.064e-3; % wavelength in mm
f1 = 80; % focal length of the 1st lens in mm
f2 = 80; % focal length of the 2nd lens in mm
N = 2^10; % number of points
Lz = 2 * (f1 + f2); % z-limit of visualization in mm
Lx = 10; %  x-limit of visualization in mm

r = 1.5; % radius of the beam in mm
n0 = 1.5; % refractive index of the axicon
a = 176; % axicon apex angle in degrees
%%%%%

k0 = 2 * pi / lambda; % wavenumber

dz = Lz / (N - 1);
dx = Lx / (N - 1);
dkx = 2 * pi / Lx;

z = 0:dz:Lz;
x = -Lx/2:dx:Lx/2;
kx = -pi/dx:dkx:pi/dx;

%% Axicon parameters %%
g = 90 - a/2;
b = asind(n0 * cosd(a / 2)) + a / 2 - 90;
p = sind(b);
%%%%%

tau = exp(1i * k0 * p * abs(x)); % phase function of the axicon
tau1 = exp(1i * k0 * x.^2 * 0.5 / f1); % phase function of the 1st lens
tau2 = exp(1i * k0 * x.^2 * 0.5 / f2); % phase function of the 2nd lens

E = exp(-(x / r).^2) .* tau; % Gaussian beam
%E = rect(x / (2 * r)) .* tau; % uniform beam

n = 1; % counter
t = 0; % current z-coordinate

I = zeros(N, N); % resulting intensity
while t <= f1
    E = exp(-1i * k0 * dz) .* ifft(fftshift(exp(-1i * sqrt(k0^2 - kx.^2) * dz) .* fftshift(fft(E))));
    I(n, :) = abs(E).^2;
    t = t + dz;
    n = n + 1;
end

E = E .* tau1;

while t <= 2 * f1 + f2
    E = exp(-1i * k0 * dz) .* ifft(fftshift(exp(-1i * sqrt(k0^2 - kx.^2) * dz) .* fftshift(fft(E))));
    I(n, :) = abs(E).^2;
    t = t + dz;
    n = n + 1;
end

E = E.*tau2;
while t <= 2 * (f1 + f2)
    E = ifft(fftshift(exp(-1i * sqrt(k0^2 - kx.^2) * dz) .* fftshift(fft(E))));
    I(n, :) = abs(E).^2;
    t = t + dz;
    n = n + 1;
end

imagesc(z, x, rot90(I), [0 4]);
colorbar
s_1 = sprintf("Kepler telescope with {f_1 =} %.2f mm and {f_2} = %.2f mm", f1, f2);
s_2 = sprintf("Beam radius = %.1f mm. Apex angle = %.1f{\\deg}. Refractive index = %.1f", r, a, n0);
title({s_1; s_2})
xlabel("z, mm")
ylabel("x, mm")
line([f1 f1], [-Lx/2 Lx/2], 'Color','red','LineStyle','--');
line([2*f1+f2 2*f1+f2], [-Lx/2 Lx/2], 'Color','red','LineStyle','--');
