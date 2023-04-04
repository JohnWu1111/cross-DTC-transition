clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

T_max = 10000;
dt = 1e-4;
t = 0:dt:T_max;
nt = length(t);
eta = 0.1;
theta = 0e-2;
epsilon0 = 5;
delta = epsilon0;
omega = 1;

phi = zeros(1,nt);
dphi = zeros(1,nt);
% phi(floor(nx/2),1) = 1;
phi(1) = 5;

for i = 2:nt 
    epsilon = epsilon0 + delta*cos(omega*2*pi*t(i));
%     [y, dy] = Heun_step(phi(i-1), dphi(i-1), dt, eta, theta, epsilon);
    [y, dy] = myFTCS(phi(i-1), dphi(i-1), dt, eta, 0, epsilon);
    phi(i) = y;
    dphi(i) = dy;
end

cut = 10000;
phi_f = abs(fft(phi(floor(nt/2)+1:end)));
phi_f(1) = 0;
phi_f_main = phi_f(1:cut);
dw = 1/(T_max/2);
w_max = 1/dt;
w = 0:dw:w_max;
w_main = w(1:cut);
[maxtab, ~]=peakdet(phi_f_main, 1e3);
[~, peak] = max(maxtab(:,2));
peak_f = w(maxtab(peak,1));
peak_per = maxtab(peak,2)/sum(maxtab(:,2));

fprintf('eta = %g, epsilon0 = %g, peak_f = %g, peak_per = %g .\n',eta,epsilon0,peak_f,peak_per);

toc;

figure;
set(gcf, 'position', [250 70 1000 900]);
titlename = strcat('\eta = ',num2str(eta), ', \epsilon_0 = ', num2str(epsilon0), ', peak_f = ', num2str(peak_f), ', peak_{per} = ', num2str(peak_per));

plot_result = tiledlayout(2,1);
ax1 = nexttile;
plot(ax1,t(floor(99*nt/100):end), phi(floor(99*nt/100):end))
ax2 = nexttile;
plot(ax2,w_main,phi_f_main);
title(plot_result,titlename)
plot_result.TileSpacing = 'compact';


function [y, dy] = myFTCS(phi, dphi, dt, eta, field, epsilon)
% V = (phi.^3 - epsilon*phi)/2;
V = (phi.^2 - epsilon)/2;
V = V.*phi;
y = phi + dt*dphi;
dy = dphi + dt*(field - eta*dphi - V);
end

function [y, dy] = Heun_step(phi, dphi, dt, eta, theta, epsilon)
field = randn*sqrt(2*eta*theta/dt);
[z, dz] = myFTCS(phi, dphi, dt, eta, field, epsilon);
[y, dy] = myFTCS((phi+z)/2, (dphi+dz)/2, dt, eta, field, epsilon);
end