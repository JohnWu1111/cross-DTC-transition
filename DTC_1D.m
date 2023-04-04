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
dx = 1;
L = 100;
x = dx:dx:L;
nx = length(x);
eta = 0.01;
theta = 0e-2;
epsilon0 = 10;
delta = 10;
omega = 1;
lambda = dt/dx^2;

dphi = 0;
% phi(floor(nx/2),1) = 1;
phi = rand(nx,1);
order = zeros(1,nt);
order(1) = sum(phi(:,1));

for i = 2:nt 
    epsilon = epsilon0 + delta*cos(omega*2*pi*t(i));
%     [y, dy] = Heun_step(phi, dphi, dx, dt, eta, theta, epsilon);
    [y, dy] = myFTCS(phi, dphi, dx, dt, eta, 0, epsilon);
    phi = y;
    dphi = dy;
    order(i) = sum(y);
end
order = order/L;

cut = 10000;
phi_f = abs(fft(order(floor(nt/2)+1:end)));
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
plot(ax1,t(floor(99*nt/100):end), order(floor(99*nt/100):end))
ax2 = nexttile;
plot(ax2,w_main,phi_f_main);
title(plot_result,titlename)
plot_result.TileSpacing = 'compact';



function [y, dy] = myFTCS(phi, dphi, dx, dt, eta, field, epsilon)
phixx = ((circshift(phi,1) + circshift(phi,-1)) - 2*phi) / dx^2;
% V = (phi.^3 - epsilon*phi)/2;
V = (phi.^2 - epsilon)/2;
V = V.*phi;
y = phi + dt*dphi;
dy = dphi + dt*(field - eta*dphi + phixx - V);
end

function [y, dy] = Heun_step(phi, dphi, dx, dt, eta, theta, epsilon)
nx = length(phi);
field = randn(nx,1)*sqrt(2*eta*theta/dt);
[z, dz] = myFTCS(phi, dphi, dx, dt, eta, field, epsilon);
[y, dy] = myFTCS((phi+z)/2, (dphi+dz)/2, dx, dt, eta, field, epsilon);
end