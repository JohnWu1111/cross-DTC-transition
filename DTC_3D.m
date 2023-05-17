clear;
% close all;
clc;
format long
tic;

myseed = 3;
rng(myseed)

T_max = 1000;
dt = 1e-3;
M = 10;
t = 0:dt*M:T_max;
nt = length(t);
L = 20;
nx = L^3;
eta = 0.05;
theta = 35;
epsilon0 = 5;
omega = 0.5;

dphi = zeros(L,L,L);
% phi(:,:,1:L) = epsilon0;
phi = epsilon0*rand(L,L,L);

order = zeros(1,nt);
order(:,1) = sum(phi,"all");

t_it = 0;
count = 1;
for i = 2:round(T_max/dt)+1
    t_it = t_it + dt;
%     epsilon = epsilon0 + delta*cos(omega*2*pi*t(i));
%     [phi, dphi] = Heun_step(phi, dphi, t_it, dt, eta, theta, epsilon0, omega*2*pi);
    [phi, dphi] = myFTCS(phi, dphi, t_it, dt, eta, 0, epsilon0, omega*2*pi);
%     order(i) = sum(sqrt(phi(1,:,:,:).^2 + phi(2,:,:,:).^2),'all');
    if mod(i+1,M) == 0
        order(:,count) = sum(phi,"all");
        count = count + 1;
    end
end
order = order/nx;

cut = 2000;
phi_f = abs(fft(order(:,floor(nt/2)+1:end),(nt+1)/2,2));
% phi_f = abs(fft(order(1,floor(nt/2)+1:end)));
phi_f(:,1) = 0;
phi_f_main = phi_f(:,1:cut);
dw = 1/(T_max/2);
w_max = 1/(M*dt);
w = 0:dw:w_max;
w_main = w(1:cut);
[maxtab, ~]=peakdet(phi_f_main(1,:), max(phi_f_main(1,:))/3);
if ~isempty(maxtab)
    [~, peak] = max(maxtab(:,2));
    peak_f = w(maxtab(peak,1));
    peak_per = maxtab(peak,2)/sum(maxtab(:,2));
else
    peak_f = 0;
    peak_per = 0;
    warning('no peak in the inteval!')
end

fprintf('eta = %g, epsilon0 = %g, peak_f = %g, peak_per = %g .\n',eta,epsilon0,peak_f,peak_per);

toc;

figure;
set(gcf, 'position', [250 70 1500 900]);
titlename = strcat('dt = ',num2str(dt),', \eta = ',num2str(eta),', \theta = ',num2str(theta), ', \epsilon_0 = ', num2str(epsilon0), ', peak_f = ', num2str(peak_f), ', peak_{per} = ', num2str(peak_per));

plot_result = tiledlayout(2,2);
ax1 = nexttile;
plot(ax1,t, order)
ax2 = nexttile;
plot(ax2,t(floor(80*nt/100):end), order(floor(80*nt/100):end));
ax3 = nexttile;
plot(ax3,t,order);
ax4 = nexttile;
plot(ax4,w_main,phi_f_main);
title(plot_result,titlename)
plot_result.TileSpacing = 'compact';

function [y, z] = myFTCS(phi, dphi, t, dt, eta, field, epsilon, omega)
[c1, d1] = f(phi, dphi, eta, field, epsilon*(cos(omega*t)+1));
[c2, d2] = f(phi+c1*dt/2, dphi+d1*dt/2, eta, field, epsilon*(cos(omega*(t+dt/2))+1));
[c3, d3] = f(phi+c2*dt/2, dphi+d2*dt/2, eta, field, epsilon*(cos(omega*(t+dt/2))+1));
[c4, d4] = f(phi+c3*dt, dphi+d3*dt, eta, field, epsilon*(cos(omega*(t+dt))+1));
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
z = dphi + dt*(d1+2*d2+2*d3+d4)/6;
end

function [y1, y2] = Heun_step(phi, dphi, t, dt, eta, theta, epsilon, omega)
L = length(phi);
field = randn(L,L,L)*sqrt(2*eta*theta/dt);
[z1, z2] = myFTCS(phi, dphi, t, dt, eta, field, epsilon, omega);
[y1, y2] = myFTCS((phi+z1)/2, (dphi+z2)/2, t, dt, eta, field, epsilon, omega);
end

function [y, z] = f(x, dx, eta, field, epsilon)

phixx = ((circshift(x,1) + circshift(x,-1) + circshift(x,1,2) ...
    + circshift(x,-1,2) + circshift(x,1,3) + circshift(x,-1,3)) - 6*x);
V = (x.^2 - epsilon)/2;
V = V.*x;
y = dx;
z = field - eta*dx + phixx - V;
end