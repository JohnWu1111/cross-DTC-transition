clear;
% close all;
clc;
format long
tic;

myseed = 5;
rng(myseed)

T_max = 1000;
dt = 1e-3;
M = 10;
t = 0:dt*M:T_max;
nt = length(t);
L = 12;
nx = L^3;
N_comp = 2;
eta = 0.05;
theta = 5;
epsilon0 = 5;
omega = 0.75;

dphi = zeros(L,L,L,N_comp);
% phi = epsilon0*rand(L,L,L,N_comp);
% phi = epsilon0*(1+0.1*rand(L,L,L,N_comp));
% phi = epsilon0*ones(L,L,L,N_comp);
phi = epsilon0*randn(L,L,L,N_comp);

order_total = zeros(1,nt);
order = zeros(N_comp,nt);
% order(1) = sum(sqrt(phi(1,:,:,:).^2 + phi(2,:,:,:).^2),'all');
order(:,1) = sum(phi,1:3);
order_total(1) = sqrt(sum(order(:,1).^2));

t_it = 0;
count = 2;
for i = 2:round(T_max/dt)+1
    t_it = t_it + dt;
%     epsilon = epsilon0 + delta*cos(omega*2*pi*t(i));
%     [phi, dphi] = Heun_step(phi, dphi, t_it, dt, eta, theta, epsilon0, omega*2*pi, N_comp);
    [phi, dphi] = myFTCS(phi, dphi, t_it, dt, eta, 0, epsilon0, omega*2*pi);
%     order(i) = sum(sqrt(phi(1,:,:,:).^2 + phi(2,:,:,:).^2),'all');
    if mod(i+1,M) == 0
        order(:,count) = sum(phi,1:3);
        order_total(count) = sqrt(sum(order(:,count).^2));
        count = count + 1;
    end
end
order_total = order_total/nx;
order = order/nx;

% cut = round(nt/3);
cut = 1000;
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
plot(ax1,t, order_total)
ax2 = nexttile;
plot(ax2,t(floor(80*nt/100):end), order_total(floor(80*nt/100):end));
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

function [y1, y2] = Heun_step(phi, dphi, t, dt, eta, theta, epsilon, omega, N_comp)
L = length(phi);
field = randn(L,L,L, N_comp)*sqrt(2*eta*theta/dt);
[z1, z2] = myFTCS(phi, dphi, t, dt, eta, field, epsilon, omega);
[y1, y2] = myFTCS((phi+z1)/2, (dphi+z2)/2, t, dt, eta, field, epsilon, omega);
end

function [y, z] = f(x, dx, eta, field, epsilon)

phi_l = circshift(x,1,1);
phi_r = circshift(x,-1,1);
phi_u = circshift(x,1,2);
phi_d = circshift(x,-1,2);
phi_a = circshift(x,1,3);
phi_b = circshift(x,-1,3);

phixx = phi_l + phi_r + phi_u + phi_d + phi_a + phi_b - 6*x;
phi_norm = sum(x.^2,4);
V = (phi_norm - epsilon)/2;
V = V.*x;
y = dx;
z = field - eta*dx + phixx - V;
end