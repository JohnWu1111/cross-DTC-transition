clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

T_max = 1000;
dt = 1e-3;
t = 0:dt:T_max;
nt = length(t);
L = 8;
nx = L^3;
N_comp = 1;
eta = 0.05;
theta = 2e-2;
epsilon0 = 5;
omega = 1;

phi = zeros(L,L,2*L,N_comp);
% phi(floor(nx/2),1) = 1;
% phi(:,:,1:L,:) = epsilon0*(2*rand(L,L,L,N_comp)-1);
phi(:,:,1:L,:) = epsilon0;
order_total = zeros(1,nt);
order = zeros(N_comp,nt);
% order(1) = sum(sqrt(phi(1,:,:,:).^2 + phi(2,:,:,:).^2),'all');
order_total(1) = sum(phi,'all');
order(:,1) = sum(phi,1:3);

for i = 2:nt 
    phi = Heun_step(phi, t(i), dt, eta, theta, epsilon0, omega*2*pi, N_comp);
%     phi = myFTCS(phi, t(i), dt, eta, 0, epsilon0, omega*2*pi);
%     order(i) = sum(sqrt(phi(1,:,:,:).^2 + phi(2,:,:,:).^2),'all');
    order_total(i) = sum(phi,'all');
    order(:,i) = sum(phi,1:3);
end
order_total = order_total/nx;
order = order/nx;

cut = 1000;
phi_f = abs(fft(order_total(floor(nt/2)+1:end)));
phi_f(1) = 0;
phi_f_main = phi_f(1:cut);
dw = 1/(T_max/2);
w_max = 1/dt;
w = 0:dw:w_max;
w_main = w(1:cut);
[maxtab, ~]=peakdet(phi_f_main, max(phi_f_main)/3);
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
titlename = strcat('dt = ',num2str(dt),', \eta = ',num2str(eta), ', \epsilon_0 = ', num2str(epsilon0), ', peak_f = ', num2str(peak_f), ', peak_{per} = ', num2str(peak_per));

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



% function y = myFTCS(phi, t, dt, eta, field, epsilon0, omega)
% c1 = f(phi, eta, field, epsilon0*(cos(omega*t)+1));
% c2 = f(phi+c1*dt/2, eta, field, epsilon0*(1+cos(omega*(t+dt/2))));
% c3 = f(phi+c2*dt/2, eta, field, epsilon0*(1+cos(omega*(t+dt/2))));
% c4 = f(phi+c3*dt, eta, field, epsilon0*(1+cos(omega*(t+dt))));
% y = phi + dt*(c1+2*c2+2*c3+c4)/6;
% end

function y = myFTCS(phi, t, dt, eta, field, epsilon0, omega)
fact_dt = [1/5 3/10 4/5 8/9 1];
fact_para = [1/5 0 0 0 0;
    3/40 9/40 0 0 0;
    44/45 -56/15 32/9 0 0;
    19372/6561 -25360/2187 64448/6561 -212/729 0;
    9017/3168 -355/33 46732/5247 49/176 -5103/18656];
fact_sum = [35/384 0 500/1113 125/192 -2187/6784 11/84];
c1 = f(phi, eta, field, epsilon0*(cos(omega*t)+1));
c2 = f(phi+dt*c1/5, eta, field, epsilon0*(1+cos(omega*(t+dt/5))));
c3 = f(phi+dt*(c1*3/40+c2*9/40), eta, field, epsilon0*(1+cos(omega*(t+dt*3/10))));
c4 = f(phi+dt*(c1*44/45-c2*56/15+c3*32/9), eta, field, epsilon0*(1+cos(omega*(t+dt*4/5))));
c5 = f(phi+dt*(c1*19372/6561-c2*25360/2187+c3*64448/6561-c4*212/729), eta, field, epsilon0*(1+cos(omega*(t+dt*8/9))));
c6 = f(phi+dt*(c1*9017/3168-c2*355/33+c3*46732/5247+c4*49/176-c5*5103/18656), eta, field, epsilon0*(1+cos(omega*(t+dt))));
y = phi + dt*(c1*35/384+c3*500/1113+c4*125/192-c5*2187/6784+c6*11/84);
end

function y = Heun_step(phi, t, dt, eta, theta, epsilon0, omega, N_comp)
L = length(phi)/2;
field = randn(L,L,L,N_comp)*sqrt(2*eta*theta/dt);
z = myFTCS(phi, t, dt, eta, field, epsilon0, omega);
y = myFTCS((phi+z)/2, t, dt, eta, field, epsilon0, omega);
end

function y = f(x, eta, field, epsilon)
L = length(x)/2;
y = x;
phi = x(:,:,1:L,:);
% dphi = x(:,:,:,L+1:end);

phixx = ((circshift(phi,1) + circshift(phi,-1) + circshift(phi,1,2) ...
    + circshift(phi,-1,2) + circshift(phi,1,3) + circshift(phi,-1,3)) - 6*phi);
% V = (phi.^3 - epsilon*phi)/2;
% phi_norm = phi(1,:,:,:).^2 + phi(2,:,:,:).^2;
phi_norm = sum(phi.^2,4);
V = (phi_norm - epsilon)/2;
V = V.*phi;
% y(:,:,:,1:L) = dphi;
y = circshift(y,L,3);
y(:,:,L+1:end,:) = field - eta*y(:,:,1:L,:) + phixx - V;
end