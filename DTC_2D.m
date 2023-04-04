clear;
% close all;
clc;
format long
tic;

myseed = 2;
rng(myseed)

T_max = 1000;
dt = 1e-4;
t = 0:dt:T_max;
nt = length(t);
L = 10;
eta = 1;
theta = 0e-2;
epsilon0 = 5;
delta = 5;
omega = 1;

dphi = zeros(L);
% phi(floor(nx/2),1) = 1;
phi = ones(L);
order = zeros(1,nt);
order(1) = sum(phi.^2,"all");
order2 = zeros(1,nt);
order2(1) = sum(phi.^2,"all");
phit = zeros(1,nt);
phit(1) = phi(1);

for i = 2:nt 
    epsilon = epsilon0 + delta*cos(omega*2*pi*t(i));
    [y, dy] = Heun_step(phi, dphi, dt, eta, theta, epsilon);
%     [y, dy] = myFTCS(phi, dphi, dt, eta, 0, epsilon);
    phi = y;
    dphi = dy;
    order(i) = sum(phi,"all");
    order2(i) = sum(phi.^2,"all");
    phit(i) = phi(1);
end
order = order/L^2;
order2 = sqrt(order2)/L;

% cut = 10000;
% dw = 1/(T_max/2);
% w_max = 1/dt;
% w = 0:dw:w_max;
% w_main = w(1:cut);
% 
% order_f = abs(fft(order(floor(nt/2)+1:end)));
% order_f(1) = 0;
% order_f_main = order_f(1:cut);
% 
% [maxtab, ~]=peakdet(order_f_main, 1e3);
% [~, peak] = max(maxtab(:,2));
% peak_f = w(maxtab(peak,1));
% peak_per = maxtab(peak,2)/sum(maxtab(:,2));
% 
% fprintf('eta = %g, epsilon0 = %g, peak_f = %g, peak_per = %g .\n',eta,epsilon0,peak_f,peak_per);
% 
% toc;
% 
% titlename = strcat('L = ',num2str(L),', \eta = ',num2str(eta), ', \epsilon_0 = ', num2str(epsilon0), ', peak_f = ', num2str(peak_f), ', peak_{per} = ', num2str(peak_per));
% figure('Name',titlename);
% set(gcf, 'position', [250 70 1000 900]);
% 
% subplot(2,2,1)
% plot(t(floor(95*nt/100):end), order(floor(95*nt/100):end))
% xlabel('t')
% ylabel('order')
% 
% subplot(2,2,2)
% plot(w_main,order_f_main);
% xlabel('\omega / 2\pi')
% 
% subplot(2,2,3)
% plot(t(floor(95*nt/100):end), order2(floor(95*nt/100):end));
% xlabel('t')
% ylabel('squareroot order')
% 
% subplot(2,2,4)
% plot(t(floor(95*nt/100):end), phit(floor(95*nt/100):end));
% xlabel('t')
% ylabel('one point')

function [y, dy] = myFTCS(phi, dphi, dt, eta, field, epsilon)
phi_u = circshift(phi,1,1);
phi_d = circshift(phi,-1,1);
phi_l = circshift(phi,1,2);
phi_r = circshift(phi,-1,2);
phixx = phi_u + phi_d + phi_l + phi_r - 4*phi;
% V = (phi.^3 - epsilon*phi)/2;
V = (phi.^2 - epsilon)/2;
V = V.*phi;
y = phi + dt*dphi;
dy = dphi + dt*(field - eta*dphi + phixx - V);
end

function [y, dy] = Heun_step(phi, dphi, dt, eta, theta, epsilon)
L = length(phi);
field = randn(L)*sqrt(2*eta*theta/dt);
[z, dz] = myFTCS(phi, dphi, dt, eta, field, epsilon);
[y, dy] = myFTCS((phi+z)/2, (dphi+dz)/2, dt, eta, field, epsilon);
end