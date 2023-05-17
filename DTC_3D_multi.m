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
L = 8;
nx = L^3;
N_comp = 2;
eta = 0.05;
theta = 1;
epsilon0 = 5;
omega = 1;

phi = zeros(L,L,2*L,N_comp);
% phi(floor(nx/2),1) = 1;
% phi(:,:,1:L,:) = epsilon0*(2*rand(L,L,L,N_comp)-1);
% phi(:,:,1:L,:) = epsilon0*(1+0.2*randn(L,L,L,N_comp));
% phi(:,:,1:L,1) = 4;
% phi(:,:,1:L,2) = -4;
phi(:,:,1:L,:) = epsilon0*(rand(L,L,L,N_comp));
% phi(:,:,1:L,:) = epsilon0*(randn(L,L,L,N_comp));
order_total = zeros(1,nt);
order = zeros(N_comp,nt);
% order(1) = sum(sqrt(phi(1,:,:,:).^2 + phi(2,:,:,:).^2),'all');
order(:,1) = sum(phi(:,:,1:L,:),1:3);
order_total(1) = sqrt(sum(order(:,1).^2));

t_it = 0;
count = 2;
for i = 2:round(T_max/dt)+1
    t_it = t_it + dt;
%     epsilon = epsilon0 + delta*cos(omega*2*pi*t(i));
    phi = Heun_step(phi, t_it, dt, eta, theta, epsilon0, omega*2*pi, N_comp);
%     phi = myFTCS(phi, t_it, dt, eta, 0, epsilon0, omega*2*pi);
%     order(i) = sum(sqrt(phi(1,:,:,:).^2 + phi(2,:,:,:).^2),'all');
    if mod(i+1,M) == 0
        order(:,count) = sum(phi(:,:,1:L,:),1:3);
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



function y = myFTCS(phi, t, dt, eta, field, epsilon, omega)
c1 = f(phi, eta, field, epsilon*(cos(omega*t)+1));
c2 = f(phi+c1*dt/2, eta, field, epsilon*(cos(omega*(t+dt/2))+1));
c3 = f(phi+c2*dt/2, eta, field, epsilon*(cos(omega*(t+dt/2))+1));
c4 = f(phi+c3*dt, eta, field, epsilon*(cos(omega*(t+dt))+1));
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
end

function y = Heun_step(phi, t, dt, eta, theta, epsilon, omega, N_comp)
L = length(phi)/2;
field = randn(L,L,L,N_comp)*sqrt(2*eta*theta/dt);
z = myFTCS(phi, t, dt, eta, field, epsilon, omega);
y = myFTCS((phi+z)/2, t, dt, eta, field, epsilon, omega);
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