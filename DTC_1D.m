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
L = 1000;
x = 1:L;
nx = length(x);
eta = 0.1;
theta = 0e-2;
epsilon0 = 10;
delta = 10;
omega = 1;

phi = zeros(L,2);
% phi(floor(nx/2),1) = 1;
phi(:,1) = rand(nx,1);
order = zeros(1,nt);
order(1) = sum(phi(:,1));

for i = 2:nt 
    epsilon = epsilon0 + delta*cos(omega*2*pi*t(i));
%     phi = Heun_step(phi, dt, eta, theta, epsilon);
    phi = myFTCS(phi, dt, eta, 0, epsilon);
    order(i) = sum(phi(:,1));
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
set(gcf, 'position', [250 70 1000 900]);
titlename = strcat('\eta = ',num2str(eta), ', \epsilon_0 = ', num2str(epsilon0), ', peak_f = ', num2str(peak_f), ', peak_{per} = ', num2str(peak_per));

plot_result = tiledlayout(2,1);
ax1 = nexttile;
plot(ax1,t(floor(95*nt/100):end), order(floor(95*nt/100):end))
ax2 = nexttile;
plot(ax2,w_main,phi_f_main);
title(plot_result,titlename)
plot_result.TileSpacing = 'compact';



function y = myFTCS(phi, dt, eta, field, epsilon)
c1 = f(phi, eta, field, epsilon);
c2 = f(phi+c1*dt/2, eta, field, epsilon);
c3 = f(phi+c2*dt/2, eta, field, epsilon);
c4 = f(phi+c3*dt, eta, field, epsilon);
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
end

function y = Heun_step(phi, dt, eta, theta, epsilon)
nx = length(phi);
field = randn(nx,1)*sqrt(2*eta*theta/dt);
z = myFTCS(phi, dt, eta, field, epsilon);
y = myFTCS((phi+z)/2, dt, eta, field, epsilon);
end

function y = f(x, eta, field, epsilon)
y = zeros(length(x),2);
phi = x(:,1);
dphi = x(:,2);
phixx = ((circshift(phi,1) + circshift(phi,-1)) - 2*phi);
% V = (phi.^3 - epsilon*phi)/2;
V = (phi.^2 - epsilon)/2;
V = V.*phi;
y(:,1) = dphi;
y(:,2) = field - eta*dphi + phixx - V;
end