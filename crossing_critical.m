clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

T_min = -100;
T_max = 350;
dt = 1e-2;
t = T_min:dt:T_max;
nt = length(t);
dx = 1;
L = 1000;
x = dx:dx:L;
nx = length(x);
eta = 1;
theta = 1e-2;
% epsilon = 1;
tau_Q = 64;
lambda = dt/dx^2;

phi = zeros(nx,nt);
dphi = zeros(nx,nt);
% phi(floor(nx/2),1) = 1;
phi(:,1) = rand(nx,1);

for i = 2:nt 
   [y, dy] = Heun_step(phi(:,i-1), dphi(:,i-1), dx, dt, eta, theta, min(1,t(i)/tau_Q));
   phi(:,i) = y;
   dphi(:,i) = dy;
end

toc;

titlename = strcat('L = ',num2str(L), ', \tau_Q = ', num2str(tau_Q), ', seed = ', num2str(myseed), ', dt = ', num2str(dt));
figure;
set(gcf, 'position', [250 70 1000 900]);

plot_result = tiledlayout(5,1);

ax1 = nexttile;
plot(ax1,x,phi(:,round((-80-T_min)/dt)))

ax1 = nexttile;
plot(ax1,x,phi(:,round((7.5-T_min)/dt)))

ax1 = nexttile;
plot(ax1,x,phi(:,round((32.5-T_min)/dt)))

ax1 = nexttile;
plot(ax1,x,phi(:,round((45-T_min)/dt)))

ax1 = nexttile;
plot(ax1,x,phi(:,round((333-T_min)/dt)))

title(plot_result,titlename)
xlabel(plot_result,'x')
ylabel(plot_result,'\phi')
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