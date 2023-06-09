clear;
% close all;
clc;
format long
tic;

No = 26;

rng(No * 1e3)

T_max = 3000;
dt = 1e-3;
M = 100;
t = 0:dt * M:T_max;
nt = length(t);
L_all = 25:5:50;
nL = length(L_all);
eta = 0.05;
theta = 10;
epsilon0 = 5;
omega = 1;
N_comp = 2;

order = zeros(N_comp, nL, nt);
order_total = zeros(nL, nt);

for n = 1:nL
    L = L_all(n);
    nx = L ^ 2;

    dphi = zeros(L, L, N_comp);
    phi = epsilon0 * ones(L, L, N_comp);
    % phi = epsilon0*randn(L);
    order(:, n, 1) = sum(phi, 1:2) / nx;
    order_total(n, 1) = sqrt(sum(order(:, n, 1) .^ 2));

    t_it = 0;
    count = 2;

    for i = 2:round(T_max / dt) + 1
        t_it = t_it + dt;
        %     epsilon = epsilon0 + delta*cos(omega*2*pi*t(i));
        [phi, dphi] = Heun_step(phi, dphi, t_it, dt, eta, theta, epsilon0, omega * 2 * pi, N_comp);
        %     [phi, dphi] = myFTCS(phi, dphi, t_it, dt, eta, 0, epsilon0, omega*2*pi);
        %     order(i) = sum(sqrt(phi(1,:,:,:).^2 + phi(2,:,:,:).^2),'all');
        if mod(i + 1, M) == 0
            order(:, n, count) = sum(phi, 1:2) / nx;
            order_total(n, count) = sqrt(sum(order(:, n, count) .^ 2));
            count = count + 1;
        end

    end

end

order_amp = sqrt(mean(order(:, :, floor(nt / 2) + 1:end) .^ 2, 3));
[~, order_rank] = max(order_amp, [], 1);
order_domi = zeros(nL, nt);
for j = 1:nL
    order_domi(j, :) = order(order_rank(1, j), j, :);
end

cut = 1000;
phi_f = abs(fft(order_domi(:,floor(nt/2)+1:end),(nt+1)/2,2));
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

fprintf('eta = %g, epsilon0 = %g, peak_f = %g, peak_per = %g .\n', eta, epsilon0, peak_f, peak_per);

clear i n

%save(strcat('main_dt', replace(num2str(dt), '.', '_'), '_theta', replace(num2str(theta), '.', '_'), '_No', num2str(No), '.mat'));

toc;

le = cell(1, nL);

for ii = 1:nL
    le{ii} = strcat('L = ', num2str(L_all(ii)));
end

%{
figure;
plot(t, order_total)
xlabel('t')
ylabel('order total')
saveas(gcf, strcat('order_total_dt', replace(num2str(dt),'.','_'), '_theta', replace(num2str(theta),'.','_'), '_No', num2str(No), '.fig'));
print(strcat('order_total_dt', replace(num2str(dt),'.','_'), '_theta', replace(num2str(theta),'.','_'), '_No', num2str(No)), '-dpng', '-r0');


figure;
plot(t(floor(80 * nt / 100):end), order_total(:, floor(80 * nt / 100):end))
xlabel('t')
ylabel('order')
legend(le);
saveas(gcf, strcat('order_part_dt', replace(num2str(dt), '.', '_'), '_theta', replace(num2str(theta), '.', '_'), '_No', num2str(No), '.fig'));
print(strcat('order_part_dt', replace(num2str(dt), '.', '_'), '_theta', replace(num2str(theta), '.', '_'), '_No', num2str(No)), '-dpng', '-r0');

figure;
plot(t, order_total)
xlabel('t')
ylabel('order')
legend(le);
saveas(gcf, strcat('order_dt', replace(num2str(dt), '.', '_'), '_theta', replace(num2str(theta), '.', '_'), '_No', num2str(No), '.fig'));
print(strcat('order_dt', replace(num2str(dt), '.', '_'), '_theta', replace(num2str(theta), '.', '_'), '_No', num2str(No)), '-dpng', '-r0');


figure;
plot(t, order_domi)
xlabel('t')
ylabel('order_{dominant}')
legend(le);
saveas(gcf, strcat('order_domi_dt', replace(num2str(dt), '.', '_'), '_theta', replace(num2str(theta), '.', '_'), '_No', num2str(No), '.fig'));
print(strcat('order_domi_dt', replace(num2str(dt), '.', '_'), '_theta', replace(num2str(theta), '.', '_'), '_No', num2str(No)), '-dpng', '-r0');

figure;
plot(w_main, phi_f_main)
xlabel('\omega')
ylabel('order_f')
legend(le);
saveas(gcf, strcat('order_f_dt', replace(num2str(dt), '.', '_'), '_theta', replace(num2str(theta), '.', '_'), '_No', num2str(No), '.fig'));
print(strcat('order_f_dt', replace(num2str(dt), '.', '_'), '_theta', replace(num2str(theta), '.', '_'), '_No', num2str(No)), '-dpng', '-r0');
%}

function [y, z] = myFTCS(phi, dphi, t, dt, eta, field, epsilon, omega)
    [c1, d1] = f(phi, dphi, eta, field, epsilon * (cos(omega * t) + 1));
    [c2, d2] = f(phi + c1 * dt / 2, dphi + d1 * dt / 2, eta, field, epsilon * (cos(omega * (t + dt / 2)) + 1));
    [c3, d3] = f(phi + c2 * dt / 2, dphi + d2 * dt / 2, eta, field, epsilon * (cos(omega * (t + dt / 2)) + 1));
    [c4, d4] = f(phi + c3 * dt, dphi + d3 * dt, eta, field, epsilon * (cos(omega * (t + dt)) + 1));
    y = phi + dt * (c1 + 2 * c2 + 2 * c3 + c4) / 6;
    z = dphi + dt * (d1 + 2 * d2 + 2 * d3 + d4) / 6;
end

function [y1, y2] = Heun_step(phi, dphi, t, dt, eta, theta, epsilon, omega, N_comp)
    L = length(phi);
    field = randn(L, L, N_comp) * sqrt(2 * eta * theta / dt);
    [z1, z2] = myFTCS(phi, dphi, t, dt, eta, field, epsilon, omega);
    [y1, y2] = myFTCS((phi + z1) / 2, (dphi + z2) / 2, t, dt, eta, field, epsilon, omega);
end

function [y, z] = f(x, dx, eta, field, epsilon)

    phi_l = circshift(x, 1, 1);
    phi_r = circshift(x, -1, 1);
    phi_u = circshift(x, 1, 2);
    phi_d = circshift(x, -1, 2);

    phixx = phi_l + phi_r + phi_u + phi_d - 4 * x;

    phi_norm = sum(x .^ 2, 3);
    V = (phi_norm - epsilon) / 2;
    V = V .* x;
    y = dx;
    z = field - eta * dx + phixx - V;
end
