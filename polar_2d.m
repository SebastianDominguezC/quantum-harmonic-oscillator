clear all; clc; 

global M h a w s;

M = 2;
h = 1;
a = 1;
w = 1;
s = sqrt(2 * h / M / w);

max = 3;
min = -3;
delta = 0.01;
x = min:delta:max;
y = min:delta:max;

% 2.2.1
[X, Y] = meshgrid(x, y);

[R, T] = get_polar(X, Y);

psi_2_1 = hidrogen_nm(2, 1, R, T);
density_2_1 = conj(psi_2_1) .* psi_2_1;

psi_3_2 = hidrogen_nm(3, 2, R, T);
density_3_2 = conj(psi_3_2) .* psi_3_2;

axis equal;

figure(1);
imagesc(x, y, density_2_1);
xlabel('x');
ylabel('y');
title('Estado polar (2, 1)');

figure(2);
imagesc(x, y, density_3_2);
xlabel('x');
ylabel('y');
title('Estado polar (3, 2)');

figure(3);
imagesc(x, y, angle(psi_2_1));
xlabel('x');
ylabel('y');
title('Fase polar (2, 1)');

figure(4);
imagesc(x, y, angle(psi_3_2));
xlabel('x');
ylabel('y');
title('Fase polar (3, 2)');

% 2.2.2
dt = 0.1;
t = 0 : dt : 2 * 2 * pi / w;

C = (469 / 900)^(-1/2);

psi_0_2 = hidrogen_nm(0, 2, R, T);
psi_1_m1 = hidrogen_nm(1, -1, R, T);
psi_2_0 = hidrogen_nm(2, 0, R, T);

psi_0_2_t = propagator(t, x, y, 0, 2, psi_0_2);
psi_1_m1_t = propagator(t, x, y, 1, -1, psi_1_m1);
psi_2_0_t = propagator(t, x, y, 2, 0, psi_2_0);

superposition = (1/3) * psi_0_2_t + (-2i/5) * psi_1_m1_t + (1/2) * psi_2_0_t;
psi = C * superposition;

phi_density_t = conj(psi) .* psi;


figure(5);
imagesc(x, y, phi_density_t(:, :, 1));
axis equal;
xlabel('x');
ylabel('y');
title('Densidad de superposicion polar');

figure(6);
animation = imagesc(x, y, phi_density_t(:, :, 1));
axis equal;
xlabel('x');
ylabel('y');
title('Densidad de superposicion polar');

for i = 1:length(t)
    animation.CData = phi_density_t(:, :, i);
    pause(0.1);
end

figure(6);
axis equal;
imagesc(x, y, angle(psi(:, :, 1)));

function [r, t] = get_polar(x, y)
    r = sqrt(x .^ 2 + y .^ 2);
    t = atan2(y, x);
end

function state = hidrogen_nm(n, m, r, theta)
    global s;
    p = 2 * r .^ 2 / s ^ 2;

    c1 = 2 * factorial(n);
    c2 = pi * factorial(n + abs(m)) * s ^ 2;
    c = sqrt(c1 / c2);

    R = radial(r, m);
    Y = spherical_harmonic(n, m, p, theta);

    state = c .* R .* Y;
end

function Y = spherical_harmonic(n, m, p, theta)
    harmonic = associated_laguerre_poly(n, abs(m), p);
    imaginary = exp(1i * m .* theta);
    Y = harmonic .* imaginary;
end 

function R = radial(r, m)
    global s;
    power = (sqrt(2) * r / s) .^ abs(m);
    exponent = exp(-r .^ 2 / s ^ 2);

    R = power .* exponent;
end

function wiggle = wiggle_factor(n, m, t)
    global h w;
    q = 2 * n + abs(m);
    E = (q + 1) * h * w;
    wiggle = exp(-1i * E * t / h);
end

function time_states = propagator(t, r, theta, n, m, phi_nx_ny)
    time_states = zeros(length(r), length(theta), length(t));
    for i = 1:length(t)
        time_states(:, :, i) = phi_nx_ny .* wiggle_factor(n, m, t(i));
    end
end