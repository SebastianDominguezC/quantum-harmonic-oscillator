clear all; clc;

global M h w s;

M = 2;
h = 1;
w = 1;
s = sqrt(2 * h / M / w);

% 2.1.1
delta = 0.01;
min = -2.2;
max = 2.2;
x = min : delta : max;
y = min : delta : max;

[X,Y] = meshgrid(x, y);
phi_1_2 = nx_ny_eigenstate(1, 2, X, Y);
phi_3_2 = nx_ny_eigenstate(3, 2, X, Y);

figure(1);
imagesc(x, y, conj(phi_1_2) .* phi_1_2);
xlabel('x');
ylabel('y');
title('Densidad probabilidad nx = 1, ny = 2');

figure(2);
imagesc(x, y, conj(phi_3_2) .* phi_3_2);
xlabel('x');
ylabel('y');
title('Densidad probabilidad nx = 3, ny = 2');

% 2.1.2
C = (469 / 900)^(-1/2);
phi_0_1 = nx_ny_eigenstate(0, 1, X, Y);
phi_3_1 = nx_ny_eigenstate(3, 1, X, Y);

superposition = (1/3) * phi_0_1 + (-2i/5) * phi_3_1 + (1/2) * phi_1_2;
phi = C * superposition;

density = conj(phi) .* phi;

figure(3);
imagesc(x, y, density);
xlabel('x');
ylabel('y');
title('Densidad probabilidad de la superposicion');

phase = angle(phi);
figure(4);
imagesc(x, y, phase);
xlabel('x');
ylabel('y');
title('Fase de la superposicion');

% 2.1.3
[phi_grad_x, phi_grad_y] = gradient(phi);
Jx = h / M * imag(conj(phi) .* phi_grad_x);
Jy = h / M * imag(conj(phi) .* phi_grad_y);

figure(5);
quiver(X(1:21:end, 1:21:end) , Y(1:21:end, 1:21:end), Jx(1:21:end, 1:21:end), Jy(1:21:end, 1:21:end), 2);
xlabel('x');
ylabel('y');
title('Corriente de probabilidad');

% 2.1.4
dt = 0.1;
t = 0 : dt : 2 * 2 * pi / w;

phi_0_1_t = propagator(t, x, y, 0, 1, phi_0_1);
phi_3_1_t = propagator(t, x, y, 3, 1, phi_3_1);
phi_1_2_t = propagator(t, x, y, 1, 2, phi_1_2);

superposition_t = (1/3) * phi_0_1_t + (-2i/5) * phi_3_1_t + (1/2) * phi_1_2_t;

phi_t = C * superposition_t;

phi_density_t = conj(phi_t) .* phi_t;

figure(6);
animation = imagesc(x, y, phi_density_t(:, :, 1));
xlabel('x');
ylabel('y');
title('Densidad de probabilidad');
axis equal;

for i = 1:length(t)
    animation.CData = phi_density_t(:, :, i);
    pause(0.1);
end

function phi_nx_ny = nx_ny_eigenstate(nx, ny, x, y)
    global s;
    r = sqrt(x .^ 2 + y .^ 2);

    c = 1 / (s * sqrt( 2 ^ (nx + ny - 1) * pi * factorial(nx) * factorial(ny)));
    exponent = exp(-r .^ 2 / s .^ 2);
    hermite_nx = polyval(HermitePoly(nx), sqrt(2) .* x / s);
    hermite_ny = polyval(HermitePoly(ny), sqrt(2) .* y / s);

    phi_nx_ny = c .* exponent .* hermite_nx .* hermite_ny;
end

function wiggle = wiggle_factor(nx, ny, t)
    global h w;
    q = nx + ny;
    E = (q + 1) * h * w;
    wiggle = exp(-1i * E * t / h);
end

function time_states = propagator(t, x, y, nx, ny, phi_nx_ny)
    time_states = zeros(length(x), length(y), length(t));
    for i = 1:length(t)
        time_states(:, :, i) = phi_nx_ny .* wiggle_factor(nx, ny, t(i));
    end
end

