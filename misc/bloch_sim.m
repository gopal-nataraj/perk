% Bloch Equation Simulator

% Parameters
T1 = 1000;
T2 = 100;
TR = 5000; 
TE = 5;
flip = 20;

% Magnetization
T = 10000; 
dt = 0.01; 
t = [0:dt:T]';
nt = length(t);

% Rotation Matrix
R = [1 0 0; 0 cos(flip) sin(flip); 0 -sin(flip) cos(flip)];

% Relaxation Matrix
E = diag([exp(-dt/T2) exp(-dt/T2) exp(-dt/T1)]);

M = NaN(3, nt);
M(:,1) = [1 0 0]';

for i = 2:nt
    if rem(t(i), TR) == 0
        M(:,i) = R*M(:,i-1);
    else
        M(:,i) = E*M(:,i-1) + (1-exp(-dt/T1))*[0 0 1]';
    end
end

% Pool transverse magnetization
Mxy = M(1,:) + 1j*M(2,:);

figure(1); 
subplot(1,2,1); plot(t, abs(Mxy)); grid on;
xlabel('Time, t (ms)'); ylabel('M_x_y');
title('Transverse Magnetization, M_x_y');

subplot(1,2,2); plot(t, M(3,:)); grid on;
xlabel('Time, t (ms)'); ylabel('M_z');
title('Longitudinal Magnetization, M_z');

    
