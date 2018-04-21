%{
% DESS SSFP Integrations
syms B;
syms E1 E2 a positive;

p = 1 - E1*cos(a) - E2^2*(E1-cos(a));
q = E2 * (1-E1) * (1+cos(a));
f = (1 - E2*exp(1i*B)) / (p - q*cos(B));
int(f, B, 0, 2*pi)
int(cos(B)/(p-q*cos(B)), B, 0, 2*pi);
int(sin(B)/(p-q*cos(B)), B, 0, 2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Symbolic Calculations
syms E1 E2 a positive;
assume(abs(E1) < 1);
assume(abs(E2) < 1);
assume(abs(a) <= pi/2); 

% a = pi/6;
% TR = 10;
% T1 = 1000;
% T2 = 100;
% E1 = exp(-TR./T1);
% E2 = exp(-TR./T2);

p = 1 - E1*cos(a) - E2^2*(E1-cos(a));
q = E2 * (1-E1) * (1+cos(a));
r = (1-E2^2)./sqrt(p.^2-q.^2);
Splus_1 = tan(a/2).*(1 - r.*(E1-cos(a)));
Sminus_1 = tan(a/2).*(1 - r.*(1 - E1.*cos(a)));
Smom_1 = Sminus_1 ./ Splus_1;

v1 = (1-E1.*cos(a)) ./ (E1-cos(a));
Splus_2 = tan(a/2).*(1 - sqrt((1-E2^2)./(v1.^2-E2^2)));
%Sminus_2 = tan(a/2).*(1 - v1.*sqrt((1-E2^2)./(v1.^2-E2^2)));
Sminus_2 = tan(a/2).*(1 - sqrt((1-E2^2)./(1-E2^2./v1.^2)));
Smom_2 = Sminus_2 ./ Splus_2;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculations for varying alpha
M0 = 1;
a = linspace(-pi,pi,1000);
TR = 10;
T1 = 500; 
T2 = 100;
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);

p = 1 - E1*cos(a) - E2^2*(E1-cos(a));
q = E2 * (1-E1) * (1+cos(a));
r = (1-E2^2)./sqrt(p.^2-q.^2);
Splus_1 = 1i*M0 * tan(a/2).*(1 - r.*(E1-cos(a)));
Sminus_1 = 1i*M0 * tan(a/2).*(1 - r.*(1 - E1.*cos(a)));
Smom_1 = Sminus_1 ./ Splus_1;

v1 = (1-E1.*cos(a)) ./ (E1-cos(a));
R = sqrt((1-E2^2)./(1-E2^2./v1.^2));
Splus_2 = 1i*M0 * tan(a/2).*(1 - (R./v1));
Sminus_2 = 1i*M0 * tan(a/2).*(1 - R);
Smom_2 = Sminus_2 ./ Splus_2;

% figure; plot(180*a/pi, R); 
% figure; plot(180*a/pi, 1./v1);
% figure; plot(180*a/pi, R./v1);

% Plots for varying alpha
figure;
subplot 321, plot(a, imag(Splus_1)); axis([-pi pi -0.2 0.2]);
title('|S^+| vs. alpha with Heule Expression'); 
subplot 322, plot(a, imag(Splus_2)); axis([-pi pi -0.2 0.2]);
title('|S^+| vs. alpha with Our Expression'); 
subplot 323, plot(a, imag(Sminus_1)); axis([-pi pi -0.2 0.2]);
title('|S^-| vs. alpha with Heule Expression'); 
subplot 324, plot(a, imag(Sminus_2)); axis([-pi pi -0.2 0.2]);
title('|S^-| vs. alpha with Our Expression'); 
subplot 325, plot(a, abs(Smom_1)); axis([-pi pi 0 1]);
title('|S^-/S^+| vs. alpha with Heule Expression');
subplot 326, plot(a, abs(Smom_2)); axis([-pi pi 0 1]);
title('|S^-/S^+| vs. alpha with Our Expression');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculations for varying E1

M0 = 1;
a = pi/6;
TR = 10;
T1 = logspace(-1,10,10000); %T1 = 1000;
T2 = 100;
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);

p = 1 - E1*cos(a) - E2^2*(E1-cos(a));
q = E2 * (1-E1) * (1+cos(a));
r = (1-E2^2)./sqrt(p.^2-q.^2);
Splus_1 = 1i*M0 * tan(a/2).*(1 - r.*(E1-cos(a)));
Sminus_1 = 1i*M0 * tan(a/2).*(1 - r.*(1 - E1.*cos(a)));
Smom_1 = Sminus_1 ./ Splus_1;

v1 = (1-E1.*cos(a)) ./ (E1-cos(a));
Splus_2 = 1i*M0 * tan(a/2).*(1 - (1./v1).*sqrt((1-E2^2)./(1-E2^2./v1.^2)));
Sminus_2 = 1i*M0 * tan(a/2).*(1 - sqrt((1-E2^2)./(1-E2^2./v1.^2)));
Smom_2 = Sminus_2 ./ Splus_2;

% Plots for varying E1
figure;
subplot 321, plot(E1, imag(Splus_1)); axis([0 1 0 0.4]);
title('|S^+| vs. E1 with Heule Expression'); 
subplot 322, plot(E1, imag(Splus_2));axis([0 1 0 0.4]);
title('|S^+| vs. E1 with Our Expression'); 
subplot 323, plot(E1, imag(Sminus_1)); axis([0 1 0 0.2]);
title('|S^-| vs. E1 with Heule Expression'); 
subplot 324, plot(E1, imag(Sminus_2)); axis([0 1 0 0.2]);
title('|S^-| vs. E1 with Our Expression'); 
subplot 325, plot(E1, imag(Smom_1)); axis([0 1 0 1]);
title('|S^-/S^+| vs. E1 with Heule Expression');
subplot 326, plot(E1, imag(Smom_2)); axis([0 1 0 1]);
title('|S^-/S^+| vs. E1 with Our Expression');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Pretty Figures for Quals II Report/Talk
a = 0:0.001:pi/2; %a = pi/6;
TR = 10;
T1 = 1000;  
T2 = 100;
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);
v1 = (1-E1.*cos(a)) ./ (E1-cos(a));
Splus_2 = tan(a/2).*(1 - (1./v1).*sqrt((1-E2.^2)./(1-E2.^2./v1.^2)));
Sminus_2 = tan(a/2).*(1 - sqrt((1-E2.^2)./(1-E2.^2./v1.^2)));
Smom_2 = Sminus_2 ./ Splus_2;

% subplot 431, plot(180*a/pi, Splus_2); axis([0 90 0 0.15]); grid on;
% title('|S^+|, varying alpha'); xlabel('alpha (degrees)'); ylabel('|S^+|');
% subplot 432, plot(180*a/pi, Sminus_2); axis([0 90 0 0.15]); grid on;
% title('|S^-|, varying alpha'); xlabel('alpha (degrees)'); ylabel('|S^-|');
% subplot 433, plot(180*a/pi, Smom_2, 'b', 180*a/pi, E2^2, 'r');
% axis([0 90 0 1]); legend('|S^-/S^+|', 'MOM', 'Location', 'SE'); grid on;
% title('|S^-/S^+|, varying alpha'); xlabel('alpha (degrees)'); ylabel('|S^-/S^+|');

figure; plot(180*a/pi, Splus_2); axis([0 90 0 0.15]); grid on;
% print -depsc splus_vs_alpha.eps
figure; plot(180*a/pi, Sminus_2); axis([0 90 0 0.15]); grid on;
% print -depsc sminus_vs_alpha.eps
figure; plot(180*a/pi, Smom_2, 'b', 180*a/pi, E2.^2, 'r'); axis([0 90 0 1]); grid on;
% print -depsc ratio_vs_alpha.eps

figure; plot(180*a/pi, Splus_2, 'LineWidth', 2); axis([0 90 0 0.15]); grid on;
% print -depsc splus_vs_alpha.eps
figure; plot(180*a/pi, Sminus_2, 'LineWidth', 2); axis([0 90 0 0.15]); grid on;
% print -depsc sminus_vs_alpha.eps
figure; plot(180*a/pi, Smom_2, 'b', 180*a/pi, E2.^2, 'r', 'LineWidth', 2); axis([0 90 0 1]); grid on;
% print -depsc ratio_vs_alpha.eps

a = pi/6; 
TR = logspace(0, 2, 1000); %TR = 10;
T1 = 1000;
T2 = 100;
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);
v1 = (1-E1.*cos(a)) ./ (E1-cos(a));
Splus_2 = tan(a/2).*(1 - (1./v1).*sqrt((1-E2.^2)./(1-E2.^2./v1.^2)));
Sminus_2 = tan(a/2).*(1 - sqrt((1-E2.^2)./(1-E2.^2./v1.^2)));
Smom_2 = Sminus_2 ./ Splus_2;

% subplot 434, semilogx(TR, Splus_2); axis([1 100 0 0.25]); grid on;
% title('|S^+|, varying TR'); xlabel('TR (ms)'); ylabel('|S^+|');
% subplot 435, semilogx(TR, Sminus_2); axis([1 100 0 0.25]); grid on;
% title('|S^-|, varying TR'); xlabel('TR (ms)'); ylabel('|S^-|');
% subplot 436, semilogx(TR, Smom_2); axis([1 100 0 1]); grid on;
% title('|S^-/S^+|, varying TR'); xlabel('TR (ms)'); ylabel('|S^-/S^+|');

figure; semilogx(TR, Splus_2); axis([1 100 0 0.25]); grid on;
% print -depsc splus_vs_TR.eps
figure; semilogx(TR, Sminus_2); axis([1 100 0 0.25]); grid on;
% print -depsc sminus_vs_TR.eps
figure; semilogx(TR, Smom_2); axis([1 100 0 1]); grid on;
% print -depsc ratio_vs_TR.eps

a = pi/6; 
TR = 10;
T1 = linspace(0.001, 3000, 1000); %T1 = 1000;
T2 = 100;
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);
v1 = (1-E1.*cos(a)) ./ (E1-cos(a));
Splus_2 = tan(a/2).*(1 - (1./v1).*sqrt((1-E2.^2)./(1-E2.^2./v1.^2)));
Sminus_2 = tan(a/2).*(1 - sqrt((1-E2.^2)./(1-E2.^2./v1.^2)));
Smom_2 = Sminus_2 ./ Splus_2;

% subplot 437, plot(T1, Splus_2); axis([0 3000 0 0.25]); grid on;
% title('|S^+|, varying T1'); xlabel('T1 (ms)'); ylabel('|S^+|');
% subplot 438, plot(T1, Sminus_2); axis([0 3000 0 0.25]); grid on;
% title('|S^-|, varying T1'); xlabel('T1 (ms)'); ylabel('|S^-|');
% subplot 439, plot(T1, Smom_2); axis([0 3000 0 1]); grid on;
% title('|S^-/S^+|, varying T1'); xlabel('T1 (ms)'); ylabel('|S^-/S^+|');

figure; plot(T1, Splus_2); axis([0 3000 0 0.25]); grid on;
% print -depsc splus_vs_T1.eps
figure; plot(T1, Sminus_2); axis([0 3000 0 0.25]); grid on;
% print -depsc sminus_vs_T1.eps
figure; plot(T1, Smom_2, 'b', T1, E2.^2, 'r'); axis([0 3000 0 1]); grid on;
% print -depsc ratio_vs_T1.eps

a = pi/6; 
TR = 10;
T1 = 1000;
T2 = linspace(0, 300, 1000); %T2 = 100;
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);
v1 = (1-E1.*cos(a)) ./ (E1-cos(a));
Splus_2 = tan(a/2).*(1 - (1./v1).*sqrt((1-E2.^2)./(1-E2.^2./v1.^2)));
Sminus_2 = tan(a/2).*(1 - sqrt((1-E2.^2)./(1-E2.^2./v1.^2)));
Smom_2 = Sminus_2 ./ Splus_2;

% subplot(4,3,10), plot(T2, Splus_2); axis([0 300 0 0.25]); grid on;
% title('|S^+|, varying T2'); xlabel('T2 (ms)'); ylabel('|S^+|');
% subplot(4,3,11), plot(T2, Sminus_2); axis([0 300 0 0.25]); grid on;
% title('|S^-|, varying T2'); xlabel('T2 (ms)'); ylabel('|S^-|');
% subplot(4,3,12), plot(T2, Smom_2); axis([0 300 0 1]); grid on;
% title('|S^-/S^+|, varying T2'); xlabel('T2 (ms)'); ylabel('|S^-/S^+|');
% print -depsc dess_behavior.eps

figure; plot(T2, Splus_2); axis([0 300 0 0.25]); grid on;
% print -depsc splus_vs_T2.eps
figure; plot(T2, Sminus_2); axis([0 300 0 0.25]); grid on;
% print -depsc sminus_vs_T2.eps
figure; plot(T2, Smom_2); axis([0 300 0 1]); grid on;
title('|S^-/S^+|, varying T2'); xlabel('T2 (ms)'); ylabel('|S^-/S^+|');
% print -depsc ratio_vs_T2.eps


