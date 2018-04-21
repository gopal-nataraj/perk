% script SE_fun_v4_test.m
% test script for SE_fun_v4.m

% vary TE
nTE = 500;
nT2 = 5;

M0  = ones(nT2,1) * exp(1i * pi/2);
T1  = ones(nT2,1) * 1200;
T2  = logspace(1,3,nT2)';

TR  = 2000 * ones(nTE,1);
TE  = linspace(1,2000-1,nTE);

s_se = NaN(nT2,nTE);
for i = 1:nTE
  s_se(:,i) = SE_fun_v4(M0, T1, T2, TR(i), TE(i));
end

figure, plot(TE, real(s_se));
figure, plot(TE, imag(s_se));