% script IR_fun_v2_test.m
% test script for IR_fun_v2.m

% vary TI and T1
nTI = 500;
nT1 = 7;

M0  = ones(nT1,1) * exp(1i * 0);
T1  = logspace(1,4,nT1);
T2  = ones(nT1,1) * 80;

TR  = 1400 * ones(nTI,1);
TI  = linspace(1,1400-1,nTI);
TE  = ones(nTI,1) * 14;

s_ir = NaN(nT1,nTI);
for i = 1:nTI
  s_ir(:,i) = IR_fun_v2(M0, T1, T2, TR(i), TI(i), TE(i));
end

figure, plot(TI, real(s_ir));
figure, plot(TI, imag(s_ir));
