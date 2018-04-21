% script IR_fun_v3_test.m
% test script for IR_fun_v3.m

%% vary TI and T1
nTI = 500;
nT1 = 7;

M0  = ones(nT1,1) * exp(1i * 0);
T1  = logspace(2.7,3.1,nT1);
T2  = ones(nT1,1) * 80;

TR  = 1400 * ones(nTI,1);
TI  = linspace(1,1400-1,nTI);
TE  = ones(nTI,1) * 14;

sp  = NaN(nT1,nTI);
sm  = NaN(nT1,nTI);
for i = 1:nTI
  [sp(:,i), sm(:,i)] = ...
    IR_fun_v3(M0, T1, T2, TR(i), TI(i), TE(i));
end

for i = 1:nT1
  str{i} = sprintf('T1=%4.3fms', T1(i));
end
% figure, plot(TI, real(sp));
% figure, plot(TI, real(sm));
figure, plot(TI, imag(sp)); legend(str{:}, 'Location', 'SE'); xlabel('TI (ms)');
figure, plot(TI, imag(sm)); legend(str{:}, 'Location', 'NE'); xlabel('TI (ms)');


%% vary TI and kap
nTI = 500;
nk  = 7;

M0  = ones(nk,1) * exp(1i * 0);
T1  = ones(nk,1) * median(T1);
T2  = ones(nk,1) * 80;

TR  = 1400 * ones(nTI,1);
TI  = linspace(1,1400-1,nTI);
TE  = ones(nTI,1) * 14;

kap = linspace(0.5,1.5,nk)';

sp  = NaN(nk,nTI);
sm  = NaN(nk,nTI);
for i = 1:nTI
  [sp(:,i), sm(:,i)] = ...
    IR_fun_v3(M0, T1, T2, TR(i), TI(i), TE(i), 'kap', kap);
end

for i = 1:nk
  str{i} = sprintf('kap=%2.1f', kap(i));
end
figure; plot(TI, imag(sp)); legend(str{:}, 'Location', 'SE'); xlabel('TI (ms)');
figure; plot(TI, imag(sm)); legend(str{:}, 'Location', 'NE'); xlabel('TI (ms)');
% figure; plot(TI, abs(sp)-abs(sm));