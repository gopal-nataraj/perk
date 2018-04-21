% script IR_fun_v4_test.m
% test script for IR_fun_v4.m

%% vary TI and inveff
nTI = 500;
nIv = 7;

M0  = ones(nIv,1) * exp(1i * 0);
T1  = ones(nIv,1) * 1000;
T2  = ones(nIv,1) * 80;

TR  = 1400 * ones(nTI,1);
TI  = linspace(1,1400-1,nTI);
TE  = ones(nTI,1) * 14;

inveff = linspace(0.7,1.3,nIv)';

sp  = NaN(nIv,nTI);
sn  = NaN(nIv,nTI);
for i = 1:nTI
  [sp(:,i), sn(:,i)] = ...
    IR_fun_v4(M0, T1, T2, TR(i), TI(i), TE(i), 'inveff', inveff);
end

for i = 1:nIv
  str{i} = sprintf('inveff=%2.1f', inveff(i));
end
figure; plot(TI, imag(sp)); legend(str{:}, 'Location', 'SE'); xlabel('TI (ms)');
figure; plot(TI, imag(sn)); legend(str{:}, 'Location', 'NE'); xlabel('TI (ms)');


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
sn  = NaN(nT1,nTI);
for i = 1:nTI
  [sp(:,i), sn(:,i)] = ...
    IR_fun_v4(M0, T1, T2, TR(i), TI(i), TE(i));
end

for i = 1:nT1
  str{i} = sprintf('T1=%4.3fms', T1(i));
end
% figure, plot(TI, real(sp));
% figure, plot(TI, real(sm));
figure, plot(TI, imag(sp)); legend(str{:}, 'Location', 'SE'); xlabel('TI (ms)');
figure, plot(TI, imag(sn)); legend(str{:}, 'Location', 'NE'); xlabel('TI (ms)');


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
sn  = NaN(nk,nTI);
for i = 1:nTI
  [sp(:,i), sn(:,i)] = ...
    IR_fun_v4(M0, T1, T2, TR(i), TI(i), TE(i), 'kap', kap);
end

for i = 1:nk
  str{i} = sprintf('kap=%2.1f', kap(i));
end
figure; plot(TI, imag(sp)); legend(str{:}, 'Location', 'SE'); xlabel('TI (ms)');
figure; plot(TI, imag(sn)); legend(str{:}, 'Location', 'NE'); xlabel('TI (ms)');