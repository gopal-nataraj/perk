% script test.m
% validation of dess models, gradients, hessians
% 
% copyright 2018, gopal nataraj, university of michigan
%
% version control
%   1.1     2018-01-06      original; validated magnitude models only

%% vary m0
m0 = col(linspace(eps, 1, 1000));
t1 = ones(size(m0))*1000;
t2 = ones(size(m0))*50;
kap = ones(size(m0))*1;
flip = 20 * (pi/180);
tr = 15;
tep = 4;
tem = 4;

bool.mag = true;
[sp, sm] = dess_fun_v2(...
  m0, t1, t2, flip, tr, tep, tem,...
  'kap', kap,...
  'mag', bool.mag);
[spgradx, smgradx] = dess_gradx(...
  m0, t1, t2, flip, tr, tep, tem,...
  'kap', kap,...
  'mag', bool.mag);
[sphessx, smhessx] = dess_hessx(...
  m0, t1, t2, flip, tr, tep, tem,...
  'kap', kap,...
  'mag', bool.mag);

figure; 
subplot(3,1,1);
plot(m0, sp, m0, sm);
ylabel('mag sig (a.u.)');
legend('dofocusing','refocusing');
subplot(3,1,2);
plot(m0, spgradx(:,1), m0, smgradx(:,1));
ylabel('1st deriv of mag sig w.r.t. abs(m0) (a.u.)');
subplot(3,1,3);
plot(m0, sphessx(:,1,1), m0, smhessx(:,1,1));
ylabel('2nd deriv of mag sig w.r.t. abs(m0) (a.u.)');
xlabel('m0 (a.u.)');

%% vary t1
t1 = col(logspace(log10(10^-5),log10(100), 1000));
m0 = ones(size(t1));
t2 = ones(size(t1))*50;
kap = ones(size(t1))*1;
flip = 20 * (pi/180);
tr = 15;
tep = 4;
tem = 4;

bool.mag = true;
[sp, sm] = dess_fun_v2(...
  m0, t1, t2, flip, tr, tep, tem,...
  'kap', kap,...
  'mag', bool.mag);
[spgradx, smgradx] = dess_gradx(...
  m0, t1, t2, flip, tr, tep, tem,...
  'kap', kap,...
  'mag', bool.mag);
[sphessx, smhessx] = dess_hessx(...
  m0, t1, t2, flip, tr, tep, tem,...
  'kap', kap,...
  'mag', bool.mag);

figure;
subplot(3,1,1);
plot(t1, sp, t1, sm);
ylabel('mag sig (a.u.)');
legend('defocusing','refocusing');
subplot(3,1,2);
plot(t1, spgradx(:,2), t1, smgradx(:,2));
ylabel('1st deriv of mag sig w.r.t. t1 (1/ms)');
subplot(3,1,3);
plot(t1, sphessx(:,2,2), t1, smhessx(:,2,2));
ylabel('2nd deriv of mag sig w.r.t. t1 (1/ms^2)');
xlabel('t1 (ms)');

%% vary t2
t2 = col(logspace(log10(10^-4), log10(30), 100));
m0 = ones(size(t2));
t1 = ones(size(t2))*1000;
kap = ones(size(t2))*1;
flip = 20 * (pi/180);
tr = 15;
tep = 4;
tem = 4;

bool.mag = true;
[sp, sm] = dess_fun_v2(...
  m0, t1, t2, flip, tr, tep, tem,...
  'kap', kap,...
  'mag', bool.mag);
[spgradx, smgradx] = dess_gradx(...
  m0, t1, t2, flip, tr, tep, tem,...
  'kap', kap,...
  'mag', bool.mag);
[sphessx, smhessx] = dess_hessx(...
  m0, t1, t2, flip, tr, tep, tem,...
  'kap', kap,...
  'mag', bool.mag);

figure;
subplot(3,1,1);
plot(t2, sp, t2, sm);
ylabel('mag sig (a.u.)');
legend('defocusing','refocusing');
subplot(3,1,2);
plot(t2, spgradx(:,3), t2, smgradx(:,3));
ylabel('1st deriv of mag sig w.r.t. t2 (1/ms)');
subplot(3,1,3);
plot(t2, sphessx(:,3,3), t2, smhessx(:,3,3));
ylabel('2nd deriv of mag sig w.r.t. t2 (1/ms^2)');
xlabel('t2 (ms)');
