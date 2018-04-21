% script test.m
% validation of spgr model, gradient, hessian
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
te = 4;

bool.mag = true;
s = spgr_fun_v2(...
  m0, t1, t2, flip, tr, te,...
  'kap', kap,...
  'mag', bool.mag);
sgradx = spgr_gradx(...
  m0, t1, t2, flip, tr, te,...
  'kap', kap,...
  'mag', bool.mag);
shessx = spgr_hessx(...
  m0, t1, t2, flip, tr, te,...
  'kap', kap,...
  'mag', bool.mag);

figure; 
subplot(3,1,1);
plot(m0, s);
ylabel('mag signal (a.u.)');
subplot(3,1,2);
plot(m0, sgradx(:,1));
ylabel('1st deriv of mag sig w.r.t. abs(m0) (a.u.)');
subplot(3,1,3);
plot(m0, shessx(:,1,1));
ylabel('2nd deriv of mag sig w.r.t. abs(m0) (a.u.)');
xlabel('m0 (a.u.)');

%% vary t1
t1 = col(logspace(log10(10^-5),log10(1000), 1000));
m0 = ones(size(t1));
t2 = ones(size(t1))*50;
kap = ones(size(t1))*1;
flip = 20 * (pi/180);
tr = 15;
te = 4;

bool.mag = true;
s = spgr_fun_v2(...
  m0, t1, t2, flip, tr, te,...
  'kap', kap,...
  'mag', bool.mag);
sgradx = spgr_gradx(...
  m0, t1, t2, flip, tr, te,...
  'kap', kap,...
  'mag', bool.mag);
shessx = spgr_hessx(...
  m0, t1, t2, flip, tr, te,...
  'kap', kap,...
  'mag', bool.mag);

figure; 
subplot(3,1,1);
plot(t1, s);
ylabel('mag signal (a.u.)');
subplot(3,1,2);
plot(t1, sgradx(:,2));
ylabel('1st deriv of mag sig w.r.t. t1 (1/ms)');
subplot(3,1,3);
plot(t1, shessx(:,2,2));
ylabel('2nd deriv of mag sig w.r.t. abs(m0) (1/ms^2)');
xlabel('t1 (ms)');

%% vary t2
t2 = col(logspace(log10(10^-4), log10(300), 100));
m0 = ones(size(t2));
t1 = ones(size(t2))*1000;
kap = ones(size(t2))*1;
flip = 20 * (pi/180);
tr = 15;
te = 4;

bool.mag = true;
s = spgr_fun_v2(...
  m0, t1, t2, flip, tr, te,...
  'kap', kap,...
  'mag', bool.mag);
sgradx = spgr_gradx(...
  m0, t1, t2, flip, tr, te,...
  'kap', kap,...
  'mag', bool.mag);
shessx = spgr_hessx(...
  m0, t1, t2, flip, tr, te,...
  'kap', kap,...
  'mag', bool.mag);

figure; 
subplot(3,1,1);
plot(t2, s);
ylabel('mag signal (a.u.)');
subplot(3,1,2);
plot(t2, sgradx(:,3));
ylabel('1st deriv of mag sig w.r.t. t2 (1/ms)');
subplot(3,1,3);
plot(t2, shessx(:,3,3));
ylabel('2nd deriv of mag sig w.r.t. t2 (1/ms^2)');
xlabel('t2 (ms)');