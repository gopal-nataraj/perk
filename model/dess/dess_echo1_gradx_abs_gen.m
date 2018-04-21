function mag_spxy_tep_gradx = dess_echo1_gradx_abs_gen(M0,T1,T2,kap,R2p,flip,TR,TEp)
%DESS_ECHO1_GRADX_ABS_GEN
%    MAG_SPXY_TEP_GRADX = DESS_ECHO1_GRADX_ABS_GEN(M0,T1,T2,KAP,R2P,FLIP,TR,TEP)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    06-Jan-2018 11:49:11

t2 = 1.0./T1;
t8 = TR.*t2;
t3 = exp(-t8);
t4 = flip.*kap;
t5 = cos(t4);
t6 = 1.0./T2;
t16 = TR.*t6.*2.0;
t7 = exp(-t16);
t9 = -t3+t5;
t10 = t3.*t5;
t11 = t10-1.0;
t12 = flip.*kap.*(1.0./2.0);
t13 = tan(t12);
t14 = abs(t13);
t15 = 1.0./t11;
t17 = -t7+1.0;
t18 = sqrt(t17);
t19 = t9.^2;
t20 = 1.0./t11.^2;
t21 = t7.*t19.*t20;
t22 = t21-1.0;
t23 = 1.0./t22;
t24 = sqrt(-t23);
t25 = -t9.*t15.*t18.*t24+1.0;
t26 = T1.*4.0;
t27 = t5.^2;
t28 = exp(t16);
t29 = T1.*2.0;
t30 = T2+t29;
t31 = R2p+t6;
t32 = exp(-TEp.*t31);
t33 = t7-1.0;
t34 = t23.*t33;
t35 = sqrt(t34);
t36 = -t9.*t15.*t35+1.0;
t37 = 1.0./T2.^2;
mag_spxy_tep_gradx = [t14.*exp(-TEp.*t6.*(R2p.*T2+1.0)).*abs(t25),M0.*1.0./T1.^2.*TR.*t14.*1.0./sqrt(t17).*1.0./sqrt(-t23).*exp(-t6.*(TEp+TR.*2.0+R2p.*T2.*TEp)).*sign(t25).*(t27-1.0).*(t28-1.0).*(exp(TR.*t2.*t6.*(T2.*3.0+t26))-t5.*exp(TR.*t2.*t6.*t30.*2.0).*2.0+t27.*exp(TR.*t2.*t6.*(T2+t26))).*1.0./(exp(TR.*t2.*t6.*(T1+T2).*2.0)+t27.*t28-t5.*exp(TR.*t2.*t6.*t30).*2.0-t27.*exp(TR.*t2.*2.0)+t5.*exp(t8).*2.0-1.0).^2,M0.*TEp.*t14.*t32.*t37.*abs(t36)-M0.*t9.*t14.*t15.*t32.*1.0./sqrt(t34).*sign(t36).*(TR.*t7.*t23.*t37.*2.0-TR.*t7.*t19.*t20.*1.0./t22.^2.*t33.*t37.*2.0).*(1.0./2.0)];
