function smxy_tem_gradx = dess_echo2_gradx_cmplx_gen(M0,T1,T2,kap,dw,R2p,flip,TR,TEm)
%DESS_ECHO2_GRADX_CMPLX_GEN
%    SMXY_TEM_GRADX = DESS_ECHO2_GRADX_CMPLX_GEN(M0,T1,T2,KAP,DW,R2P,FLIP,TR,TEM)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    06-Jan-2018 11:47:35

t2 = 1.0./T2;
t17 = TR.*t2.*2.0;
t3 = exp(-t17);
t5 = 1.0./T1;
t6 = TR.*t5;
t7 = exp(-t6);
t8 = flip.*kap;
t9 = cos(t8);
t4 = -t7+t9;
t10 = R2p.*T2;
t11 = T2.*dw.*1i;
t12 = t10+t11-1.0;
t27 = TEm.*t2.*t12;
t13 = exp(-t27);
t14 = flip.*kap.*(1.0./2.0);
t15 = tan(t14);
t16 = exp(t6);
t18 = -t3+1.0;
t19 = sqrt(t18);
t20 = t4.^2;
t21 = t7.*t9;
t22 = t21-1.0;
t23 = 1.0./t22.^2;
t24 = t3.*t20.*t23;
t25 = t24-1.0;
t26 = 1.0./t25;
t28 = sqrt(-t26);
t29 = t19.*t28;
t30 = t29-1.0;
t31 = 1.0./T2.^2;
t32 = t9-t16;
t33 = t9.^2;
t34 = t33-1.0;
t35 = 1.0./sqrt(-t26);
t36 = TR.*t5.*2.0;
t37 = exp(t36);
t38 = T1.*2.0;
smxy_tem_gradx = [t13.*t15.*t30.*1i,M0.*1.0./T1.^2.*TR.*t13.*t15.*t19.*1.0./t25.^2.*1.0./t32.^3.*t34.*t35.*exp(TR.*t2.*t5.*(T2-t38)).*(t9.*t16-1.0).*-1i,M0.*TEm.*t13.*t15.*t30.*t31.*-1i+M0.*TR.*t15.*1.0./sqrt(t18).*t31.*t32.^2.*t34.*t35.*exp(t2.*(TEm+TR.*2.0-R2p.*T2.*TEm-T2.*TEm.*dw.*1i)).*(t37-1.0).*1.0./(exp(TR.*t2.*t5.*(T1+T2).*2.0)+t9.*t16.*2.0-t33.*t37+t33.*exp(t17)-t9.*exp(TR.*t2.*t5.*(T2+t38)).*2.0-1.0).^2.*1i];