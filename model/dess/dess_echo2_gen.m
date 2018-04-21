function smxy_tem = dess_echo2_gen(M0,T1,T2,kap,dw,R2p,flip,TR,TEm)
%DESS_ECHO2_GEN
%    SMXY_TEM = DESS_ECHO2_GEN(M0,T1,T2,KAP,DW,R2P,FLIP,TR,TEM)

%    This function was generated by the Symbolic Math Toolbox version 6.1.
%    16-Jun-2016 12:45:48

t2 = 1.0./T2;
t3 = exp(TR.*t2.*-2.0);
t5 = 1.0./T1;
t6 = TR.*t5;
t7 = exp(-t6);
t8 = flip.*kap;
t9 = cos(t8);
t4 = -t7+t9;
smxy_tem = M0.*exp(TEm.*dw.*-1i).*exp(-TEm.*(R2p-t2)).*tan(flip.*kap.*(1.0./2.0)).*(sqrt((t3-1.0)./(t3.*t4.^2.*1.0./(t7.*t9-1.0).^2-1.0))-1.0).*1i;