function sTE_xy_gradx = SE_shortTR_gradx_gen(M0,T1,T2,kap,TR,TE,flip_ex,flip_ref,dw)
%SE_SHORTTR_GRADX_GEN
%    STE_XY_GRADX = SE_SHORTTR_GRADX_GEN(M0,T1,T2,KAP,TR,TE,FLIP_EX,FLIP_REF,DW)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    30-May-2016 08:38:14

t2 = 1.0./T2;
t3 = TE.*t2;
t4 = exp(t3);
t5 = 1.0./T1;
t6 = TE.*t5.*(1.0./2.0);
t7 = exp(t6);
t13 = TR.*2.0;
t8 = TE-t13;
t9 = t5.*t8.*(1.0./2.0);
t10 = exp(t9);
t11 = flip_ex.*kap;
t12 = flip_ref.*kap;
t14 = sin(t11);
t15 = TE.*t2.*(1.0./2.0);
t16 = exp(t15);
t17 = TE.*dw.*(1.0./2.0);
t18 = sin(t12);
t19 = sin(t17);
t20 = cos(t12);
t21 = cos(t17);
t22 = t21.^2;
t23 = cos(t11);
t24 = T1.^2;
t25 = T1.*TE.*3.0;
t26 = T2.*TE.*2.0;
t27 = T1.*TE.*2.0;
t28 = T2.*TE.*3.0;
t29 = T1.*TE;
t37 = T2.*TR.*2.0;
t30 = t28+t29-t37;
t31 = t2.*t5.*t30.*(1.0./2.0);
t32 = exp(t31);
t33 = T1.*2.0;
t34 = T2+t33;
t35 = TE.*t2.*t5.*t34.*(1.0./2.0);
t36 = exp(t35);
t38 = T2.*TE;
t78 = T2.*TR;
t39 = t29+t38-t78;
t40 = t2.*t5.*t39;
t41 = exp(t40);
t42 = t27-t37+t38;
t43 = t2.*t5.*t42.*(1.0./2.0);
t44 = exp(t43);
t45 = t20.^2;
t46 = t26+t29-t37;
t47 = t2.*t5.*t46.*(1.0./2.0);
t48 = exp(t47);
t49 = t23.^2;
t50 = t18.^2;
t51 = t19.^2;
t52 = t14.^2;
t59 = T2.*TR.*4.0;
t53 = t28+t29-t59;
t54 = t2.*t5.*t53.*(1.0./2.0);
t55 = exp(t54);
t56 = TE.*3.0;
t57 = t13-t56;
t64 = t5.*t57.*(1.0./2.0);
t58 = exp(-t64);
t60 = TR.*4.0;
t61 = t56-t60;
t62 = t5.*t61.*(1.0./2.0);
t63 = exp(t62);
t65 = t21.*1i;
t66 = t19+t65;
t67 = T1+T2;
t68 = TE.*t2.*t5.*t67;
t69 = exp(t68);
t70 = T2.*TE.*4.0;
t71 = T2.*3.0;
t72 = T1+t71;
t73 = TE.*t2.*t5.*t72.*(1.0./2.0);
t74 = exp(t73);
t75 = t29-t37+t70;
t76 = t2.*t5.*t75.*(1.0./2.0);
t77 = exp(t76);
t79 = t27+t28-t37;
t80 = t2.*t5.*t79.*(1.0./2.0);
t81 = exp(t80);
t82 = TE.*2.0;
t83 = t13-t82;
t94 = t5.*t83;
t84 = exp(-t94);
t85 = TE.*5.0;
t86 = t60-t85;
t88 = t5.*t86.*(1.0./2.0);
t87 = exp(-t88);
t89 = t29-t59+t70;
t90 = t2.*t5.*t89.*(1.0./2.0);
t91 = exp(t90);
t92 = TR-t82;
t93 = exp(-t5.*t92);
sTE_xy_gradx = [(t10.*t14.*1i-t7.*t10.*t14.*1i-t7.*t14.*t20.*1i+t7.*t14.*t22.*1i-t10.*t14.*t22.*1i+t16.*t18.*t19.*1i-t16.*t18.*t21+t7.*t10.*t14.*t20.*1i-t7.*t16.*t18.*t19.*1i+t7.*t14.*t19.*t21+t7.*t16.*t18.*t21+t7.*t14.*t20.*t22.*1i-t10.*t14.*t19.*t21-t10.*t14.*t20.*t22.*1i-t16.*t18.*t19.*t23.*1i+t16.*t18.*t21.*t23+t7.*t14.*t19.*t20.*t21-t10.*t14.*t19.*t20.*t21+t10.*t16.*t18.*t19.*t23.*1i-t10.*t16.*t18.*t21.*t23)./(t4.*t7-t4.*t10.*t20.*t23+t7.*t10.*t14.*t16.*t18.*t19),(M0.*t66.*(TE.*t18.*t36.*1i+TE.*t14.*t21.*t32-TE.*t18.*t23.*t36.*1i-TR.*t14.*t21.*t32.*2.0+TR.*t18.*t23.*t41.*2.0i-TE.*t14.*t19.*t20.*t32.*1i-TE.*t14.*t20.*t21.*t32-TE.*t18.*t20.*t23.*t44.*1i+TE.*t14.*t19.*t32.*t45.*1i+TE.*t14.*t19.*t48.*t50.*2.0i+TE.*t18.*t20.*t44.*t49.*1i+TR.*t14.*t19.*t20.*t32.*2.0i+TR.*t14.*t20.*t21.*t32.*2.0-TR.*t18.*t20.*t23.*t41.*2.0i-TR.*t14.*t20.*t21.*t48.*2.0+TR.*t18.*t20.*t23.*t44.*2.0i-TR.*t14.*t19.*t32.*t45.*2.0i+TR.*t14.*t19.*t45.*t48.*2.0i-TR.*t18.*t20.*t44.*t49.*2.0i-TE.*t14.*t20.*t21.*t23.*t55-TE.*t14.*t19.*t23.*t48.*t50.*2.0i+TE.*t14.*t19.*t23.*t45.*t55.*1i+TE.*t14.*t21.*t23.*t45.*t55+TE.*t14.*t19.*t23.*t50.*t55.*1i+TE.*t18.*t19.*t21.*t52.*t58-TE.*t18.*t20.*t51.*t52.*t58.*1i+TE.*t18.*t45.*t51.*t52.*t63.*1i+TE.*t18.*t50.*t51.*t52.*t63.*1i+TR.*t14.*t20.*t21.*t23.*t48.*2.0-TR.*t14.*t19.*t23.*t45.*t48.*2.0i+TR.*t14.*t19.*t23.*t48.*t50.*2.0i-TR.*t18.*t19.*t21.*t52.*t58.*2.0+TR.*t18.*t20.*t51.*t52.*t58.*2.0i-TE.*t14.*t19.*t20.*t23.*t45.*t55.*1i-TE.*t14.*t19.*t20.*t23.*t50.*t55.*1i-TE.*t18.*t19.*t20.*t21.*t52.*t63).*(1.0./2.0))./(t24.*exp(TE.*t2.*t5.*(T1.*3.0+T2.*2.0).*(1.0./2.0))-t20.*t23.*t24.*exp(t2.*t5.*(t25+t26-T2.*TR.*2.0).*(1.0./2.0)).*2.0+t24.*t45.*t49.*exp(t2.*t5.*(t25+t26-T2.*TR.*4.0).*(1.0./2.0))+t14.*t18.*t19.*t24.*exp(t2.*t5.*(t27+t28-T2.*TR.*2.0).*(1.0./2.0)).*2.0+t24.*t50.*t51.*t52.*exp(t2.*t5.*(t29+t70-T2.*TR.*4.0).*(1.0./2.0))-t14.*t18.*t19.*t20.*t23.*t24.*exp(t2.*t5.*(t27+t28-T2.*TR.*4.0).*(1.0./2.0)).*2.0),M0.*1.0./T2.^2.*TE.*t66.*exp(TE.*t2.*t5.*t67.*(-1.0./2.0)).*1.0./(exp(TE.*t2.*t5.*t67.*(1.0./2.0))-t20.*t23.*exp(t2.*t5.*(t29-t37+t38).*(1.0./2.0))+t14.*t18.*t19.*exp(t5.*(TE-TR))).^2.*(t18.*exp(TE.*t2.*t5.*(t33+t71).*(1.0./2.0)).*1i-t18.*t69.*1i-t14.*t21.*t74.*2.0+t18.*t23.*t69.*1i+t14.*t21.*t77.*2.0-t18.*t23.*t81.*1i+t14.*t20.*t21.*t32.*2.0+t18.*t20.*t23.*t41.*1i-t14.*t19.*t32.*t45.*2.0i-t14.*t19.*t32.*t50.*2.0i+t14.*t19.*t20.*t74.*2.0i-t18.*t20.*t41.*t49.*1i-t14.*t19.*t20.*t77.*2.0i-t14.*t20.*t21.*t77.*2.0-t18.*t20.*t23.*t81.*1i+t14.*t19.*t45.*t77.*2.0i+t14.*t19.*t50.*t77.*2.0i+t18.*t20.*t49.*exp(t2.*t5.*(t27+t28-t59).*(1.0./2.0)).*1i+t14.*t20.*t21.*t23.*t32.*2.0-t14.*t19.*t23.*t32.*t45.*2.0i-t14.*t21.*t23.*t45.*t55.*2.0-t14.*t20.*t21.*t23.*t91.*2.0+t14.*t19.*t23.*t45.*t91.*2.0i+t14.*t21.*t23.*t45.*t91.*2.0+t18.*t19.*t21.*t52.*t87-t18.*t19.*t21.*t52.*t93-t18.*t20.*t51.*t52.*t87.*1i+t18.*t20.*t51.*t52.*t93.*1i-t18.*t45.*t51.*t52.*t84.*1i+t18.*t45.*t51.*t52.*t87.*1i-t18.*t50.*t51.*t52.*t84.*1i+t18.*t50.*t51.*t52.*t87.*1i+t14.*t19.*t20.*t23.*t45.*t55.*2.0i+t14.*t19.*t20.*t23.*t50.*t55.*2.0i-t14.*t19.*t20.*t23.*t45.*t91.*2.0i+t18.*t19.*t20.*t21.*t52.*t84-t14.*t19.*t20.*t23.*t50.*t91.*2.0i-t18.*t19.*t20.*t21.*t52.*t87).*(-1.0./2.0)];