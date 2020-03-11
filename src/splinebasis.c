#include <math.h>

double N0 (double u)
{
  if (u < 0.0e0)
    return(0.0e0);
  else if (u < 0.1e1)
    return(u);
  else if (u < 0.2e1)
    return(0.2e1 - u);
  else if (0.2e1 <= u)
    return(0.0e0);
  else
    return(0.0e0);
}


double N1 (double u)
{
  double t11;
  double t3;
  double t7;
  t3 = u * u;
  t7 = pow(u - 0.1e1, 0.2e1);
  t11 = pow(u - 0.2e1, 0.2e1);
  if (u < 0.0e0)
    return(0.0e0);
  else if (u < 0.1e1)
    return(t3 / 0.2e1);
  else if (u < 0.2e1)
    return(-0.1e1 / 0.2e1 + u - t7);
  else if (u < 0.3e1)
    return(0.5e1 / 0.2e1 - u + t11 / 0.2e1);
  else if (0.3e1 <= u)
    return(0.0e0);
  else
    return(0.0e0);
}

double N2 (double u)
{
  double t15;
  double t16;
  double t21;
  double t22;
  double t3;
  double t7;
  double t8;
  double t9;
  t3 = u * u;
  t7 = u / 0.2e1;
  t8 = u - 0.1e1;
  t9 = t8 * t8;
  t15 = u - 0.2e1;
  t16 = t15 * t15;
  t21 = u - 0.3e1;
  t22 = t21 * t21;
  if (u < 0.0e0)
    return(0.0e0);
  else if (u < 0.1e1)
    return(t3 * u / 0.6e1);
  else if (u < 0.2e1)
    return(-0.1e1 / 0.3e1 + t7 + t9 / 0.2e1 - t9 * t8 / 0.2e1);
  else if (u < 0.3e1)
    return(0.2e1 / 0.3e1 - t16 + t16 * t15 / 0.2e1);
  else if (u < 0.4e1)
    return(0.5e1 / 0.3e1 - t7 + t22 / 0.2e1 - t22 * t21 / 0.6e1);
  else if (0.4e1 <= u)
    return(0.0e0);
  else
    return(0.0e0);
}
double N3 (double u)
{
  double t13;
  double t17;
  double t18;
  double t19;
  double t23;
  double t27;
  double t28;
  double t3;
  double t32;
  double t36;
  double t37;
  double t4;
  double t41;
  double t7;
  double t8;
  double t9;
  t3 = u * u;
  t4 = t3 * t3;
  t7 = u / 0.6e1;
  t8 = u - 0.1e1;
  t9 = t8 * t8;
  t13 = t9 * t9;
  t17 = u / 0.2e1;
  t18 = u - 0.2e1;
  t19 = t18 * t18;
  t23 = t19 * t19;
  t27 = u - 0.3e1;
  t28 = t27 * t27;
  t32 = t28 * t28;
  t36 = u - 0.4e1;
  t37 = t36 * t36;
  t41 = t37 * t37;
  if (u < 0.0e0)
    return(0.0e0);
  else if (u < 0.1e1)
    return(t4 / 0.24e2);
  else if (u < 0.2e1)
    return(-0.1e1 / 0.8e1 + t7 + t9 / 0.4e1 + t9 * t8 / 0.6e1 - t13 / 0.6e1);
  else if (u < 0.3e1)
    return(-0.13e2 / 0.24e2 + t17 - t19 / 0.4e1 - t19 * t18 / 0.2e1 + t23 / 0.4e1);
  else if (u < 0.4e1)
    return(0.47e2 / 0.24e2 - t17 - t28 / 0.4e1 + t28 * t27 / 0.2e1 - t32 / 0.6e1);
  else if (u < 0.5e1)
    return(0.17e2 / 0.24e2 - t7 + t37 / 0.4e1 - t36 * t37 / 0.6e1 + t41 / 0.24e2);
  else if (0.5e1 <= u)
    return(0.0e0);
  else
    return(0.0e0);
}
double N4 (double u)
{
  double t10;
  double t14;
  double t20;
  double t21;
  double t22;
  double t26;
  double t3;
  double t32;
  double t33;
  double t35;
  double t4;
  double t41;
  double t42;
  double t46;
  double t52;
  double t53;
  double t57;
  double t8;
  double t9;
  t3 = u * u;
  t4 = t3 * t3;
  t8 = u / 0.24e2;
  t9 = u - 0.1e1;
  t10 = t9 * t9;
  t14 = t10 * t10;
  t20 = 0.5e1 / 0.12e2 * u;
  t21 = u - 0.2e1;
  t22 = t21 * t21;
  t26 = t22 * t22;
  t32 = u - 0.3e1;
  t33 = t32 * t32;
  t35 = t33 * t33;
  t41 = u - 0.4e1;
  t42 = t41 * t41;
  t46 = t42 * t42;
  t52 = u - 0.5e1;
  t53 = t52 * t52;
  t57 = t53 * t53;
  if (u < 0.0e0)
    return(0.0e0);
  else if (u < 0.1e1)
    return(t4 * u / 0.120e3);
  else if (u < 0.2e1)
    return(-0.1e1 / 0.30e2 + t8 + t10 / 0.12e2 + t10 * t9 / 0.12e2 + t14 / 0.24e2 - t14 * t9 / 0.24e2);
  else if (u < 0.3e1)
    return(-0.37e2 / 0.60e2 + t20 + t22 / 0.6e1 - t22 * t21 / 0.6e1 - t26 / 0.6e1 + t26 * t21 / 0.12e2);
  else if (u < 0.4e1)
    return(0.11e2 / 0.20e2 - t33 / 0.2e1 + t35 / 0.4e1 - t35 * t32 / 0.12e2);
  else if (u < 0.5e1)
    return(0.113e3 / 0.60e2 - t20 + t42 / 0.6e1 + t42 * t41 / 0.6e1 - t46 / 0.6e1 + t46 * t41 / 0.24e2);
  else if (u < 0.6e1)
    return(0.13e2 / 0.60e2 - t8 + t53 / 0.12e2 - t53 * t52 / 0.12e2 + t57 / 0.24e2 - t57 * t52 / 0.120e3);
  else if (0.6e1 <= u)
    return(0.0e0);
  else
    return(0.0e0);
}
double N5 (double u)
{
  double t10;
  double t14;
  double t22;
  double t23;
  double t24;
  double t28;
  double t3;
  double t36;
  double t37;
  double t38;
  double t4;
  double t42;
  double t50;
  double t51;
  double t55;
  double t63;
  double t64;
  double t68;
  double t76;
  double t77;
  double t8;
  double t81;
  double t9;
  t3 = u * u;
  t4 = t3 * t3;
  t8 = u / 0.120e3;
  t9 = u - 0.1e1;
  t10 = t9 * t9;
  t14 = t10 * t10;
  t22 = 0.5e1 / 0.24e2 * u;
  t23 = u - 0.2e1;
  t24 = t23 * t23;
  t28 = t24 * t24;
  t36 = u / 0.3e1;
  t37 = u - 0.3e1;
  t38 = t37 * t37;
  t42 = t38 * t38;
  t50 = u - 0.4e1;
  t51 = t50 * t50;
  t55 = t51 * t51;
  t63 = u - 0.5e1;
  t64 = t63 * t63;
  t68 = t64 * t64;
  t76 = u - 0.6e1;
  t77 = t76 * t76;
  t81 = t77 * t77;
  if (u < 0.0e0)
    return(0.0e0);
  else if (u < 0.1e1)
    return(t4 * t3 / 0.720e3);
  else if (u < 0.2e1)
    return(-0.1e1 / 0.144e3 + t8 + t10 / 0.48e2 + t10 * t9 / 0.36e2 + t14 / 0.48e2 + t14 * t9 / 0.120e3 - t14 * t10 / 0.120e3);
  else if (u < 0.3e1)
    return(-0.27e2 / 0.80e2 + t22 + 0.3e1 / 0.16e2 * t24 + t24 * t23 / 0.36e2 - t28 / 0.16e2 - t28 * t23 / 0.24e2 + t28 * t24 / 0.48e2);
  else if (u < 0.4e1)
    return(-0.209e3 / 0.360e3 + t36 - 0.5e1 / 0.24e2 * t38 - 0.2e1 / 0.9e1 * t38 * t37 + t42 / 0.24e2 + t42 * t37 / 0.12e2 - t42 * t38 / 0.36e2);
  else if (u < 0.5e1)
    return(0.631e3 / 0.360e3 - t36 - 0.5e1 / 0.24e2 * t51 + 0.2e1 / 0.9e1 * t51 * t50 + t55 / 0.24e2 - t55 * t50 / 0.12e2 + t55 * t51 / 0.48e2);
  else if (u < 0.6e1)
    return(0.269e3 / 0.240e3 - t22 + 0.3e1 / 0.16e2 * t64 - t64 * t63 / 0.36e2 - t68 / 0.16e2 + t68 * t63 / 0.24e2 - t68 * t64 / 0.120e3);
  else if (u < 0.7e1)
    return(0.37e2 / 0.720e3 - t8 + t77 / 0.48e2 - t77 * t76 / 0.36e2 + t81 / 0.48e2 - t81 * t76 / 0.120e3 + t81 * t77 / 0.720e3);
  else if (0.7e1 <= u)
    return(0.0e0);
  else
    return(0.0e0);
}
double N6 (double u)
{
  double t10;
  double t100;
  double t102;
  double t11;
  double t13;
  double t15;
  double t25;
  double t26;
  double t27;
  double t29;
  double t3;
  double t31;
  double t40;
  double t41;
  double t42;
  double t44;
  double t46;
  double t5;
  double t56;
  double t57;
  double t59;
  double t68;
  double t69;
  double t71;
  double t73;
  double t83;
  double t84;
  double t86;
  double t88;
  double t9;
  double t97;
  double t98;
  t3 = u * u;
  t5 = t3 * t3;
  t9 = u / 0.720e3;
  t10 = u - 0.1e1;
  t11 = t10 * t10;
  t13 = t11 * t10;
  t15 = t11 * t11;
  t25 = 0.7e1 / 0.90e2 * u;
  t26 = u - 0.2e1;
  t27 = t26 * t26;
  t29 = t27 * t26;
  t31 = t27 * t27;
  t40 = 0.49e2 / 0.144e3 * u;
  t41 = u - 0.3e1;
  t42 = t41 * t41;
  t44 = t42 * t41;
  t46 = t42 * t42;
  t56 = u - 0.4e1;
  t57 = t56 * t56;
  t59 = t57 * t57;
  t68 = u - 0.5e1;
  t69 = t68 * t68;
  t71 = t69 * t68;
  t73 = t69 * t69;
  t83 = u - 0.6e1;
  t84 = t83 * t83;
  t86 = t84 * t83;
  t88 = t84 * t84;
  t97 = u - 0.7e1;
  t98 = t97 * t97;
  t100 = t98 * t97;
  t102 = t98 * t98;
  if (u < 0.0e0)
    return(0.0e0);
  else if (u < 0.1e1)
    return(t5 * t3 * u / 0.5040e4);
  else if (u < 0.2e1)
    return(-0.1e1 / 0.840e3 + t9 + t11 / 0.240e3 + t13 / 0.144e3 + t15 / 0.144e3 + t15 * t10 / 0.240e3 + t15 * t11 / 0.720e3 - t15 * t13 / 0.720e3);
  else if (u < 0.3e1)
    return(-0.83e2 / 0.630e3 + t25 + t27 / 0.10e2 + t29 / 0.18e2 - t31 * t26 / 0.60e2 - t31 * t27 / 0.120e3 + t31 * t29 / 0.240e3);
  else if (u < 0.4e1)
    return(-0.659e3 / 0.840e3 + t40 + t42 / 0.16e2 - 0.19e2 / 0.144e3 * t44 - t46 / 0.16e2 + t46 * t41 / 0.48e2 + t46 * t42 / 0.48e2 - t46 * t44 / 0.144e3);
  else if (u < 0.5e1)
    return(0.151e3 / 0.315e3 - t57 / 0.3e1 + t59 / 0.9e1 - t59 * t57 / 0.36e2 + t59 * t57 * t56 / 0.144e3);
  else if (u < 0.6e1)
    return(0.4883e4 / 0.2520e4 - t40 + t69 / 0.16e2 + 0.19e2 / 0.144e3 * t71 - t73 / 0.16e2 - t73 * t68 / 0.48e2 + t73 * t69 / 0.48e2 - t73 * t71 / 0.240e3);
  else if (u < 0.7e1)
    return(0.103e3 / 0.210e3 - t25 + t84 / 0.10e2 - t86 / 0.18e2 + t88 * t83 / 0.60e2 - t88 * t84 / 0.120e3 + t88 * t86 / 0.720e3);
  else if (u < 0.8e1)
    return(0.5e1 / 0.504e3 - t9 + t98 / 0.240e3 - t100 / 0.144e3 + t102 / 0.144e3 - t102 * t97 / 0.240e3 + t102 * t98 / 0.720e3 - t100 * t102 / 0.5040e4);
  else if (0.8e1 <= u)
    return(0.0e0);
  else
    return(0.0e0);
}


int dN0 (double u)
{
  if (u < 0.0e0)
    return(0);
  else if (u < 0.1e1)
    return(1);
  else if (u < 0.2e1)
    return(-1);
  else if (0.2e1 < u)
    return(0);
  else
    return(0);
}


double dN1 (double u)
{
  if (u <= 0.0e0)
    return(0.0e0);
  else if (u <= 0.1e1)
    return(u);
  else if (u <= 0.2e1)
    return(0.3e1 - 0.2e1 * u);
  else if (u <= 0.3e1)
    return(-0.3e1 + u);
  else if (0.3e1 < u)
    return(0.0e0);
  else
    return(0.0e0);
}


double dN2 (double u)
{
  double t13;
  double t18;
  double t3;
  double t7;
  t3 = u * u;
  t7 = pow(u - 0.1e1, 0.2e1);
  t13 = pow(-0.2e1 + u, 0.2e1);
  t18 = pow(-0.3e1 + u, 0.2e1);
  if (u <= 0.0e0)
    return(0.0e0);
  else if (u <= 0.1e1)
    return(t3 / 0.2e1);
  else if (u <= 0.2e1)
    return(-0.1e1 / 0.2e1 + u - 0.3e1 / 0.2e1 * t7);
  else if (u <= 0.3e1)
    return(0.4e1 - 0.2e1 * u + 0.3e1 / 0.2e1 * t13);
  else if (u <= 0.4e1)
    return(-0.7e1 / 0.2e1 + u - t18 / 0.2e1);
  else if (0.4e1 < u)
    return(0.0e0);
  else
    return(0.0e0);
}


double dN3 (double u)
{
  double t15;
  double t16;
  double t21;
  double t22;
  double t28;
  double t29;
  double t3;
  double t7;
  double t8;
  double t9;
  t3 = u * u;
  t7 = u / 0.2e1;
  t8 = u - 0.1e1;
  t9 = t8 * t8;
  t15 = -0.2e1 + u;
  t16 = t15 * t15;
  t21 = -0.3e1 + u;
  t22 = t21 * t21;
  t28 = -0.4e1 + u;
  t29 = t28 * t28;
  if (u <= 0.0e0)
    return(0.0e0);
  else if (u <= 0.1e1)
    return(t3 * u / 0.6e1);
  else if (u <= 0.2e1)
    return(-0.1e1 / 0.3e1 + t7 + t9 / 0.2e1 - 0.2e1 / 0.3e1 * t9 * t8);
  else if (u <= 0.3e1)
    return(0.3e1 / 0.2e1 - t7 - 0.3e1 / 0.2e1 * t16 + t16 * t15);
  else if (u <= 0.4e1)
    return(0.1e1 - t7 + 0.3e1 / 0.2e1 * t22 - 0.2e1 / 0.3e1 * t22 * t21);
  else if (u <= 0.5e1)
    return(-0.13e2 / 0.6e1 + t7 - t29 / 0.2e1 + t29 * t28 / 0.6e1);
  else if (0.5e1 < u)
    return(0.0e0);
  else
    return(0.0e0);
}


double dN4 (double u)
{
  double t13;
  double t17;
  double t18;
  double t19;
  double t23;
  double t27;
  double t28;
  double t3;
  double t30;
  double t34;
  double t35;
  double t39;
  double t4;
  double t43;
  double t44;
  double t48;
  double t7;
  double t8;
  double t9;
  t3 = u * u;
  t4 = t3 * t3;
  t7 = u / 0.6e1;
  t8 = u - 0.1e1;
  t9 = t8 * t8;
  t13 = t9 * t9;
  t17 = u / 0.3e1;
  t18 = -0.2e1 + u;
  t19 = t18 * t18;
  t23 = t19 * t19;
  t27 = -0.3e1 + u;
  t28 = t27 * t27;
  t30 = t28 * t28;
  t34 = -0.4e1 + u;
  t35 = t34 * t34;
  t39 = t35 * t35;
  t43 = -0.5e1 + u;
  t44 = t43 * t43;
  t48 = t44 * t44;
  if (u <= 0.0e0)
    return(0.0e0);
  else if (u <= 0.1e1)
    return(t4 / 0.24e2);
  else if (u <= 0.2e1)
    return(-0.1e1 / 0.8e1 + t7 + t9 / 0.4e1 + t9 * t8 / 0.6e1 - 0.5e1 / 0.24e2 * t13);
  else if (u <= 0.3e1)
    return(-0.1e1 / 0.4e1 + t17 - t19 / 0.2e1 - 0.2e1 / 0.3e1 * t19 * t18 + 0.5e1 / 0.12e2 * t23);
  else if (u <= 0.4e1)
    return(0.3e1 - u + t28 * t27 - 0.5e1 / 0.12e2 * t30);
  else if (u <= 0.5e1)
    return(-0.7e1 / 0.4e1 + t17 + t35 / 0.2e1 - 0.2e1 / 0.3e1 * t35 * t34 + 0.5e1 / 0.24e2 * t39);
  else if (u <= 0.6e1)
    return(-0.7e1 / 0.8e1 + t7 - t44 / 0.4e1 + t44 * t43 / 0.6e1 - t48 / 0.24e2);
  else if (0.6e1 < u)
    return(0.0e0);
  else
    return(0.0e0);
}


double dN5 (double u)
{
  double t10;
  double t14;
  double t20;
  double t21;
  double t22;
  double t26;
  double t3;
  double t32;
  double t33;
  double t34;
  double t38;
  double t4;
  double t44;
  double t45;
  double t49;
  double t55;
  double t56;
  double t60;
  double t66;
  double t67;
  double t71;
  double t8;
  double t9;
  t3 = u * u;
  t4 = t3 * t3;
  t8 = u / 0.24e2;
  t9 = u - 0.1e1;
  t10 = t9 * t9;
  t14 = t10 * t10;
  t20 = 0.3e1 / 0.8e1 * u;
  t21 = -0.2e1 + u;
  t22 = t21 * t21;
  t26 = t22 * t22;
  t32 = 0.5e1 / 0.12e2 * u;
  t33 = -0.3e1 + u;
  t34 = t33 * t33;
  t38 = t34 * t34;
  t44 = -0.4e1 + u;
  t45 = t44 * t44;
  t49 = t45 * t45;
  t55 = -0.5e1 + u;
  t56 = t55 * t55;
  t60 = t56 * t56;
  t66 = -0.6e1 + u;
  t67 = t66 * t66;
  t71 = t67 * t67;
  if (u <= 0.0e0)
    return(0.0e0);
  else if (u <= 0.1e1)
    return(t4 * u / 0.120e3);
  else if (u <= 0.2e1)
    return(-0.1e1 / 0.30e2 + t8 + t10 / 0.12e2 + t10 * t9 / 0.12e2 + t14 / 0.24e2 - t14 * t9 / 0.20e2);
  else if (u <= 0.3e1)
    return(-0.13e2 / 0.24e2 + t20 + t22 / 0.12e2 - t22 * t21 / 0.4e1 - 0.5e1 / 0.24e2 * t26 + t26 * t21 / 0.8e1);
  else if (u <= 0.4e1)
    return(0.19e2 / 0.12e2 - t32 - 0.2e1 / 0.3e1 * t34 + t34 * t33 / 0.6e1 + 0.5e1 / 0.12e2 * t38 - t38 * t33 / 0.6e1);
  else if (u <= 0.5e1)
    return(0.4e1 / 0.3e1 - t32 + 0.2e1 / 0.3e1 * t45 + t45 * t44 / 0.6e1 - 0.5e1 / 0.12e2 * t49 + t49 * t44 / 0.8e1);
  else if (u <= 0.6e1)
    return(-0.25e2 / 0.12e2 + t20 - t56 / 0.12e2 - t56 * t55 / 0.4e1 + 0.5e1 / 0.24e2 * t60 - t60 * t55 / 0.20e2);
  else if (u <= 0.7e1)
    return(-0.31e2 / 0.120e3 + t8 - t67 / 0.12e2 + t67 * t66 / 0.12e2 - t71 / 0.24e2 + t71 * t66 / 0.120e3);
  else if (0.7e1 < u)
    return(0.0e0);
  else
    return(0.0e0);
}


double dN6 (double u)
{
  double t10;
  double t14;
  double t22;
  double t23;
  double t24;
  double t26;
  double t3;
  double t34;
  double t35;
  double t36;
  double t4;
  double t40;
  double t49;
  double t50;
  double t53;
  double t60;
  double t61;
  double t65;
  double t73;
  double t74;
  double t76;
  double t8;
  double t84;
  double t85;
  double t89;
  double t9;
  t3 = u * u;
  t4 = t3 * t3;
  t8 = u / 0.120e3;
  t9 = u - 0.1e1;
  t10 = t9 * t9;
  t14 = t10 * t10;
  t22 = u / 0.5e1;
  t23 = -0.2e1 + u;
  t24 = t23 * t23;
  t26 = t24 * t24;
  t34 = u / 0.8e1;
  t35 = -0.3e1 + u;
  t36 = t35 * t35;
  t40 = t36 * t36;
  t49 = -0.4e1 + u;
  t50 = t49 * t49;
  t53 = t50 * t50;
  t60 = -0.5e1 + u;
  t61 = t60 * t60;
  t65 = t61 * t61;
  t73 = -0.6e1 + u;
  t74 = t73 * t73;
  t76 = t74 * t74;
  t84 = -0.7e1 + u;
  t85 = t84 * t84;
  t89 = t85 * t85;
  if (u <= 0.0e0)
    return(0.0e0);
  else if (u <= 0.1e1)
    return(t4 * t3 / 0.720e3);
  else if (u <= 0.2e1)
    return(-0.1e1 / 0.144e3 + t8 + t10 / 0.48e2 + t10 * t9 / 0.36e2 + t14 / 0.48e2 + t14 * t9 / 0.120e3 - 0.7e1 / 0.720e3 * t14 * t10);
  else if (u <= 0.3e1)
    return(-0.29e2 / 0.90e2 + t22 + t24 / 0.6e1 - t26 / 0.12e2 - t26 * t23 / 0.20e2 + 0.7e1 / 0.240e3 * t26 * t24);
  else if (u <= 0.4e1)
    return(-0.5e1 / 0.144e3 + t34 - 0.19e2 / 0.48e2 * t36 - t36 * t35 / 0.4e1 + 0.5e1 / 0.48e2 * t40 + t40 * t35 / 0.8e1 - 0.7e1 / 0.144e3 * t40 * t36);
  else if (u <= 0.5e1)
    return(0.8e1 / 0.3e1 - 0.2e1 / 0.3e1 * u + 0.4e1 / 0.9e1 * t50 * t49 - t53 * t49 / 0.6e1 + 0.7e1 / 0.144e3 * t53 * t50);
  else if (u <= 0.6e1)
    return(-0.139e3 / 0.144e3 + t34 + 0.19e2 / 0.48e2 * t61 - t61 * t60 / 0.4e1 - 0.5e1 / 0.48e2 * t65 + t65 * t60 / 0.8e1 - 0.7e1 / 0.240e3 * t61 * t65);
  else if (u <= 0.7e1)
    return(-0.23e2 / 0.18e2 + t22 - t74 / 0.6e1 + t76 / 0.12e2 - t76 * t73 / 0.20e2 + 0.7e1 / 0.720e3 * t76 * t74);
  else if (u <= 0.8e1)
    return(-0.43e2 / 0.720e3 + t8 - t85 / 0.48e2 + t85 * t84 / 0.36e2 - t89 / 0.48e2 + t89 * t84 / 0.120e3 - t89 * t85 / 0.720e3);
  else if (0.8e1 < u)
    return(0.0e0);
  else
    return(0.0e0);
}
