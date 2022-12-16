#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
 *    This file is part of TEMPO2. 
 * 
 *    TEMPO2 is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    TEMPO2 is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing 
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tempo2.h"
#include "DDSTGmodel.h"


/* Timing model      */
/* Based on bnryddgr.f */

/* Damour & Deruelle model assuming the general theory of relativity */
/* Computes the pulsar orbit time, torb, at the time of observation, t =ct(n)-pepoch */
/* Pulsar proper time is TP = T + TORB */
/* Units are c=g=1. */

void mass2ddstg(double am,double am2,double x,double ecc,double an,double *arr,double *ar,
        double *xk,double *si,double *gamma,double *pbdot,double *dr,double *dth, DEF_factors *factors, DEF_PPK_derivatives *derivatives);


double DDSTGmodel(pulsar *psr,int p,int ipos,int param)
{
    double an,afac;
    double pb,k;
    double rad2deg = 180.0/M_PI;
    double SUNMASS = 4.925490947e-6, BARE_SUNMASS;
    double m2,tt0,t0,x,ecc,er,xdot,edot,dr,dth,eth,am2,ct,f0,cm;
    double pbdot,xpbdot,phase,u,du,gamma,m,m1,arr,ar,xk,fac,xomdot;
    double dtdm2,dgmdm2,dsidm2,dkdm2,dthdm2,dpbdm2,darrdm2,dnum,denom;
    double orbits,a0aligned;
    int norbits;
    double  cu,onemecu,ae,omega,omz,sw,cw,alpha,beta,bg,dre,drep,drepp,anhat,su;
    double sqr1me2,cume,brace,si,dlogbr,ds,da,a0,b0,d2bar,torb,am;
    double csigma,ce,cx,comega,cgamma,cdth,cm2,csi;
    double fact1,fact2,fact3,fact4,fact5,fact6,fact7,fact8,denumm,denomm,darrdm,ck,dkdm,cdr;
    double ddrdm,cpbdot,dpbdm,csini,dsidm,an0,dgamdm,dthdm,ddrdm2;
    const char *CVS_verNum = "$Id$";

    DEF_tables tables;
    DEF_factors factors;
    DEF_PPK_derivatives derivatives;
    double alpha0, beta0, mA, alphaA, betaA, kA, mB, alphaB, betaB, kB;
    double bare_m, bare_m1, bare_m2;
    double dk_dm, dk_dm2, dgamma_dm, dgamma_dm2, dsi_dm, dsi_dm2, ddr_dm, ddr_dm2, ddth_dm, ddth_dm2, dpbdot_dm, dpbdot_dm2, ddt_dm, ddt_dm2;
    double epsNum, tt, frb, delta, delta_old, diff, ae1, cae, sae, psi, cpsi, spsi, dRoe, dEin, dSha, dAbe;
    int loop_counter1, loop_counter2;

    if (displayCVSversion == 1) CVSdisplayVersion("DDSTGmodel.C","DDSTGmodel()",CVS_verNum);

    t0 = psr[p].param[param_t0].val[0];
    ct = psr[p].obsn[ipos].bbat;    

    tt0 = (ct-t0)*SECDAY;

    f0 = psr[p].param[param_f].val[0];

    xomdot = 0.0;  /* WHAT SHOULD THIS BE??? */
    afac = 0.0;    /* WHAT SHOULD THIS BE??? */

    if (psr[p].param[param_sini].paramSet[0]==1) si = getParameterValue(&psr[p],param_sini,0);
    else si = 0.0;

    if (psr[p].param[param_m2].paramSet[0]==1) am2 = psr[p].param[param_m2].val[0];
    else am2 = 0.0;

    pb = psr[p].param[param_pb].val[0]*SECDAY;
    frb = 1.0/pb;
    an = 2.0*M_PI/pb;

    if (psr[p].param[param_mtot].paramSet[0]==1) am = psr[p].param[param_mtot].val[0];
    else am = 0.0;

    m  = am*SUNMASS;
    m2 = am2*SUNMASS;
    m1 = m-m2;

    if (psr[p].param[param_om].paramSet[0]==1) omz = psr[p].param[param_om].val[0];
    else omz = 0.0;

    if (psr[p].param[param_a1dot].paramSet[0]==1) xdot  = psr[p].param[param_a1dot].val[0];
    else xdot  = 0.0;

    if (psr[p].param[param_pbdot].paramSet[0] == 1) pbdot = psr[p].param[param_pbdot].val[0];
    else pbdot = 0.0;

    if (psr[p].param[param_edot].paramSet[0] == 1) edot = psr[p].param[param_edot].val[0];
    else edot = 0.0;

    if (psr[p].param[param_xpbdot].paramSet[0] == 1) xpbdot = psr[p].param[param_xpbdot].val[0];
    else xpbdot = 0.0;

    x = psr[p].param[param_a1].val[0]+xdot*tt0;
    ecc = psr[p].param[param_ecc].val[0]+edot*tt0;

    /* Given system masses m,m2 interpolate DDSTG parameters*/

    factors.alpha0 = psr->param[param_alpha0].val[0];
    factors.beta0 = psr->param[param_beta0].val[0];
    factors.companion_type = psr->companion_type;
    factors.mA = am - am2;
    factors.mB = am2;

    alpha0 = factors.alpha0;
    beta0 = factors.betaA;
    mA = factors.mA;
    mB = factors.mB;

    BARE_SUNMASS = SUNMASS/(1+pow(alpha0,2));
    bare_m = am*BARE_SUNMASS;
    bare_m2 = am2*BARE_SUNMASS;
    bare_m1 = bare_m-bare_m2;

    if (psr->ddstg_init == 0)
    {
        puts("DDSTG tables are not loaded. Starting loading.");
        psr->ddstg_init = 1;
        tables.eosName = psr->eos_name;
        printf("EOS: %s\n", tables.eosName);
        printf("%s\n", getenv(TEMPO2_ENVIRON));
        DEF_tables_load(&tables);
        allocate_DEF_arrays(psr, &tables);
        interpolate_DEF_arrays(psr, &tables);
        DEF_tables_free(&tables);
        puts("DDSTG tables were loaded and arrays were interpolated.");
    //    interpolate_DEF_params(psr, mA, &alphaA, &betaA, &kA);
    //    printf("mA, alphaA, betaA, kA: %lf %lf %lf %lf\n", mA, alphaA, betaA, kA);

        interpolate_DEF_factors(psr, &factors);
        printf("alpha0, beta0: %lf %lf\n", factors.alpha0, factors.beta0);
        printf("mA, alphaA, betaA, kA: %lf %lf %lf %lf\n", factors.mA, factors.alphaA, factors.betaA, factors.kA);
        printf("Companion type: %s\n", factors.companion_type);
        printf("mB, alphaB, betaB, kB: %lf %lf %lf %lf\n", factors.mB, factors.alphaB, factors.betaB, factors.kB);
    }

    interpolate_DEF_factors(psr, &factors);

    alphaA = factors.alphaA;
    betaA = factors.betaA;
    kA = factors.kA;
    alphaB = factors.alphaB;
    betaB = factors.betaB;
    kB = factors.kB;



    /* Given system masses m,m2 and keplerian parameters x,ecc,an, calculate the values
     * of arr,ar,si,gamma,pbdot under general relativity */

    mass2ddstg(am,am2,x,ecc,an,&arr,&ar,&xk,&si,&gamma,&pbdot,&dr,&dth, &factors, &derivatives);

    k=xk;

//    printf("dr value in DDSTG vs DDGR: %.17g ", dr);
//    dr = (3.0*pow(m1,2) + 6.0*m1*m2 + 2.0*pow(m2,2))/(arr*m); //GR
//    printf("%.17g\n", dr);
    er = ecc*(1.0+dr);
//    printf("dth value in DDSTG vs DDGR: %.17g ", dth);
//    dth = (3.5*m1*m1 + 6*m1*m2 + 2*m2*m2)/(arr*m); //GR
//    printf("%.17g\n", dth);
    eth = ecc*(1.0+dth);

//----------------------------------------------------------------------------------------
/*
    orbits = tt0/pb - 0.5*(pbdot+xpbdot)*pow(tt0/pb,2);
    //  printf("xpbdot = %.14g %.14g\n",(double)xpbdot/1e-12,(double)orbits);
    norbits = (int)orbits;
    if (orbits<0.0) norbits--;
    phase=2.0*M_PI*(orbits-norbits);
    //  Compute eccentric anomaly u by iterating Kepler's equation.
    u=phase+ecc*sin(phase)*(1.0+ecc*cos(phase));

    do {
        
    	fac = 1.0/(1.0-ecc*cos(u));  // NOTE COULD BE WRONG IN DDmodel - SEE USE OF FAC !!!!
	du=(phase-(u-ecc*sin(u)))*fac; 
        u=u+du;
   	
    } while (fabs(du)>1.0e-14);  // 1e-12 in DDmodel

    //  DD equations 17a, 29
    ae = 2.0*atan(sqrt((1+ecc)/(1-ecc))*tan(0.5*u));
    if(ae<0.0) ae=ae+2.0*M_PI;
    ae = 2.0*M_PI*orbits + ae-phase;
    omega=omz/rad2deg + (k+xomdot/(an*rad2deg*365.25*86400.0))*ae;
    // DD equations 46 through 52
    su=sin(u);
    cu=cos(u);
    sw=sin(omega);
    cw=cos(omega);
    alpha=x*sw;
    beta=x*sqrt(1-pow(eth,2))*cw;
    bg=beta+gamma;
    dre=alpha*(cu-er) + bg*su;
    drep=-alpha*su + bg*cu;
    drepp=-alpha*cu - bg*su;
    onemecu=1.0-ecc*cu;
    anhat=an/onemecu;

    // DD equations 26,27,57

    cume=cu-ecc;
    sqr1me2=sqrt(1-pow(ecc,2));
    brace=onemecu-si*(sw*cume+sqr1me2*cw*su);
    if (brace<=0)
    {
        printf("ERROR: In DDSTG model, brace < 0\n");
        exit(1);
    }
    dlogbr=log(brace);
    ds=-2*bare_m2*dlogbr; //change m2 to bare_m2

    // These will be different if spin axis not aligned -- IS THIS AN ASSUMPTION OF THE MODEL?
    a0aligned = an*ar/(2.0*M_PI*f0*si*sqr1me2);
    a0 = afac*a0aligned;
    b0 = 0.0;
    da = a0*(sin(omega+ae)+ecc*sw) + b0*(cos(omega+ae) + ecc*cw);


    //  Now compute d2bar, the orbital time correction in DD equation 42.
    d2bar=dre*(1-anhat*drep+(pow(anhat,2))*(pow(drep,2) + 0.5*dre*drepp -
                +    0.5*ecc*su*dre*drep/onemecu)) + ds + da;
    torb=-d2bar;


    printf("Comparison inversion DDGR vs DDSTG");
    printf("DDGR ae: \n%.17g\n", ae);
    printf("DDGR omega: \n%.17g\n", omega);
    printf("DDGR torb: \n%.17g\n", torb);
*/
//----------------------------------------------------------------------------------------

    epsNum = 1.0e-10;
    delta  = 0.0;
    loop_counter1 = 0;
    do
    {   
        loop_counter1 += 1;
        delta_old = delta;        
        tt  = tt0 - delta;
        orbits = tt*frb - 0.5*(pbdot+xpbdot)*pow(tt*frb,2);
        norbits = (int)orbits;
        if (orbits < 0.0) norbits--;
        phase = 2.0*M_PI*(orbits - norbits);

//  Compute eccentric anomaly u by iterating Keplers equation.
        u=phase+ecc*sin(phase)*(1+ecc*cos(phase));
        loop_counter2 = 0;
        do {
            loop_counter2 += 1;
            du=(phase-(u-ecc*sin(u)))/(1.0-ecc*cos(u));
            u=u+du;
        } while (fabs(du)>1.0e-14 && loop_counter2 < 30);  // 1e-12 in DDmodel

//  DD equations 17b, 17c, 29, and 46 through 52
        su=sin(u);
        cu=cos(u);
        onemecu=1.0-ecc*cu;
        cume = cu - ecc;

        cae=(cu-ecc)/onemecu;
        sae=sqrt(1.0-pow(ecc,2))*su/onemecu;
        ae1=atan2(sae,cae);
        if(ae1 < 0.0) ae1=ae1+2.0*M_PI;
        ae=ae1 + 2.0*M_PI*norbits;

        omega=omz/rad2deg + (k+xomdot/(an*rad2deg*365.25*86400.0))*ae;

        sw=sin(omega);
        cw=cos(omega);

        psi  = omega + ae1; // angle w.r.t. ascending node
        spsi = sin(psi);
        cpsi = cos(psi);

//  Roemer delay (DD)
        alpha = x*sw;
        beta  = x*sqrt(1.0 - pow(eth,2))*cw;
        dRoe  = alpha*(cu - er) + beta*su;

//  Einstein delay (DD)
        dEin  = gamma * su;

//  Shapiro delay

        sqr1me2 = sqrt(1.0 - pow(ecc,2));

        brace  = onemecu - si*(sw*cume + sqr1me2*cw*su);
        if (brace<=0)
        {
            printf("ERROR: In DDSTG model, brace < 0\n");
            exit(1);
        }
        dlogbr = log(brace);
        dSha   = -2.0*bare_m2*dlogbr;                                   //change m2 to bare_m2

//  Aberration
        a0aligned=an*ar/(2.0*M_PI*f0*si*sqr1me2);                         //added a0 and b0
        a0=afac*a0aligned;
        b0=0.0;

        dAbe = a0*(spsi+ecc*sw)+b0*(cpsi+ecc*cw);

        delta = dRoe + dEin + dSha + dAbe;

        diff  = fabs(delta - delta_old);
    } while( (diff > epsNum) && loop_counter1 < 30 );
//  Inversion of timing model by iteration: end of loop

/*
    printf("\nDDSTG ae: \n%.17g\n", ae);
    printf("DDSTG omega: \n%.17g\n", omega);
    printf("DDSTG torb: \n%.17g\n", torb);

    printf("\ndSha: %.17g %.17g %.17g\n", dSha, ds, dSha/ds);
    printf("dAbe: %.17g %.17g %.17g\n", dAbe, da, dAbe/da);
    printf("DDGR / DDSTG: %.17g\n", -torb / delta);
*/
    torb = -delta;


//----------------------------------------------------------------------------------------


    if (param==-1)  return torb;

    /* Now get partial derivatives */

    dk_dm = derivatives.dk_dm;
    dk_dm2 = derivatives.dk_dm2;
    dgamma_dm = derivatives.dgamma_dm;
    dgamma_dm2 = derivatives.dgamma_dm2;
    dsi_dm = derivatives.dsi_dm;
    dsi_dm2 = derivatives.dsi_dm2;
    ddr_dm = derivatives.ddr_dm;
    ddr_dm2 = derivatives.ddr_dm2;
    ddth_dm = derivatives.ddth_dm;
    ddth_dm2 = derivatives.ddth_dm2;
    dpbdot_dm = derivatives.dpbdot_dm;
    dpbdot_dm2 = derivatives.dpbdot_dm2;

//----------------------------------------------------------------------------------------
//  Here non needed

    an0 = sqrt(m/pow(arr,3));

    csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu;
    ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2;
    cx=sw*cume+sqr1me2*cw*su;
    comega=x*(cw*cume-sqr1me2*sw*su);
    cgamma=su;
    cm2=-2*dlogbr;
    csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace; 

    fact1=(m/(2*arr)) * ((m-m2)*m2/pow(m,2) - 9);
    fact2=(3*m/(2*pow(arr,4))) * (1.0 + fact1);
    fact3=(m/2*arr) * (m2/pow(m,2)-2*(m-m2)*m2/pow(m,3));
    fact4=(1+fact1)*3*m/(2*pow(arr,4)*an0);
    fact5=an0*fact1/arr;

    denumm=(1+fact1)/(2*pow(arr,3)*an0) + an0*(fact1/m+fact3);
    denomm=fact4+fact5;
    darrdm=denumm/denomm;

    dnum = an0*(m-2*m2)/(2*arr*m);
    denom = an0*fact1/arr + fact2/an0;

    darrdm2 = dnum/denom;

    dgmdm2 = ((m+2*m2)/arr - (m2*(m+m2)*darrdm2/pow(arr,2)))*ecc/(an*m);
    cdth=-ecc*ecc*x*cw*su/sqr1me2;
    dthdm2 = -dth*darrdm2/arr - (m+m2)/(arr*m);

    dkdm = k/m - k*darrdm/arr;
    dsidm2 = -(m*x/(arr*m2))*(1.0/m2+darrdm2/arr);
    ck = ae*comega;
    dkdm2 = -k*darrdm2/arr;
    cdr = -ecc*x*sw;
    ddrdm2 = -dr*darrdm2/arr - 2*m2/(arr*m);
    dtdm2 = -2*dlogbr;

    csini = 2*m2*(sw*cume+sqr1me2*cw*su)/brace;

    dsidm=-(m*x/(arr*m2))*(-1.0/m+darrdm/arr);
    dpbdm = pbdot/(m-m2) - pbdot/(3*m);
    cpbdot = -csigma*an*pow(tt0,2)/(2*pb);
    ddrdm = -dr/m - dr*darrdm/arr + 6/arr;

    dpbdm2 = pbdot/m2 - pbdot/(m-m2);

    cm2 = dtdm2+cgamma*dgmdm2+csini*dsidm2+ck*dkdm2+cdr*ddrdm2+cdth*dthdm2+cpbdot*dpbdm2; //DDGR
    fact6=1.0/(arr*m);
    fact7=-(m+m2)/(arr*pow(m,2));
    fact8=-(m+m2)*darrdm/(pow(arr,2)*m);
    dgamdm = (ecc*m2/an)*(fact6+fact7+fact8);

    dthdm=-dth/m - dth*darrdm/arr + (7*m-m2)/(arr*m);
    cm = ck*dkdm+cgamma*dgamdm+cdr*ddrdm+cdth*dthdm+cpbdot*dpbdm+csini*dsidm; //DDGR

//    printf("cm and cm2 values:\nDDGR: %.17g %.17g\n", cm, cm2);


//----------------------------------------------------------------------------------------
    an0 = sqrt(m/pow(arr,3));

    csigma = x*(-sw*su+sqr1me2*cw*cu)/onemecu;
    ce = su*csigma-x*sw-ecc*x*cw*su/sqr1me2;
    cx = sw*cume+sqr1me2*cw*su;
    comega = x*(cw*cume-sqr1me2*sw*su);
    cgamma = su;
    csini = 2*bare_m2*(sw*cume+sqr1me2*cw*su)/brace; // change m2 to bare_m2
    ck = ae*comega;
    cdr = -ecc*x*sw;
    cdth = -ecc*ecc*x*cw*su/sqr1me2;
    cpbdot = -csigma*an*pow(tt0,2)/(2*pb);

    ddt_dm2 = -2*dlogbr * bare_m2/m2;


/*
    printf("ddt_dm2 value in DDSTG vs DDGR: %.17g %.17g\n", ddt_dm2, dtdm2);
    printf("dgamma_dm2 value in DDSTG vs DDGR: %.17g %.17g\n", dgamma_dm2, dgmdm2);
    printf("dsi_dm2 value in DDSTG vs DDGR: %.17g %.17g\n", dsi_dm2, dsidm2);
    printf("dk_dm2 value in DDSTG vs DDGR: %.17g %.17g\n", dk_dm2, dkdm2);
    printf("ddr_dm2 value in DDSTG vs DDGR: %.17g %.17g\n", ddr_dm2, ddrdm2);
    printf("ddth_dm2 value in DDSTG vs DDGR: %.17g %.17g\n", ddth_dm2, dthdm2);
    printf("dpbdot_dm2 value in DDSTG vs DDGR: %.17g %.17g\n", dpbdot_dm2, dpbdm2);
    printf("cm2 value in DDSTG vs DDGR: %.17g %.17g\n\n", cm2, ddt_dm2 + cgamma*dgamma_dm2 + csini*dsi_dm2 + ck*dk_dm2 + cdr*ddr_dm2 + cdth*ddth_dm2 + cpbdot*dpbdot_dm2);
    printf("cm2 value in DDSTG / DDGR: %.17g\n\n", cm2 / (ddt_dm2 + cgamma*dgamma_dm2 + csini*dsi_dm2 + ck*dk_dm2 + cdr*ddr_dm2 + cdth*ddth_dm2 + cpbdot*dpbdot_dm2));

    printf("dgamma_dm value in DDSTG vs DDGR: %.17g %.17g\n", dgamma_dm, dgamdm);
    printf("dsi_dm value in DDSTG vs DDGR: %.17g %.17g\n", dsi_dm, dsidm);
    printf("dk_dm value in DDSTG vs DDGR: %.17g %.17g\n", dk_dm, dkdm);
    printf("ddr_dm value in DDSTG vs DDGR: %.17g %.17g\n", ddr_dm, ddrdm);
    printf("ddth_d value in DDSTG vs DDGR: %.17g %.17g\n", ddth_dm, dthdm);
    printf("dpbdot_dm2 value in DDSTG vs DDGR: %.17g %.17g\n", dpbdot_dm, dpbdm);
    printf("cm value in DDSTG vs DDGR: %.17g %.17g\n\n", cm, ck*dk_dm + cgamma*dgamma_dm + cdr*ddr_dm + cdth*ddth_dm + cpbdot*dpbdot_dm + csini*dsi_dm);
    printf("cm value in DDSTG / DDGR: %.17g\n\n", cm / (ck*dk_dm + cgamma*dgamma_dm + cdr*ddr_dm + cdth*ddth_dm + cpbdot*dpbdot_dm + csini*dsi_dm));
*/



    cm2 = ddt_dm2 + cgamma*dgamma_dm2 + csini*dsi_dm2 + ck*dk_dm2 + cdr*ddr_dm2 + cdth*ddth_dm2 + cpbdot*dpbdot_dm2;

    cm = ck*dk_dm + cgamma*dgamma_dm + cdr*ddr_dm + cdth*ddth_dm + cpbdot*dpbdot_dm + csini*dsi_dm;

//    printf("DDSTG: %.17g %.17g\n", cm, cm2);


    if (param==-2)  /* Set derived parameters */
    {
        /* calculated values (assuming GR) */
        psr[p].param[param_sini].paramSet[0]=1;
        psr[p].param[param_sini].val[0]=si;

        /* Should be xomdot??? */
        psr[p].param[param_omdot].paramSet[0]=1;
        psr[p].param[param_omdot].val[0]=360.0*365.25*xk/(pb/SECDAY);

        psr[p].param[param_gamma].paramSet[0]=1;
        psr[p].param[param_gamma].val[0]=gamma;

        psr[p].param[param_pbdot].paramSet[0]=1;
        psr[p].param[param_pbdot].val[0]=pbdot;

        psr[p].param[param_dtheta].paramSet[0]=1;
        psr[p].param[param_dtheta].val[0]=dth;

        psr[p].param[param_dr].paramSet[0]=1;
        psr[p].param[param_dr].val[0]=dr;

        return 0;

    }

    if (param==param_pb)
        return -csigma*an*SECDAY*tt0/(pb*SECDAY); 
    else if (param==param_a1)
        return cx;
    else if (param==param_ecc)
        return ce;
    else if (param==param_om)
        return comega;
    else if (param==param_t0)
        return -csigma*an*SECDAY;
    else if (param==param_pbdot)
        return 0.5*tt0*(-csigma*an*SECDAY*tt0/(pb*SECDAY));
    else if (param==param_xpbdot)
        return 0.5*tt0*(-csigma*an*SECDAY*tt0/(pb*SECDAY));
    else if (param==param_sini)
        return csi;
    else if (param==param_m2)
        return cm2*SUNMASS;
    else if (param==param_mtot)
        return cm*SUNMASS;
    else if (param==param_a1dot) /* Also known as xdot */
        return cx*tt0;

    return 0.0;
}

void updateDDSTG(pulsar *psr,double val,double err,int pos)
{
    if (pos==param_pb)
    {
        psr->param[param_pb].val[0] += val/SECDAY;
        psr->param[param_pb].err[0]  = err/SECDAY;
    }
    else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_sini || pos==param_m2
            || pos==param_mtot)
    {
        psr->param[pos].val[0] += val;
        psr->param[pos].err[0]  = err;
    }
    else if (pos==param_om)
    {
        psr->param[pos].val[0] += val*180.0/M_PI;
        psr->param[pos].err[0]  = err*180.0/M_PI;
    }
    else if (pos==param_pbdot || pos==param_xpbdot)
    {
        psr->param[pos].val[0] += val;
        psr->param[pos].err[0]  = err;
    }
    else if (pos==param_a1dot)
    {
        psr->param[pos].val[0] += val;
        psr->param[pos].err[0]  = err;
    }
    else if (pos==param_omdot)
    {
        psr->param[pos].val[0] += val*(SECDAY*365.25)*180.0/M_PI;
        psr->param[pos].err[0]  = err*(SECDAY*365.25)*180.0/M_PI;
    }
}

/* Given system masses of m,m2 and Keplerian parameters x,ecc and an, this 
 * routine calculates values of arr, ar, si, gamma and pbdot under GR */

void mass2ddstg(double am,double am2,double x,double ecc,double an,double *arr,double *ar,
        double *xk,double *si,double *gamma,double *pbdot,double *dr,double *dth, DEF_factors *factors, DEF_PPK_derivatives *derivatives)
{
    double SUNMASS = 4.925490947e-6, BARE_SUNMASS;
    double ARRTOL = 1.0e-10;

    double m,m2,m1,arr0,arrold;

    double alpha0, beta0, mA, alphaA, betaA, kA, mB, alphaB, betaB, kB;
    double GAB, XA, XB, nu, gamma_AB, eps;
    double pbdot_m, pbdot_d, pbdot_qphi, pbdot_qg;
    double W_mult;

    alpha0 = factors->alpha0;
    beta0 = factors->betaA;
    mA = factors->mA;
    mB = factors->mB;
    alphaA = factors->alphaA;
    betaA = factors->betaA;
    kA = factors->kA;
    alphaB = factors->alphaB;
    betaB = factors->betaB;
    kB = factors->kB;

    BARE_SUNMASS = SUNMASS / (1 + pow(alpha0,2));

    m=am*SUNMASS;
    m2 = am2*SUNMASS;
    m1 = m-m2;

    XA = m1/m;
    XB = m2/m;
    GAB = (1.0 + alphaA*alphaB) / (1.0 + pow(alpha0,2));
    nu = XA * XB;

    if (m<0)
    {
        printf("ERROR: problem in DDSTG model (mtot < 0)\n");
        exit(1);
    }

    gamma_AB = 1.0 - 2*alphaA*alphaB / (1 + alphaA*alphaB);         //Damour2009 76
    eps = 2*gamma_AB + 1.0;                                          //Damour Tailor1992 near 3.4

    *arr = pow(GAB*m/pow(an,2), 1.0/3.0) * (1 - 1.0/6.0*(5*eps + 3.0 - 2*nu) * pow(GAB*m*an, 2.0/3.0));  //Damour Tailor1992 3.9 correct 

/*
//    printf("arr value in DDSTG vs DDGR: %.17g ", *arr);
    arr0 = pow(m/(an*an),1.0/3.0);
    *arr = arr0;
    do {
        arrold = *arr;
        *arr = arr0*pow(1.0+(m1*m2/pow(m,2) - 9.0)*0.5*m/(*arr),2.0/3.0);
    } while (fabs(((*arr)-arrold)/(*arr)) > ARRTOL);
    *arr = arr0*pow(1.0+(m1*m2/pow(m,2) - 9.0)*0.5*m/(*arr),2.0/3.0); //GR
//    printf("%.17g\n", *arr);
*/

    *ar = (*arr)*m2/m;        //correct


    *si=x/(*ar); //the exact value

    *xk = 3.0/(1-pow(ecc,2))*pow(GAB*m*an, 2.0/3.0) * ((1-1.0/3*alphaA*alphaB)/(1+alphaA*alphaB) - (XA*betaB*pow(alphaA,2) + XB*betaA*pow(alphaB,2)) / (6.0 * pow(1+alphaA*alphaB,2))  ); //Damour1996 5.3

 
//    *xk = 3.0*m/((*arr)*(1-pow(ecc,2))) * ((1-1.0/3*alphaA*alphaB)/(1+alphaA*alphaB) - (XA*betaB*pow(alphaA,2) + XB*betaA*pow(alphaB,2)) / (6.0 * pow(1+alphaA*alphaB,2))  ); // modified

    W_mult = 1.0 / (1 - 1.0/6.0*(5*eps + 3.0 - 2*nu) * pow(GAB*m*an, 2.0/3.0));

//    *xk *= W_mult;

//    printf("k value in DDSTG vs DDGR: %.17g ", *xk);
//    *xk=3.0*m/((*arr)*(1.0-ecc*ecc)); //GR
//    printf("%.17g\n", *xk);

    *gamma = (ecc*pow(an*GAB*m,2.0/3.0)*m2*(m + alphaB*kA*m + m2 + alphaA*alphaB*m2))/((an + alphaA*alphaB*an)*pow(m,2)); //Damour1996 4.9

//    *gamma *= W_mult;

//    printf("gamma value in DDSTG vs DDGR: %.17g ", *gamma);
//    *gamma = ecc*m2*(m1+2*m2)/(an*(*arr)*m); //GR
//    printf("%.17g\n", *gamma);

    pbdot_m = (-3*pow(ecc,2)*(4 + pow(ecc,2))*pow(an*GAB*m,5.0/3.0)*m1*m2*
     pow(alphaA + alphaB*(5.0/3.0 + betaA/(1 + alphaA*alphaB)) + (alphaA*betaB)/(1 + alphaA*alphaB) + 
       (2*(alphaA - alphaB)*m2)/(3.*m),2)*M_PI)/(4.*(1 + alphaA*alphaB)*pow(1 - pow(ecc,2),3.5)*pow(m,2)); //Damour1992 6.52a

    pbdot_d = (2*m1*m2*((pow(alphaA - alphaB,2)*an*(-2 + pow(ecc,2) + pow(ecc,4))*GAB*m)/2. - 
       (2*(alphaA - alphaB)*pow(an*GAB*m,5.0/3.0)*
          (((1 + alphaA*alphaB)*(32 + 124*pow(ecc,2) + 19*pow(ecc,4))*(m - 2*m2)*(alphaA*m1 + alphaB*m2))/4. + 
            5*(1 + 3*pow(ecc,2) + (3*pow(ecc,4))/8.)*m*(alphaA*betaB*m1 - alphaB*betaA*m2)))/
        (5.*(1 + alphaA*alphaB)*pow(m,2)))*M_PI)/((1 + alphaA*alphaB)*pow(1 - pow(ecc,2),3.5)*pow(m,2)); //Damour1992 6.52b

    pbdot_d = -2*M_PI/(1+alphaA*alphaB)*nu*(GAB*m*an) * 
      (1+pow(ecc,2)/2) * pow(1-pow(ecc,2), -2.5)*pow(alphaA - alphaB,2) -
      4*M_PI/(1+alphaA*alphaB)*nu*pow(GAB*m*an, 5.0/3.0) * 
      pow(1-pow(ecc,2), -3.5) * (8.0/5 *(1 + (31.0/8)*pow(ecc,2) + 
       (19.0/32)*pow(ecc,4)) * (alphaA-alphaB) * (alphaA*XA + alphaB*XB) *
      (XA - XB) + (1 + 3*pow(ecc,2) + (3.0/8)*pow(ecc,4)) * (alphaA-alphaB) * 
      (betaB*alphaA*XA - betaA*alphaA*XB) / (1+alphaA*alphaB));

    pbdot_qphi = -((96 + 292*pow(ecc,2) + 37*pow(ecc,4))*pow(an*GAB*m,5.0/3.0)*m1*m2*pow(alphaB*m1 + alphaA*m2,2)*M_PI)/
    (15.*(1 + alphaA*alphaB)*pow(1 - pow(ecc,2),3.5)*pow(m,4)); //Damour1992 6.52c

    pbdot_qg = (-2*(96 + 292*pow(ecc,2) + 37*pow(ecc,4))*pow(an*GAB*m,5.0/3.0)*m1*m2*M_PI)/
    (5.*pow(1 - pow(ecc,2),3.5)*m*(m + alphaA*alphaB*m)); //Damour1992 6.52d

    *pbdot = pbdot_m + pbdot_d + pbdot_qphi + pbdot_qg;

//    printf("pbdot value in DDGR / DDSTG: %.17g \n", -(96.0*2.0*M_PI/5.0)*pow(an,5.0/3.0)*pow(1.0-pow(ecc,2),-3.5) * (1+(73.0/24)*pow(ecc,2) + (37.0/96)*pow(ecc,4)) * m1*m2*pow(m,-1.0/3.0) / (*pbdot));
//    printf("pbdot value in DDSTG vs DDGR: %.17g ", *pbdot);
//    *pbdot = -(96.0*2.0*M_PI/5.0)*pow(an,5.0/3.0)*pow(1.0-pow(ecc,2),-3.5) * (1+(73.0/24)*pow(ecc,2) + (37.0/96)*pow(ecc,4)) * m1*m2*pow(m,-1.0/3.0); //GR
//    printf("%.17g\n", *pbdot);

    *dr = -((pow(an*GAB*m,2.0/3.0)*((-3 + alphaA*alphaB)*pow(m,2) + alphaB*(-alphaA + kA)*m*m2 + 
         (1 + alphaA*alphaB)*pow(m2,2)))/((1 + alphaA*alphaB)*pow(m,2)));

//    *dr *= W_mult;

    *dth = -(pow(an*GAB*m,2.0/3.0)*((-7 + alphaA*alphaB)*pow(m,2) + (1 + alphaA*alphaB)*pow(m2,2) + 
        2*m*(m2 + alphaB*kA*m2)))/(2.*(1 + alphaA*alphaB)*pow(m,2));

//    *dth *= W_mult;

    derivatives->dk_dm2 = 0.0;
    derivatives->dgamma_dm2 = (*gamma) * (1.0/m2 + 1.0/(m+m2));
    derivatives->dsi_dm2 = -*si/m2;
    derivatives->ddr_dm2 = *dr * 2*m2/(-3*pow(m,2) + pow(m2,2));
    derivatives->ddth_dm2 = *dth * 2*(m+m2)/(-7*pow(m,2) + 2*m*m2 + pow(m2,2));
    derivatives->dpbdot_dm2 = *pbdot * (1.0/m2 + 1.0/(-m+m2));

    derivatives->dk_dm = *xk * 2.0/(3*m);
    derivatives->dgamma_dm = *gamma * (-4.0/(3*m) + 1.0/(m+m2));
    derivatives->ddr_dm = *dr * (-4.0/(3*m) - 6*m/(-3*pow(m,2) + pow(m2,2)));
    derivatives->ddth_dm = *dth * (-4.0/(3*m) - 2*(7*m-m2)/(-7*pow(m,2) + 2*m*m2 + pow(m2,2)));
    derivatives->dpbdot_dm = *pbdot * (2*m+m2)/(3*pow(m,2) - 3*m*m2);
    derivatives->dsi_dm = *si * 2.0/(3*m);

}

void DEF_tables_load(DEF_tables *tables)
{
    FILE *file;
    char fileName[256], tempstr[256];
    int N_mA, N_alpha0, N_beta0;

    printf("eosName = %s\n", tables->eosName);

    sprintf(fileName,"%s/data_ddstg/%s/dims.dat", getenv(TEMPO2_ENVIRON), tables->eosName);

    printf("Reading from %s\n", fileName);
    file = fopen(fileName, "r");
    fgets(tempstr, 256, file);
    fscanf(file, "%i %i %i", &(tables->N_mA), &(tables->N_alpha0), &(tables->N_beta0));
    fclose(file);
    printf("N_mA, N_alpha0, N_beta0 are %i %i %i\n", tables->N_mA, tables->N_alpha0, tables->N_beta0);

    N_mA = tables->N_mA;
    N_alpha0 = tables->N_alpha0; 
    N_beta0 = tables->N_beta0;

    tables->alpha0_table = (double *)malloc(N_alpha0 * sizeof(double));
    tables->beta0_table = (double *)malloc(N_beta0 * sizeof(double));

    sprintf(fileName,"%s/data_ddstg/%s/alpha0.dat", getenv(TEMPO2_ENVIRON), tables->eosName);
    printf("Reading from %s\n", fileName);
    file = fopen(fileName, "r");
    for (int i = 0; i < N_alpha0; i++) {
        fscanf(file, "%lf,", &(tables->alpha0_table[i]) );
    }
    fclose(file);

    sprintf(fileName,"%s/data_ddstg/%s/beta0.dat", getenv(TEMPO2_ENVIRON), tables->eosName);
    printf("Reading from %s\n", fileName);
    file = fopen(fileName, "r");
    for (int i = 0; i < N_beta0; i++) {
        fscanf(file, "%lf,", &(tables->beta0_table[i]) );
    }
    fclose(file);

    printf("%lf\n", tables->beta0_table[11]);

    tables->mA_table = allocate_3Darray(N_alpha0, N_beta0, N_mA);
    tables->alphaA_table = allocate_3Darray(N_alpha0, N_beta0, N_mA);
    tables->betaA_table = allocate_3Darray(N_alpha0, N_beta0, N_mA);
    tables->kA_table = allocate_3Darray(N_alpha0, N_beta0, N_mA);

//  free_3Darray(tables->mA_table, N_alpha0, N_beta0, N_mA);
//  free_3Darray(tables->alphaA_table, N_alpha0, N_beta0, N_mA);
//  free_3Darray(tables->betaA_table, N_alpha0, N_beta0, N_mA);
//  free_3Darray(tables->kA_table, N_alpha0, N_beta0, N_mA);

    sprintf(fileName,"%s/data_ddstg/%s/massA.dat", getenv(TEMPO2_ENVIRON), tables->eosName);
    read_3Darray(fileName, tables->mA_table, N_alpha0, N_beta0, N_mA);

    sprintf(fileName,"%s/data_ddstg/%s/alphaA.dat", getenv(TEMPO2_ENVIRON), tables->eosName);
    read_3Darray(fileName, tables->alphaA_table, N_alpha0, N_beta0, N_mA);

    sprintf(fileName,"%s/data_ddstg/%s/betaA.dat", getenv(TEMPO2_ENVIRON), tables->eosName);
    read_3Darray(fileName, tables->betaA_table, N_alpha0, N_beta0, N_mA);

    sprintf(fileName,"%s/data_ddstg/%s/kA.dat", getenv(TEMPO2_ENVIRON), tables->eosName);
    read_3Darray(fileName, tables->kA_table, N_alpha0, N_beta0, N_mA);

//  for (int i = 0; i < N_mA; i++) {
//      printf("%lf ", tables->alphaA_table[1][1][i]);
//  }
}

void DEF_tables_free(DEF_tables *tables)
{
    free(tables->alpha0_table);
    free(tables->beta0_table);

    free_3Darray(tables->mA_table, tables->N_alpha0, tables->N_beta0, tables->N_mA);
    free_3Darray(tables->alphaA_table, tables->N_alpha0, tables->N_beta0, tables->N_mA);
    free_3Darray(tables->betaA_table, tables->N_alpha0, tables->N_beta0, tables->N_mA);
    free_3Darray(tables->kA_table, tables->N_alpha0, tables->N_beta0, tables->N_mA);

    puts("The memory for DEF tables is freed");
}


void allocate_DEF_arrays(pulsar *psr, DEF_tables *tables)
{
    psr->N_mA = tables->N_mA;
    psr->mA_array = (double *)malloc(psr->N_mA * sizeof(double));
    psr->alphaA_array = (double *)malloc(psr->N_mA * sizeof(double));
    psr->betaA_array = (double *)malloc(psr->N_mA * sizeof(double));
    psr->kA_array = (double *)malloc(psr->N_mA * sizeof(double));
}

void interpolate_DEF_arrays(pulsar *psr, DEF_tables *tables)
{
    double alpha0 = psr->param[param_alpha0].val[0], beta0 = psr->param[param_beta0].val[0];
    int i_alpha0, i_beta0;
    int N_mA = tables->N_mA, N_alpha0 = tables->N_alpha0, N_beta0 = tables->N_beta0;
    double x_alpha0, x_beta0, x_logalpha0;

    double *mA_array, *alphaA_array, *betaA_array, *kA_array;
    double *alpha0_table, *beta0_table;
    double ***mA_table, ***alphaA_table, ***betaA_table, ***kA_table;

    alpha0_table = tables->alpha0_table;
    beta0_table = tables->beta0_table;

    mA_table = tables->mA_table;
    alphaA_table = tables->alphaA_table;
    betaA_table = tables->betaA_table;
    kA_table = tables->kA_table;

    mA_array = psr->mA_array;
    alphaA_array = psr->alphaA_array;
    betaA_array = psr->betaA_array;
    kA_array = psr->kA_array;

    i_alpha0 = locate_in_array(alpha0, alpha0_table, N_alpha0);
    i_beta0 = locate_in_array(beta0, beta0_table, N_beta0);

    printf("i_alpha0, i_beta0: %i %i\n", i_alpha0, i_beta0);
    printf("alpha0: %lf < %lf < %lf\n", alpha0_table[i_alpha0], alpha0, alpha0_table[i_alpha0+1]);
    printf("beta0: %lf < %lf < %lf\n",  beta0_table[i_beta0], beta0, beta0_table[i_beta0+1]);

    x_beta0 = (beta0 - beta0_table[i_beta0]) / (beta0_table[i_beta0+1] - beta0_table[i_beta0]);
    x_alpha0 = (alpha0 - alpha0_table[i_alpha0]) / (alpha0_table[i_alpha0+1] - alpha0_table[i_alpha0]);
    x_logalpha0 = log(alpha0 / alpha0_table[i_alpha0]) / log(alpha0_table[i_alpha0+1] / alpha0_table[i_alpha0]);

    printf("x_alpha0, x_beta0, x_logalpha0: %lf %lf %lf\n", x_alpha0, x_beta0, x_logalpha0);

    for (int i = 0; i < N_mA; i++) {
        mA_array[i] = bilinear_interpolation(mA_table[i_alpha0][i_beta0][i], mA_table[i_alpha0+1][i_beta0][i], mA_table[i_alpha0][i_beta0+1][i], mA_table[i_alpha0+1][i_beta0+1][i], x_alpha0, x_beta0);
        alphaA_array[i] = bilinear_interpolation(alphaA_table[i_alpha0][i_beta0][i], alphaA_table[i_alpha0+1][i_beta0][i], alphaA_table[i_alpha0][i_beta0+1][i], alphaA_table[i_alpha0+1][i_beta0+1][i], x_alpha0, x_beta0);
        betaA_array[i] = bilinear_interpolation(betaA_table[i_alpha0][i_beta0][i], betaA_table[i_alpha0+1][i_beta0][i], betaA_table[i_alpha0][i_beta0+1][i], betaA_table[i_alpha0+1][i_beta0+1][i], x_alpha0, x_beta0);
        kA_array[i] = bilinear_interpolation(kA_table[i_alpha0][i_beta0][i], kA_table[i_alpha0+1][i_beta0][i], kA_table[i_alpha0][i_beta0+1][i], kA_table[i_alpha0+1][i_beta0+1][i], x_alpha0, x_beta0);
    }

    printf("\nalphaA interpolation box:\n");
    printf("%lf %lf %lf\n", alphaA_table[i_alpha0][i_beta0+1][50], x_alpha0, alphaA_table[i_alpha0+1][i_beta0+1][50]);
    printf("%lf %lf %lf\n", x_beta0, alphaA_array[50], x_beta0);
    printf("%lf %lf %lf\n\n", alphaA_table[i_alpha0][i_beta0][50], x_alpha0, alphaA_table[i_alpha0+1][i_beta0][50]);
}

void interpolate_DEF_params(pulsar *psr, double mA, double *alphaA, double *betaA, double *kA)
{
    int i_mA, N_mA = psr->N_mA;
    double *mA_array, *alphaA_array, *betaA_array, *kA_array;
    mA_array = psr->mA_array;
    alphaA_array = psr->alphaA_array;
    betaA_array = psr->betaA_array;
    kA_array = psr->kA_array;

    i_mA = locate_in_array(mA, mA_array, N_mA);

    *alphaA = interpolate_linear(mA, i_mA, mA_array, alphaA_array);
    *betaA = interpolate_linear(mA, i_mA, mA_array, betaA_array);
    *kA = interpolate_linear(mA, i_mA, mA_array, kA_array);

    printf("i_mA = %i\n", i_mA);
    printf("mA    : %lf < %lf < %lf\n", mA_array[i_mA], mA, mA_array[i_mA+1]);
    printf("alphaA: %lf < %lf < %lf\n", alphaA_array[i_mA], *alphaA, alphaA_array[i_mA+1]);
    printf("betaA : %lf < %lf < %lf\n", betaA_array[i_mA], *betaA, betaA_array[i_mA+1]);
    printf("kA    : %lf < %lf < %lf\n\n", kA_array[i_mA], *kA, kA_array[i_mA+1]);

}

void interpolate_DEF_factors(pulsar *psr, DEF_factors *factors)
{
    int i_mA, N_mA = psr->N_mA;
    double *mA_array, *alphaA_array, *betaA_array, *kA_array;
//    double alpha0, beta0, mA, alphaA, betaA, kA, mB, alphaB, betaB, kB;

    mA_array = psr->mA_array;
    alphaA_array = psr->alphaA_array;
    betaA_array = psr->betaA_array;
    kA_array = psr->kA_array;

    if ((factors->alpha0 == 0.0) && (factors->beta0 == 0.0))
    {
        factors->alphaA = 0.0;
        factors->betaA = 0.0;
        factors->kA = 0.0;
        factors->alphaB = 0.0;
        factors->betaB = 0.0;
        factors->kB = 0.0;
    }
    else
    {
        i_mA = locate_in_array(factors->mA, mA_array, N_mA);
    
        factors->alphaA = interpolate_linear(factors->mA, i_mA, mA_array, alphaA_array);
        factors->betaA = interpolate_linear(factors->mA, i_mA, mA_array, betaA_array);
        factors->kA = interpolate_linear(factors->mA, i_mA, mA_array, kA_array);
 
/*    
        printf("i_mA = %i\n", i_mA);
        printf("mA    : %lf < %lf < %lf\n", mA_array[i_mA], factors->mA, mA_array[i_mA+1]);
        printf("alphaA: %lf < %lf < %lf\n", alphaA_array[i_mA], factors->alphaA, alphaA_array[i_mA+1]);
        printf("betaA : %lf < %lf < %lf\n", betaA_array[i_mA], factors->betaA, betaA_array[i_mA+1]);
        printf("kA    : %lf < %lf < %lf\n\n", kA_array[i_mA], factors->kA, kA_array[i_mA+1]);
*/
    
        if (strcmp(factors->companion_type, "WD") == 0)
        {
            factors->alphaB = factors->alpha0;
            factors->betaB = factors->beta0;
            factors->kB = 0.0;
        }
        else if(strcmp(factors->companion_type, "NS") == 0)
        {
            i_mA = locate_in_array(factors->mA, mA_array, N_mA);
            factors->alphaB = interpolate_linear(factors->mB, i_mA, mA_array, alphaA_array);
            factors->betaB = interpolate_linear(factors->mB, i_mA, mA_array, betaA_array);
            factors->kB = interpolate_linear(factors->mB, i_mA, mA_array, kA_array);
        }
        else if(strcmp(factors->companion_type, "BH") == 0)
        {
            factors->alphaB = 0.0;
            factors->betaB = 0.0;
            factors->kB = 0.0;
        }
    }
}


double interpolate_linear(double x_value, int ind, double *x_array, double *y_array)
{
    return (y_array[ind]*(x_array[ind+1] - x_value) + y_array[ind+1]*(x_value - x_array[ind])) / (x_array[ind+1] - x_array[ind]);
}

double bilinear_interpolation(double f00, double f10, double f01, double f11, double x, double y)
{
    return f00*(1-x)*(1-y) + f10*x*(1-y) + f01*(1-x)*y + f11*x*y;
}

int locate_in_array(double value, double *array, int N)
{
    int i, i_left, i_right, i_middle;

    i_left = 0;
    i_right = N-1;

    do
    {
        i_middle = (i_right + i_left) / 2;
        if ((array[N-1] >= array[0]) == (value >= array[i_middle]))
            i_left = i_middle;
        else
            i_right = i_middle;
    } while(i_right - i_left > 1);

    if (value == array[0])
        i = 0;
    else if (value >= array[N-1])
        i = N-2;
    else
        i = i_left;

//    printf("debug %f %f \n", value, array[N-1]);

    return i;
}

void read_3Darray(char *fileName, double ***array, int dim1, int dim2, int dim3)
{
    FILE *file;
    double dummy1, dummy2;
    printf("Reading from %s\n", fileName);
    file = fopen(fileName, "r");
    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
                fscanf(file, "%lf%lf", &dummy1, &dummy2);
            for (int k = 0; k < dim3; k++) {
                if(fscanf(file, "%lf,", &(array[i][j][k])) != 1)
                {
                    fprintf(stderr, "Error reading table %s\n", fileName);
                    exit(0);
                }
            }
        }
    }
    fclose(file);
}


double*** allocate_3Darray(int dim1, int dim2, int dim3)
{
    double ***array;
    array = (double ***)malloc(dim1*sizeof(double**));
    if (array == NULL)
    {
        fprintf(stderr, "Out of memory");
        exit(0);
    }
    for (int i = 0; i< dim1; i++) {
        array[i] = (double **) malloc(dim2*sizeof(double *));
        if (array[i] == NULL)
        {
            fprintf(stderr, "Out of memory");
            exit(0);
        }
        for (int j = 0; j < dim2; j++) {
            array[i][j] = (double *)malloc(dim3*sizeof(double));
            if (array[i][j] == NULL)
            {
                fprintf(stderr, "Out of memory");
                exit(0);
            }
        }
    }
    return array;
}

void free_3Darray(double ***array, int dim1, int dim2, int dim3)
{
    for (int i = 0; i< dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}