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

/* This plugin allows the user to simulate an isotropic gravitational wave background comprised of arbitrary polarization states. */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "T2toolkit.h"
#include "GWsim.h"
#include "TKfit.h"
#include <cpgplot.h>

using namespace std;

void doPlot(pulsar *psr,int npsr,gwgeneralSrc *gw,longdouble **gwRes,longdouble timeOffset,int *numsources,longdouble tspan,gwgenSpec gwAmps);
longdouble getTspan(pulsar *psr,int npsr);
void plotResiduals(pulsar *psr,longdouble **gwRes,int p,longdouble timeOffset,int plotType);
void plotSpectrum(gwgeneralSrc *gw,int *numsources,longdouble tspan,gwgenSpec gwAmps);
void plotPosn(pulsar *psr,int npsr,gwgeneralSrc *gw,int *numsources);
void draw_grid(double start_gl,double end_gl,double start_gb,double end_gb,double gstep,double bstep,int celestialCoords);
void convertXY_celestial(double raj,double decj,double *retx,double *rety);

void help() /* Display help */
{
  printf("--------------------------------------\n");
  printf("Command line inputs:\n\n");
  printf("-addWhite add white noise based on the TOA error bar size\n");
  printf("-alphaTT GW transverse tensor spectral exponent (usually -0.666) \n");
  printf("-alphaST GW transverse scalar spectral exponent (usually -0.666) \n");
  printf("-alphaSL GW longitudinal scalar spectral exponent (usually -0.666) \n");
  printf("-alphaVL GW longitudinal vector spectral exponent (usually -0.666) \n");
  printf("-clock simulate clock errors instead of a GWB\n");
  printf("-dist  pulsar distance in kpc\n");
  printf("-f     parfile.par timfile.tim:  input par and tim files\n");
  printf("-flo   lowest GW frequency to simulate (default = 0.01/Tspan) (Hz) \n");
  printf("-fhi   highest GW frequency to simulate (default = d^-1) (Hz)\n");
  printf("-gwampTT GW transverse tensor amplitude in dimensionless units\n");
  printf("-gwampST GW transverse scalar amplitude in dimensionless units\n");
  printf("-gwampSL GW longitudinal scalar amplitude in dimensionless units\n");
  printf("-gwampVL GW longitudinal vector amplitude in dimensionless units\n");
  printf("-h     this help file\n");
  printf("-linear Use linear spacing for GW frequencies (default = logarithmic spacing)\n");
  printf("-ngw   Number of gravitational waves to simulate for each polarisation state\n");
  printf("-plot  Plot the GW positions and induced residuals\n");
  printf("-seed  Set random number seed (choose negative integer)\n");
  printf("-writebkgrd bkgrdfile.dat: write generated background out to file bkgrdfile.dat\n");
  printf("-writebkgrdid I: if specified, uses integer I to identify the background realisation in the output file. Default is 0 if unspecified.\n");
  printf("-zero  Set all the post-fit residuals to zero before continuing\n");
  printf("\n\n");
  printf("Example usage: tempo2 -gr GWgeneralbkgrd -f tt.par 0437.2048.tim -dist 1 -gwampTT 1e-15 -alphaTT -0.666666 -gwampST 1e-15 -alphaST -0.666666 -gwampSL 1e-15 -alphaSL -0.666666 -gwampVL 1e-15 -alphaVL -0.666666 -ngw 1000 -plot\n\n");
  printf("Plot options\n\n");
  printf("1 Plot sky plot with pulsar position indicated and all GW source positions\n");
  printf("2 Plot GW power versus GW frequency\n");
  printf("3 Plot the induced residuals caused by the GW background\n");
  printf("4 Plot the induced residuals caused by the GW background after quadratic removal\n");
  printf("5 Plot the pre-fit timing residuals\n");
  printf("6 Plot the post-fit timing residuals\n");
  printf("h This help file\n");
  printf("q Quit\n");
  printf("s Save residuals and site arrival times\n");
  printf("p Select which pulsar number to plot\n");
  printf("--------------------------------------\n");
  
}


/* The main function called from the TEMPO2 package is 'graphicalInterface' */
/* Therefore this function is required in all plugins                       */
extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
  char parFile[MAX_PSR][MAX_FILELEN];
  char timFile[MAX_PSR][MAX_FILELEN];
  FILE *fout,*fout2;
  int i,j,k,p;
  double globalParameter;
  int writebkgrd=0;
  char bkgrdFile[MAX_FILELEN];
  int bkgrdreal=0;
  int setgwAmpTT=0,setAlphaTT=0;
  int setgwAmpST=0,setAlphaST=0;
  int setgwAmpSL=0,setAlphaSL=0;
  int setgwAmpVL=0,setAlphaVL=0;
  int plotIt=0;
  gwgenSpec gwAmps;
  longdouble timeOffset;
  longdouble ra_p,dec_p;
  double flo=0.0,fhi=0.0;
  longdouble kp[3];            /* Vector pointing to pulsar           */
  longdouble tspan;
  longdouble time;
  longdouble **gwRes;
  longdouble dist[MAX_PSR];
  longdouble mean;
  int clock=0;
  int distNum=0;
  int logspacing=1;
  int ngw=0;
  int addWhite=0;
  long seed=TKsetSeed();
  int zeroResiduals=0;
  double scale;
  char fname[100];
  int outGW=0;
  int identicalTimes=0;
  gwgeneralSrc *gw;
  int numsources[4];

  *npsr = 0;  /* For a graphical interface that only shows results for one pulsar */

  printf("Graphical Interface: GWgeneralbkgrd\n");
  printf("Author:              J. Gair, adapted from GWbkgrd by G. Hobbs\n");
  printf("Version:             1.0\n");
  printf(" --- type 'h' for help information\n");


  /* Obtain the .par and the .tim file from the command line */
  if (argc==4) /* Only provided .tim name */
    {
      strcpy(timFile[0],argv[3]);
      strcpy(parFile[0],argv[3]);
      parFile[0][strlen(parFile[0])-3] = '\0';
      strcat(parFile[0],"par");
    }

  gwAmps.tensor_amp=gwAmps.st_amp=gwAmps.sl_amp=gwAmps.vl_amp=0.;

  /* Obtain all parameters from the command line */
  for (i=2;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)  // Read parameter file and arrival time files
	{
 	  strcpy(parFile[*npsr],argv[i+1]); 
	  strcpy(timFile[*npsr],argv[i+2]);	  
	  (*npsr)++;
	}
      else if (strcmp(argv[i],"-identicalTimes")==0)
	identicalTimes=1;
      else if (strcmp(argv[i],"-outGW")==0)
	outGW=1;
      else if (strcmp(argv[i],"-clock")==0)
	clock=1;
      else if (strcmp(argv[i],"-addWhite")==0)
	{
	  printf("Using TOA error bars to add white noise\n");
	  addWhite=1;
	}
      else if (strcmp(argv[i],"-h")==0)
	{
	  help();
	  exit(1);
	}
      else if (strcmp(argv[i],"-dist")==0) // Distance in kpc
	{
      dist[distNum] = parse_longdouble(argv[++i]);
	  dist[distNum]*=longdouble(3.086e19);
	  distNum++;
	}
      else if (strcmp(argv[i],"-gwampTT")==0)
	{sscanf(argv[++i],"%lf",&(gwAmps.tensor_amp)); setgwAmpTT=1;}
      else if (strcmp(argv[i],"-gwampST")==0)
	{sscanf(argv[++i],"%lf",&(gwAmps.st_amp)); setgwAmpST=1;}
      else if (strcmp(argv[i],"-gwampSL")==0)
	{sscanf(argv[++i],"%lf",&(gwAmps.sl_amp)); setgwAmpSL=1;}
      else if (strcmp(argv[i],"-gwampVL")==0)
	{sscanf(argv[++i],"%lf",&(gwAmps.vl_amp)); setgwAmpVL=1;}
      else if (strcmp(argv[i],"-alphaTT")==0)
	{sscanf(argv[++i],"%lf",&(gwAmps.tensor_alpha)); setAlphaTT=1;}
      else if (strcmp(argv[i],"-alphaST")==0)
	{sscanf(argv[++i],"%lf",&(gwAmps.st_alpha)); setAlphaST=1;}
      else if (strcmp(argv[i],"-alphaSL")==0)
	{sscanf(argv[++i],"%lf",&(gwAmps.sl_alpha)); setAlphaSL=1;}
      else if (strcmp(argv[i],"-alphaVL")==0)
	{sscanf(argv[++i],"%lf",&(gwAmps.vl_alpha)); setAlphaVL=1;}
      else if (strcmp(argv[i],"-ngw")==0)
	{sscanf(argv[++i],"%d",&ngw);}
      else if (strcmp(argv[i],"-plot")==0)
	plotIt=1;
      else if (strcmp(argv[i],"-flo")==0)
	sscanf(argv[++i],"%lf",&flo);
      else if (strcmp(argv[i],"-writebkgrd")==0) {
        writebkgrd=1;
        strcpy(bkgrdFile,argv[i+1]);
      } else if (strcmp(argv[i],"-writebkgrdid")==0)
        sscanf(argv[++i],"%i",&bkgrdreal);
      else if (strcmp(argv[i],"-zero")==0)
	zeroResiduals=1;
      else if (strcmp(argv[i],"-fhi")==0)
	sscanf(argv[++i],"%lf",&fhi);
      else if (strcmp(argv[i],"-seed")==0)
	sscanf(argv[++i],"%ld",&seed);
      else if (strcmp(argv[i],"-linear")==0)
	logspacing=0;
    }
  

  // Check that all the parameters are set
  if ((setgwAmpTT+setgwAmpST+setgwAmpSL+setgwAmpVL)==0)
    {
      gwAmps.tensor_amp = 1e-20;
      printf("WARNING: no GW background amplitudes specified. Automatically setting the GW transverse tensor amplitude to be 1e-20\nUse -gwampTT option to set manually.\n");
    }

  if ((setAlphaTT==0) && (setgwAmpTT==1))
    {
      gwAmps.tensor_alpha=-2.0/3.0;
      printf("WARNING: automatically setting the characteristic strain exponent for the transverse tensor background to -2/3\nUse -alphaTT option to set manually\n");
    }
  if ((setAlphaST==0) && (setgwAmpST==1))
    {
      gwAmps.st_alpha=-2.0/3.0;
      printf("WARNING: automatically setting the characteristic strain exponent for the transverse scalar background to -2/3\nUse -alphaST option to set manually\n");
    }
  if ((setAlphaSL==0) && (setgwAmpSL==1))
    {
      gwAmps.sl_alpha=-2.0/3.0;
      printf("WARNING: automatically setting the characteristic strain exponent for the longitudinal scalar background to -2/3\nUse -alphaSL option to set manually\n");
    }
  if ((setAlphaVL==0) && (setgwAmpVL==1))
    {
      gwAmps.vl_alpha=-2.0/3.0;
      printf("WARNING: automatically setting the characteristic strain exponent for the longitudinal vector background to -2/3\nUse -alphaVL option to set manually\n");
    }

  scale = pow(86400.0*365.25,gwAmps.tensor_alpha);
  gwAmps.tensor_amp *= scale;
  scale = pow(86400.0*365.25,gwAmps.st_alpha);
  gwAmps.st_amp *= scale;
  scale = pow(86400.0*365.25,gwAmps.sl_alpha);
  gwAmps.sl_amp *= scale;
  scale = pow(86400.0*365.25,gwAmps.vl_alpha);
  gwAmps.vl_amp *= scale;


  if (ngw==0)
    {
      printf("WARNING: automatically setting the number of gravitational waves to 1000\nUse -ngw option to set manually\n");
      ngw=1000;
    }

  for (int i=0;i<4;i++)
	numsources[i]=0;
  if (gwAmps.tensor_amp)
	numsources[0]=ngw;
  if (gwAmps.st_amp)
	numsources[1]=ngw;
  if (gwAmps.sl_amp)
	numsources[2]=ngw;
  if (gwAmps.vl_amp)
	numsources[3]=ngw;

  if (*npsr==0)
    {
      printf("ERROR: Must use -f option to provide a pulsar timing model\n");
      exit(1);
    }
  if (distNum!=*npsr)
    {
      printf("ERROR: Distances not provided for all the pulsars: Npsr = %d, Ndist = %d\nUse -dist to provide distances (in kpc) on the command line.\n",*npsr,distNum);
      exit(1);
    }

  int ngwtot=0;
  for (int kk=0;kk<4;kk++)
	ngwtot+=numsources[kk];
  if ((gw = (gwgeneralSrc *)malloc(sizeof(gwgeneralSrc)*ngwtot))==NULL)
    {
      printf("Unable to allocate memory for %d GW sources\n",ngwtot);
      exit(1);
    }
  gwRes = (longdouble **)malloc(MAX_PSR*sizeof(longdouble*));
  for (i=0;i<MAX_PSR;i++)
    gwRes[i] = (longdouble *)malloc(MAX_OBSN*sizeof(longdouble));

  readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
  readTimfile(psr,timFile,*npsr);         /* Load the arrival times    */
  preProcess(psr,*npsr,argc,argv);
  formBatsAll(psr,*npsr);                 /* Form the barycentric arrival times */
  formResiduals(psr,*npsr,1);           /* Form the residuals                 */

  if (zeroResiduals==1)
    {
      for (j=0;j<5;j++)
	{
	  for (p=0;p<*npsr;p++)
	    {
	      for (i=0;i<psr[p].nobs;i++)
		psr[p].obsn[i].sat -= psr[p].obsn[i].residual/SECDAY;
	    }
	  formBatsAll(psr,*npsr);         /* Form the barycentric arrival times */
	  formResiduals(psr,*npsr,1);   /* Form the residuals                 */
	}
    }


  // Set range of frequencies for GW simulation
  tspan=getTspan(psr,*npsr)*SECDAY;
  if (flo==0)
    {
      flo=0.01/tspan;
      ld_printf("flo = %.5lg, tspan = %.5Lg\n",flo,tspan);
    }
  if (fhi==0)
    {
      fhi = 1.0/(longdouble)SECDAY;
      printf("fhi = %.5lg\n",fhi);
    }

  timeOffset = psr[0].param[param_pepoch].val[0];
   GWgeneralbackground(gw,numsources,&seed,flo,fhi,gwAmps,logspacing);

  if (writebkgrd) {
	fout=fopen(bkgrdFile,"wb");
        GWgeneralbackground_write(gw,fout,ngwtot,bkgrdreal);
        fclose(fout);
  }

  i=0;
  for (p=0;p<4;p++) {
  	for (int kk=0;kk<numsources[p];kk++)
    	{
      	setupgeneralGW(&gw[i]);
    	i++;
    	}
  }
  if (clock==1)
    fout = fopen("signal.dat","w");
  for (p=0;p<*npsr;p++)
    {
      if (outGW==1)
	{
	  char str[128];
	  sprintf(str,"%s.gw",psr[p].name);
	  fout2 = fopen(str,"w");
	}

      if (clock==0)
	{
	  ra_p   = psr[p].param[param_raj].val[0];
	  dec_p  = psr[p].param[param_decj].val[0];
	}
      else
	{
	  ra_p   = psr[0].param[param_raj].val[0];
	  dec_p  = psr[0].param[param_decj].val[0];
	}
      setupPulsar_GWsim(ra_p,dec_p,kp);
      mean=0.0;
      for (i=0;i<psr[p].nobs;i++) 
	{
	  if (identicalTimes==1)
	    time = (psr[0].obsn[i].sat - timeOffset)*SECDAY;
	  else
	    time = (psr[p].obsn[i].sat - timeOffset)*SECDAY;
	  gwRes[p][i] = 0.0;
	  for (k=0;k<ngw;k++)
	    gwRes[p][i]+=calculateResidualgeneralGW(kp,&gw[k],time,dist[p]);
	  mean += gwRes[p][i];
	}
      mean /= (longdouble)psr[p].nobs;
      for (i=0;i<psr[p].nobs;i++)
	{
	  if (clock==1)
	    fprintf(fout,"%.10f %.10g\n",(double)psr[p].obsn[i].sat,(double)gwRes[p][i]);
	  if (outGW==1)
	    {
	      if (identicalTimes==1)
		fprintf(fout2,"%.10f %.10g\n",(double)psr[0].obsn[i].sat,(double)gwRes[p][i]);
	      else
		fprintf(fout2,"%.10f %.10g\n",(double)psr[p].obsn[i].sat,(double)gwRes[p][i]);
	    }
	  gwRes[p][i]-=mean;
	  psr[p].obsn[i].sat += (gwRes[p][i]/(longdouble)SECDAY);
	  if (addWhite==1)
	    psr[p].obsn[i].sat += (longdouble)TKgaussDev(&seed)*(psr[p].obsn[i].toaErr*1.0e-6)/(longdouble)SECDAY;
	}
    }
  if (clock==1)
    fclose(fout);
  if (outGW==2)
    fclose(fout2);
  for (p=0;p<*npsr;p++)
    {
      sprintf(fname,"%s.gwsim.tim",psr[p].name);
      writeTim(fname,psr+p,"tempo2");
    }
  if (plotIt==1)
    {
      formBatsAll(psr,*npsr);                 /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);           /* Form the residuals                 */
      doFitAll(psr,*npsr,0);
      formBatsAll(psr,*npsr);                 /* Form the barycentric arrival times */
      formResiduals(psr,*npsr,1);           /* Form the residuals                 */
      
      doPlot(psr,*npsr,gw,gwRes,timeOffset,numsources,tspan,gwAmps);
    }
  return 0;
}

void doPlot(pulsar *psr,int npsr,gwgeneralSrc *gw,longdouble **gwRes,longdouble timeOffset,int *numsources,longdouble tspan, gwgenSpec gwAmps)
{
  int plot=1;
  int pulsar=0;
  int endit=0;
  int i,j;
  float mx,my;
  char key;
  char fname[100];
  FILE *fout;

  cpgbeg(0,"/xs",1,1);
  cpgpap(0,0.5);
  cpgask(0);
  do {
  if (plot==3 || plot==4 || plot==5 || plot==6)
    plotResiduals(psr,gwRes,pulsar,timeOffset,plot-2);
  else if (plot==2)
    plotSpectrum(gw,numsources,tspan,gwAmps);
  else if (plot==1)
    plotPosn(psr,npsr,gw,numsources);

  cpgcurs(&mx,&my,&key);
  if (key=='q') endit=1;
  else if (key=='1') plot=1;
  else if (key=='2') plot=2;
  else if (key=='5') plot=5;
  else if (key=='6') plot=6;
  else if (key=='3')
    plot=3;
  else if (key=='4')
    plot=4;
  else if (key=='h')
    help();
  else if (key=='s') // Save residuals and site arrival times
    {
      for (j=0;j<npsr;j++)
	{
	  printf("Writing new site arrival times\n");
	  sprintf(fname,"%s.gwsim.tim",psr[j].name);
	  writeTim(fname,psr+j,"tempo2");
	  printf("Writing new residual file\n");
	  sprintf(fname,"%s.res",psr[j].name);
	  fout = fopen(fname,"w");
	  for (i=0;i<psr[j].nobs;i++)
	    {
	      ld_fprintf(fout,"%.15Lf %g %g\n",psr[j].obsn[i].sat,(double)psr[j].obsn[i].residual,(double)psr[j].obsn[i].toaErr);
	    }
	  fclose(fout);
	}
    }
  else if (key=='p') // Select pulsar
    {
      printf("Enter pulsar number (from 0 to %d)\n",npsr-1);
      scanf("%d",&pulsar);
    }

  } while (endit==0);
  cpgend();
}

void plotPosn(pulsar *psr,int npsr,gwgeneralSrc *gw,int *numsources)
{
  int i;
  int ngw=0;
  for (int j=0;j<4;j++)
	ngw+=numsources[j];
  double px[ngw],py[ngw];
  double rad2deg = 180.0/M_PI;
  float fx[ngw],fy[ngw];

  draw_grid(-180,180,-90,90,60,30,1);
  cpglab("","","Gravitational wave source and pulsar positions");
  // Plot the GW source positions
  i=0;
  for (int j=0;j<4;j++) {
	cpgsci(2+j);
  	for (int k=0;k<numsources[j];k++) {
      		//convertXY_celestial((double)(gw[i].phi_g*rad2deg)-180,
		//	  (double)(gw[i].theta_g*rad2deg)-90,&px[0],&py[0]);
      		convertXY_celestial((double)(gw[i].phi_g*rad2deg)-180,
			  90-(double)(gw[i].theta_g*rad2deg),&px[0],&py[0]);
      		fx[0] = (float)px[0];
      		fy[0] = (float)py[0];
      		cpgpt(1,fx,fy,1); 
		i++;
    	}
  }
  cpgsci(1);

  // Plot the pulsar positions  
  for (i=0;i<npsr;i++)
    {
      convertXY_celestial((double)(psr[i].param[param_raj].val[0]*rad2deg)-180,
			  (double)psr[i].param[param_decj].val[0]*rad2deg,&px[0],&py[0]);
      fx[0] = (float)px[0];
      fy[0] = (float)py[0];
      cpgsch(2); cpgsci(2); cpgpt(1,fx,fy,12); cpgsci(1); cpgsch(1);
    }
}

void plotSpectrum(gwgeneralSrc *gw,int *numsources,longdouble tspan, gwgenSpec gwAmps)
{
  float gmaxx,gminx,gmaxy,gminy;
  float maxx,minx,maxy,miny;
  float **fx,**fy,sx[2],sy[2];
  fx=(float **)malloc(4*sizeof(float *));
  fy=(float **)malloc(4*sizeof(float *));
  for (int j=0;j<4;j++) {
	fx[j]=(float *)malloc(numsources[j]*sizeof(float));
	fy[j]=(float *)malloc(numsources[j]*sizeof(float));
  }
  int i;
  i=0;
  gmaxx=gminx=gmaxy=gminy=0.;
  for (int j=0;j<4;j++) {
	if (numsources[j]) {
  		for (int k=0;k<numsources[j];k++) {
      			fx[j][k] = (float)(log10((gw[i].omega_g)/(2.0*M_PI)));
			if (j == 0)
      				fy[j][k] = (float)(log10(pow(gw[i].aplus_g,2)+pow(gw[i].across_g,2)));
			else if (j == 1)
      				fy[j][k] = (float)(log10(pow(gw[i].ast_g,2)));
			else if (j == 2)
      				fy[j][k] = (float)(log10(pow(gw[i].asl_g,2)));
			else if (j == 3)
      				fy[j][k] = (float)(log10(pow(gw[i].avx_g,2)+pow(gw[i].avy_g,2)));
			i++;
    		}
  		minx = TKfindMin_f(fx[j],numsources[j]);
  		maxx = TKfindMax_f(fx[j],numsources[j]);
  		miny = TKfindMin_f(fy[j],numsources[j]);
  		maxy = TKfindMax_f(fy[j],numsources[j]);
		if (!gminx || (minx < gminx))
			gminx=minx;
		if (!gmaxx || (maxx > gmaxx))
			gmaxx=maxx;
		if (!gminy || (miny < gminy))
			gminy=miny;
		if (!gmaxy || (maxy > gmaxy))
			gmaxy=maxy;
	}
  }
  cpgenv(gminx-1,gmaxx+1,gminy,gmaxy,0,1);
  cpglab("log\\d10\\u[f\\dg\\u] (Hz)","log\\d10\\u[A\\d+\\u(f)\\u2\\d+A\\dx\\u(f)\\u2\\d]","");
  for (int j=0;j<4;j++) {
  	cpgsci(2+j);
  	cpgpt(numsources[j],fx[j],fy[j],1);
  }
  cpgsci(1);

  sx[0] = log10(1.0/86400.0); sx[1] = sx[0];
  sy[0] = gminy; sy[1] = gmaxy;
  cpgsci(1); cpgsls(4); cpgline(2,sx,sy); cpgsci(1); cpgsls(1);
  sx[0] = log10(1.0/tspan); sx[1] = sx[0];
  sy[0] = gminy; sy[1] = gmaxy;
  cpgsci(1); cpgsls(4); cpgline(2,sx,sy); cpgsci(1); cpgsls(1);
  
  // Draw line with slope of 2*"alpha"
  //  fx[0] = minx-1;
  //  fy[0] = maxy;
  //  fx[1] = maxx+1;
  //  fy[1] = 2.0*alpha*fx[1]+(fy[0]-2.0*alpha*fx[0]);
  
  sx[0] = gminx-1;
  sx[1] = gmaxx+1;
  for (int j=0;j<4;j++) {
	if (numsources[j]) { 
	if (j == 0) {
  		sy[0] = log10(gwAmps.tensor_amp*gwAmps.tensor_amp/12.0/M_PI/M_PI*pow(pow(10,sx[0]),2*gwAmps.tensor_alpha));
  		sy[1] = log10(gwAmps.tensor_amp*gwAmps.tensor_amp/12.0/M_PI/M_PI*pow(pow(10,sx[1]),2*gwAmps.tensor_alpha));
	} else if (j == 1) {
  		sy[0] = log10(gwAmps.st_amp*gwAmps.st_amp/12.0/M_PI/M_PI*pow(pow(10,sx[0]),2*gwAmps.st_alpha));
  		sy[1] = log10(gwAmps.st_amp*gwAmps.st_amp/12.0/M_PI/M_PI*pow(pow(10,sx[1]),2*gwAmps.st_alpha));
	} else if (j == 2) {
  		sy[0] = log10(gwAmps.sl_amp*gwAmps.sl_amp/12.0/M_PI/M_PI*pow(pow(10,sx[0]),2*gwAmps.sl_alpha));
  		sy[1] = log10(gwAmps.sl_amp*gwAmps.sl_amp/12.0/M_PI/M_PI*pow(pow(10,sx[1]),2*gwAmps.sl_alpha));
	} else if (j == 3) {
  		sy[0] = log10(gwAmps.vl_amp*gwAmps.vl_amp/12.0/M_PI/M_PI*pow(pow(10,sx[0]),2*gwAmps.vl_alpha));
  		sy[1] = log10(gwAmps.vl_amp*gwAmps.vl_amp/12.0/M_PI/M_PI*pow(pow(10,sx[1]),2*gwAmps.vl_alpha));
	}
  	printf("line from %g %g to %g %g\n",sx[0],sy[0],sx[1],sy[1]);
  	cpgsci(2+j); cpgsls(3); cpgline(2,sx,sy); 
	}
  }
  cpgsci(1); cpgsls(1);
}

void plotResiduals(pulsar *psr,longdouble **gwRes,int p,longdouble timeOffset,int plotType)
{
  float px[MAX_OBSN],py[MAX_OBSN];
  float minx,maxx,miny,maxy;
  float yerr1[MAX_OBSN],yerr2[MAX_OBSN];
  int i;
  for (i=0;i<psr[p].nobs;i++)
    {
      px[i] = (float)(psr[p].obsn[i].sat - timeOffset);
      if (plotType==1 || plotType==2)
	py[i] = (float)gwRes[p][i];
      else if (plotType==3)
	py[i] = (float)psr[p].obsn[i].prefitResidual;
      else if (plotType==4)
	py[i] = (float)psr[p].obsn[i].residual;
    }
  TKremovePoly_f(px,py,psr[p].nobs,1);
  for (i=0;i<psr[p].nobs;i++)
    {
      yerr1[i] = py[i] - psr[p].obsn[i].toaErr*1e-6;
      yerr2[i] = py[i] + psr[p].obsn[i].toaErr*1e-6;
    }
  if (plotType==2)
    TKremovePoly_f(px,py,psr[p].nobs,3);
  minx = TKfindMin_f(px,psr[p].nobs);
  maxx = TKfindMax_f(px,psr[p].nobs);
  miny = TKfindMin_f(py,psr[p].nobs);
  maxy = TKfindMax_f(py,psr[p].nobs);
  cpgenv(minx,maxx,miny,maxy,0,1);
  if (plotType==1)
    cpglab("Day","GW residual (s)","Induced residuals due to GW background");
  else if (plotType==2)
    cpglab("Day","GW residual (s)","Induced residuals due to GW background after quadratic removed");
  else if (plotType==3)
    cpglab("Day","GW residual (s)","Pre-fit timing residuals");
  else if (plotType==4)
    cpglab("Day","GW residual (s)","Post-fit timing residuals");
    
  cpgpt(psr[p].nobs,px,py,9);
  if (plotType==3 || plotType==4)
    cpgerry(psr[p].nobs,px,yerr1,yerr2,1);
}

longdouble getTspan(pulsar *psr,int npsr)
{
  longdouble first,last;
  int i,p;
    
  
  first = psr[0].obsn[0].sat;
  last = psr[0].obsn[0].sat;

  for (p=0;p<npsr;p++)
    {
      for (i=0;i<psr[p].nobs;i++)
	{
	  if (first > psr[p].obsn[i].sat)
	    first = psr[p].obsn[i].sat;
	  if (last < psr[p].obsn[i].sat)
	    last = psr[p].obsn[i].sat;
	}
    }

  return last-first;
}

void draw_grid(double start_gl,double end_gl,double start_gb,double end_gb,double gstep,double bstep,int celestialCoords)
{
  double l,b,x,y;
  float plotx[1000],ploty[1000];
  int count=0;
  char str[100];
  cpgenv(start_gl,end_gl,start_gb,end_gb,0,-2);
  cpgsls(4);

  /* Plot lines of latitude */
  for (b=start_gb;b<=end_gb;b+=bstep)
    {
      count=0;
      for (l=start_gl;l<=end_gl;l=l+1.0)
	{
	  if (celestialCoords==1) convertXY_celestial(l,b,&x,&y);
	  /*get_xy(l,b,&x,&y); */
	  plotx[count] = (float)x;
	  ploty[count] = (float)y;
	  /*	  printf("%d %f %f\n",count,plotx[count],ploty[count]); */
	  count++;
	}
      cpgline(count,plotx,ploty);
    }

  /* Plot lines of longitude */
  for (l=start_gl;l<=end_gl;l+=gstep)
    {
      count=0;
      for (b=start_gb;b<=end_gb;b=b+1.0)
	{
	  if (celestialCoords==1) convertXY_celestial(l,b,&x,&y);
	  /*	  get_xy(l,b,&x,&y); */
	  plotx[count] = (float)x;
	  ploty[count] = (float)y;
	  count++;
	}
      if (l==-180 || l==180)
	cpgsls(1);
      else
	cpgsls(4);
      cpgline(count,plotx,ploty);
    }
  

  /* Label axes */
  cpgsci(3);
  for (l=0;l<360;l+=gstep)
    {
      if (celestialCoords==1) convertXY_celestial(l-180,-45,&x,&y);

      /*      if (celestialCoords==1)
	get_xy(l+180.0,-45,&x,&y);
      else
      get_xy(l,-45,&x,&y); */
      if (l!=180.0 || celestialCoords==1)
	{
	  if (celestialCoords==0 || l!=0)
	    {
	      if (celestialCoords==0) sprintf(str,"%.0f\\uo\\d",l);
	      else sprintf(str,"%.0f\\uh\\d",l/360.0*24.0);
	      cpgptxt((float)x,(float)y,0,0.5,str);
	    }
	}
    }
  for (b=-60;b<=60;b+=bstep)
    {
      if (celestialCoords==1) convertXY_celestial(-180,b,&x,&y);
      /*      get_xy(180,b,&x,&y); */
      if (b>0)
	{
	  sprintf(str,"+%.0f\\uo\\d",b);
	  cpgptxt((float)x,(float)y,0,1.0,str);
	}
      else if (b==0)
	{
	  sprintf(str,"%.0f\\uo\\d",b);
	  cpgptxt((float)x-2,(float)y,0,1.0,str);
	}
      else
	{
	  sprintf(str,"%.0f\\uo\\d",b);
	  cpgptxt((float)x,(float)y-7,0,1.0,str);
	}
    }
  cpgsci(1);
  cpgsls(1);
}

/* Convert from RAJ, DECJ to x,y using Aitoff projection */
void convertXY_celestial(double raj,double decj,double *retx,double *rety)
{
  double sa;
  double r2deg = 180.0/M_PI;
  double alpha2,delta;
  double r2,f,cgb,denom;
  double x_ret,y_ret;

  sa = raj;
  alpha2 = sa/(2*r2deg);
  delta = decj/r2deg;   

  r2 = sqrt(2.);    
  f = 2*r2/M_PI;    

  cgb = cos(delta);    
  denom =sqrt(1. + cgb*cos(alpha2));

  x_ret = cgb*sin(alpha2)*2.*r2/denom;
  y_ret = sin(delta)*r2/denom;

  x_ret = x_ret*r2deg/f;
  y_ret = y_ret*r2deg/f;

  *retx = x_ret;
  *rety = y_ret;
}
const char * plugVersionCheck = TEMPO2_h_VER;
