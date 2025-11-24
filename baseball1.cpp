///
/// Starter template for first baseball problem
/// Solve for the initial speed of the pitch given the initial parameters
/// xend : distance to home plate [18.5] m
/// z0 : height of release of ball [1.4] m
/// theta0 : angle of release above horizontal [1] degree
///
///  Do not change the interface for running the program
///  Fill in the value of vPitch in the print statement with your solution
///  at the end of main()
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

struct Params {
  double g;   // acceleration [m/s^2]
  double m;   // mass of object [kg], nb proj. In vacuum funcs do not depend on the mass
  double d;   // m diameter of ball
  double b;   // b,c params for air resistance
  double c;
};

void baseball_derivs(double t, double y[], double dydt[], void *params){
  Params *p = (Params*)params;
  double vx = y[2];
  double vz = y[3];
  double v = sqrt(vx*vx + vz*vz);

  if (v< 1e-8)v = 1e-8;

  double Fmag = (p->b * p->d)* v +( p->c* p->d* p->d) *v*v;
  double drag_x = -Fmag* (vx/v);
  double drag_z = -Fmag *(vz/v);

  dydt[0] = vx;
  dydt[1] = vz;
  dydt[2] = drag_x/ p->m;
  dydt[3] = -p->g + drag_z/p->m;
}

void RK4Step(double &t, double y[], int nvar, double dt, 
              void (*derivs)(double, double[], double[], void*), void *params) {
    double k1[nvar], k2[nvar], k3[nvar], k4[nvar];
    double ytemp[nvar];
    double dydt[nvar];
    
    derivs(t, y, dydt, params);
    for(int i=0; i<nvar; i++) {
        k1[i] = dt * dydt[i];
        ytemp[i] = y[i] + 0.5*k1[i];
    }
    
    derivs(t+0.5*dt, ytemp, dydt, params);
    for(int i=0; i<nvar; i++) {
        k2[i] = dt * dydt[i];
        ytemp[i] = y[i] + 0.5*k2[i];
    }
    
    derivs(t+0.5*dt, ytemp, dydt, params);
    for(int i=0; i<nvar; i++) {
        k3[i] = dt * dydt[i];
        ytemp[i] = y[i] + k3[i];
    }
    
    derivs(t+dt, ytemp, dydt, params);
    for(int i=0; i<nvar; i++) {
        k4[i] = dt * dydt[i];
    }
    
    for(int i=0; i<nvar; i++) {
        y[i] += (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
    }
    t += dt;
}

double compute_z_at_xend(double v0, double theta_rad, double z0, double xend, Params &pars){
  const int nvar = 4;
  double y[nvar];

  y[0] = 0.0;
  y[1] = z0;
  y[2] = v0*cos(theta_rad);
  y[3] = v0*sin(theta_rad);

  double t = 0.0;
  double dt = 0.001;

  void *p_par = (void *) &pars;

  while (y[0] < xend && y[1]>0.0){
    RK4Step(t, y, nvar, dt, baseball_derivs, p_par);
  }
  return y[1];
}

int main(int argc, char **argv){

  // examples of parameters
  Params pars;
  pars.g=9.81;
  pars.m=0.145;    
  pars.d=0.075;   
  pars.b=1.6e-4;  
  pars.c=0.25;

  double xend=18.5;       // meters to plate
  double z0=1.4;             // height of release [m]
  double theta0=1;         // angle of velocity at release (degrees)
                                      // convert to radians before using!
  bool showPlot=false;    // keep this flag false by default
  
  // allow changing the parameters from the command line
  int c;
  while ((c = getopt (argc, argv, "x:z:t:p")) != -1)
    switch (c) {
    case 'x':
      xend = atof(optarg);
      break;
    case 'z':
      z0 = atof(optarg);
      break;
    case 't':
      theta0 = atof(optarg);
      break;
    case 'p':
      showPlot=true;
      break;
    case '?':
      fprintf (stderr, "Unknown option `%c'.\n", optopt);
    }
  TApplication theApp("App", &argc, argv); // init ROOT App for displays


  // write code to solve for vPitch here
  double theta_rad = theta0 *M_PI/180.0;

  double z_target = 0.9;

  double v_low =0.0;
  double v_high =100.0;
  double tolerance = 0.001;

  while (v_high-v_low>tolerance){
    double v_mid = 0.5*(v_low+v_high);
    double z_mid = compute_z_at_xend(v_mid, theta_rad, z0, xend, pars);

    if (z_mid>z_target){
      v_high = v_mid;
    } else {
      v_low = v_mid;
    }
  }

  double vPitch = 0.5*(v_low+v_high);
  

  // do not change these lines
  printf("********************************\n");
  printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n",xend,z0,theta0);
  printf("v_pitch = %lf m/s\n",vPitch);
  printf("********************************\n");

  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}

