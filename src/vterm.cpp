///
/// @file 
/// @brief Test of Runge-Kutta solver for series of ODEs
/// @author Bob Hirosky
/// @date 31 Dec 2019 
/// 
/// Use the Rk4 solver for coupled ODEs to solve for projectile 
/// motion with air resistance
///
/// Definition of our variables
/// x    = time <br>
/// y[0] = position along i axis  ; f_ri = dri/dt => velocity along i axis  <br>
/// y[1] = velocity along i axis  ; f_vi = dvi/dt => acceleration along i axis <br>
/// y[2] = position along j axis  ; f_rj = drj/dt => velocity along j axis <br>
/// y[3] = velocity along j axis  ; f_vj = dvj/dt => acceleration along j axis <br>


#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
using namespace std;


struct Params {
  double g;   ///< acceleration [m/s^2]
  double m;   ///< mass of object [kg], nb proj. In vacuum funcs do not depend on the mass
  double air_k;  ///< constant for air resistance. mass DOES matter with air resistance
} ;

// functions to describe simple projectile motion
// here use use ri,rj,rk to define directions to prevent confusion with
// standard ODE notation, where x=independent variable, \vec y=dependent variable(s)


/// \brief Change in position along \f$\hat i\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
double f_ri(double x, const vector<double> &y, void *params=0){ 
  (void) x;   // prevent unused variable warning
  return y[1];
}

/// \brief Change in velocity along  \f$\hat i\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
double f_vi(double x, const vector<double> &y, void *params=0){ 
  (void) x;
  Params *p = (Params*)params;
  return -p->air_k * sqrt(y[1]*y[1] + y[3]*y[3]) * y[1] / p->m;
  // return 0;  // if no air, no forces/acceleration along i direction in this problem
}

/// \brief Change in position along \f$\hat j\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
///
/// Air resistance model: F= \f$k v^2\f$
///
double f_rj(double x, const vector<double> &y, void *params=0){  
  (void) x;   // prevent unused variable warning
  return y[3];
}

/// Change in velocity along  \f$\hat j\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
double f_vj(double x, const vector<double> &y, void *params=0){  
  (void) x;
  Params *p = (Params*)params;
  return -p->air_k * sqrt(y[1]*y[1] + y[3]*y[3]) * y[3] / p->m - p->g;
  // return -g;    // if no air constant acceleration along -j direction: F/m = -g
}

double f_ri_noair(double x, const vector<double> &y, void *params=0){ 
  (void) x;   // prevent unused variable warning
  (void) params;
  return y[1];
}

double f_vi_noair(double x, const vector<double> &y, void *params=0){ 
  (void) x;
  (void) y;
  (void) params;
  return 0;
}

double f_rj_noair(double x, const vector<double> &y, void *params=0){  
  (void) x;   // prevent unused variable warning
  (void) params;
  return y[3];
}

double f_vj_noair(double x, const vector<double> &y, void *params=0){  
  (void) x;
  (void) y;
  Params *p = (Params*)params;
  return -p->g;
}

/// \brief Stopping condition
/// \param[in] x independent variable
/// \param[in] y dependent variables
///
/// Returns 0(1) to flag continuation(termination) of calculation 
double f_stop(double x, const vector<double> &y, void *params=0){
  (void) x;
  (void) params;
  if (y[2]<0) return 1;  // stop calulation if the current step takes height to negative value
  return 0;  // continue calculation
}

double total_energy(const vector<double> &y, const Params &p) {
  double ke = 0.5 * p.m * (y[1]*y[1] + y[3]*y[3]);
  double pe = p.m * p.g * y[2];
  return ke + pe;
}

double find_terminal_velocity(const Params &pars, double dt=0.01) {
  vector<pfunc_t> v_fun(4);
  v_fun[0]=f_ri; v_fun[1]=f_vi; v_fun[2]=f_rj; v_fun[3]=f_vj;
  
  vector<double> y(4);
  y[0]=0; y[1]=0;
  y[2]=10000; y[3]=0;
  
  void *p_par = (void*) &pars;
  double x=0;
  double tolerance = 1e-4;
  double check_interval = 1.0;
  double v_prev = 0;
  double v_current;
  int steps_per_check = (int)(check_interval/dt);
  int max_checks = 100;
  
  for (int check = 0; check < max_checks; check++) {
    double xmax = x + check_interval;
    auto tgN = RK4SolveN(v_fun, y, steps_per_check, x, xmax, p_par);
    
    v_current = sqrt(y[1]*y[1] + y[3]*y[3]);
    
    if (check > 0 && fabs(v_current - v_prev) / v_current < tolerance) {
      return v_current;
    }
    
    v_prev = v_current;
  }
  
  return sqrt(y[1]*y[1] + y[3]*y[3]);
}

/// \brief Use RK4 method to describe simple projectile motion.
int main(int argc, char **argv){

  // setup default parameters
  Params pars;
  pars.g=9.81;
  pars.m=10.0;
  pars.air_k=0.1;
  void *p_par = (void*) &pars;

  double theta=45;   // initial angle degrees
  double v0=100;     // m/s

  bool energy_test = false;
  bool vterm_test = false;
  bool mass_study = false;
  
  int c;
  while ((c = getopt (argc, argv, "v:t:m:k:eVM")) != -1)
    switch (c) {
    case 'v':
      v0 = atof(optarg);
      break;
    case 't':
      theta = atof(optarg);
      break;
    case 'm':
      pars.m = atof(optarg);
      break;
    case 'k':
      pars.air_k = atof(optarg);
      break;
    case 'e':
      energy_test = true;
      break;
    case 'V':
      vterm_test = true;
      break;
    case 'M':
      mass_study = true;
      break;
    case '?':
      fprintf (stderr, "Unknown option `%c'.\n", optopt);
      return 1;
    }
  
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  // ******************************************************************************
  // ** this block is useful for supporting both high and std resolution screens **
  UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height  
  //UInt_t dw = gClient->GetDisplayWidth();
  UInt_t dw = 1.1*dh;
  // ******************************************************************************

  if (energy_test) {
    cout << "\n Energy Conservation Test (No Air Resistance) " << endl;
    
    vector<pfunc_t> v_fun_noair(4);
    v_fun_noair[0]=f_ri_noair;
    v_fun_noair[1]=f_vi_noair;
    v_fun_noair[2]=f_rj_noair;
    v_fun_noair[3]=f_vj_noair;
    
    vector<int> step_sizes = {50, 100, 200, 500, 1000, 2000};
    ofstream energy_file("energy_conservation.dat");
    energy_file << "# nsteps dt energy_error\n";
    
    for (int nsteps : step_sizes) {
      vector<double> y(4);
      y[0]=0; y[1]=v0*cos(theta*M_PI/180);
      y[2]=0; y[3]=v0*sin(theta*M_PI/180);
      
      double E0 = total_energy(y, pars);
      double x=0, xmax=20;
      
      auto tgN = RK4SolveN(v_fun_noair, y, nsteps, x, xmax, p_par, f_stop);
      double Ef = total_energy(y, pars);
      double dE = fabs((Ef - E0) / E0);
      double dt = xmax / nsteps;
      
      cout << "nsteps=" << nsteps << ", dt=" << dt << "s, dE/E0=" << dE << endl;
      energy_file << nsteps << " " << dt << " " << dE << "\n";
    }
    energy_file.close();
  }

   if (vterm_test) {
    cout << "\n Terminal Velocity Determination " << endl;
    cout << "Mass: " << pars.m << " kg" << endl;
    cout << "Drag coefficient: " << pars.air_k << " kg/m" << endl;
    
    double vt = find_terminal_velocity(pars);
    cout << "Terminal velocity: " << vt << " m/s" << endl;
    
    double vt_theory = sqrt(pars.m * pars.g / pars.air_k);
    cout << "Theoretical vt: " << vt_theory << " m/s" << endl;
    cout << "Difference: " << fabs(vt - vt_theory) / vt_theory * 100 << "%" << endl;
  }

   if (mass_study) {
    cout << "\n Terminal Velocity vs Mass Study " << endl;
    
    ofstream mass_file("vterm_vs_mass.dat");
    mass_file << "# mass(kg) vt_simulated(m/s) vt_theoretical(m/s)\n";
    
    vector<double> masses, vt_sim, vt_theo;
    
    for (double log_m = -3; log_m <= 1.01; log_m += 0.2) {
      double mass = pow(10, log_m);
      Params p_temp = pars;
      p_temp.m = mass;
      
      double vt = find_terminal_velocity(p_temp, 0.005);
      double vt_theory = sqrt(mass * p_temp.g / p_temp.air_k);
      double error_percent = fabs(vt-vt_theory)/vt_theory *100;
      
      cout << "m=" << mass << " kg, vt=" << vt << " m/s (theory: " 
           << vt_theory << " m/s, error" << error_percent <<"%)"<< endl;
      
      masses.push_back(mass);
      vt_sim.push_back(vt);
      vt_theo.push_back(vt_theory);
      
      mass_file << mass << " " << vt << " " << vt_theory << "\n";
    }
    mass_file.close();

     TCanvas *c1 = new TCanvas("c1","Terminal Velocity vs Mass",dw,dh);
    c1->SetLogx();
    c1->SetLogy();
    
    TGraph *gr_sim = new TGraph(masses.size(), &masses[0], &vt_sim[0]);
    TGraph *gr_theo = new TGraph(masses.size(), &masses[0], &vt_theo[0]);
    
    gr_sim->SetMarkerStyle(20);
    gr_sim->SetMarkerColor(kBlue);
    gr_sim->SetLineColor(kBlue);
    gr_sim->SetTitle("Terminal Velocity vs Mass;Mass (kg);Terminal Velocity (m/s)");
    
    gr_theo->SetMarkerStyle(24);
    gr_theo->SetMarkerColor(kRed);
    gr_theo->SetLineColor(kRed);
    gr_theo->SetLineStyle(2);
    
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr_sim);
    mg->Add(gr_theo);
    mg->SetTitle("Terminal Velocity vs Mass;Mass (kg);Terminal Velocity (m/s)");
    mg->Draw("ALP");
    
    TLegend *leg = new TLegend(0.15, 0.7, 0.45, 0.85);
    leg->AddEntry(gr_sim, "Simulated", "lp");
    leg->AddEntry(gr_theo, "Theoretical", "lp");
    leg->Draw();
    
    c1->SaveAs("vterm.pdf");
    
  }
  
  if (!energy_test && !vterm_test && !mass_study) {
    cout << "\n Standard Projectile Motion " << endl;
    vector<pfunc_t> v_fun(4);
    v_fun[0]=f_ri; v_fun[1]=f_vi; v_fun[2]=f_rj; v_fun[3]=f_vj;
    
    vector<double> y(4);
    y[0]=0; y[1]=v0*cos(theta*M_PI/180);
    y[2]=0; y[3]=v0*sin(theta*M_PI/180);
    
    cout << "Vinit: " << v0 << " m/s" << endl;
    cout << "Angle: " << theta << " deg" << endl;
    cout << "(vx,vy) " << y[1] << " , "  <<  y[3] << " m/s" << endl;
    cout << "Mass: " << pars.m << " kg" << endl;
    cout << "Drag: " << pars.air_k << " kg/m" << endl;
    
    double x=0, xmax=20;
    int nsteps=200;
    
    auto tgN = RK4SolveN(v_fun, y, nsteps, x, xmax, p_par, f_stop);

    cout << "Number of trajectory points: " << tgN[2].GetN() << endl;
    cout << "Flight time: " << x << " seconds" << endl;
    cout << "Final position: (" << y[0] << ", " << y[2] << ")" << endl;
    cout << "Final velocity = " << sqrt(y[1]*y[1]+y[3]*y[3]) << " m/s" << endl;
    

    TFile *tf=new TFile("vterm.root","recreate");
    for (unsigned i=0; i<v_fun.size(); i++){
      tgN[i].Write();
    }
    tf->Close();

  }

  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}
