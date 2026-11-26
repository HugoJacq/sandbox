/**
# Wave breaking with stratification (multilayer solver)

Initial significant wave height is maintained using a wind forcing, see
[wind_forcing.h](https://basilisk.fr/sandbox/hugoj/lib/wind_forcing.h) and
the associated [test case](https://basilisk.fr/sandbox/hugoj/test_windinput/test_forcing/linear_wave_wind_input.c).

Temperature is a passive tracer
*/

double rho = 1000;
double pw = 0.;
#if HOLD_FORCING
  //#define a_baro(f,k) (1.0*(-G*(f[k]-f[k-1])/Delta - (pw/rho)*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta)))
  // a_baro is defined in hydro.h
  #define a_baro_p(f,k) (- (pw/rho)*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta))
#endif

const double g_ = 9.81 [1,-2];        // [m.s-2] Gravity

#include "layered/hydro.h"
#include "layered/nh-staggered.nh" // Bruno's nh.h copied in src/layered
//#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files
#include "hugoj/lib/spectrum.h"       // Initial conditions generation
#include "hugoj/lib/diffusionH.h"     // Neumann at top and bot for diffusion

/**
Dimensions : [Length, Time, Temperature, Energy, mass]
*/
char file_out[20]="out.nc";

scalar T;
double strat, Tini, Trand, Ts, t_iniT;
double ppw;
double cp, betaT, diff_T;
int coeff_kpL0 = 5;               
double kp, Tp;         
double NT0, tend, tstart_diag, tend_diag, smalltime;

/** to output every dtout, a value needs to be assigned before main ! */
double dtout=1.0;

double Re, Pr, Ri;                     
double h_;
double * T_inip;

#define Tini(z) strat/(g_*betaT)*z + Ts  // Initial temperature profile

/** A function to compute the energy of the simulation */
void energy_stats (double *mean, double *variance, double *pe, double *ke)
{
  double sum = 0.;
  foreach (reduction(+:sum))
    sum += eta[];
  double m = sum/pow (N, dimension);

  double gpe = 0., ke2 = 0.;
  foreach (reduction(+:gpe) reduction(+:ke2)) {
    gpe += sq (eta[] - m)*dv();
    foreach_layer() {
      double norm2 = sq (w[]);
      foreach_dimension()
        norm2 += sq (u.x[]);
      ke2 += norm2*h[]*dv();
    }
  }
  *mean = m;
  *variance = gpe/sq(L0);  // <(eta-m)^2> since sum of dv() = L0^dimension
  *pe = 0.5*rho*G*gpe;
  *ke = 0.5*rho*ke2;
}


/** wind forcing controler is added with this [file](https://basilisk.fr/sandbox/hugoj/lib/wind_forcing.h) */
#if HOLD_FORCING
#include "hugoj/lib/wind_forcing.h"
#endif

// Note
// For now I have to use this face_fields event to add the forcing.
// When i'll be able to use a vector inside the a_baro macro I will 
// // be able to use the original machinery with the a_baro macro overloaded.
// #if FULL_FORCING
// #else
#if HOLD_FORCING 
#if 1
event face_fields (i++,first) {
  foreach_face (x) {
    double ax = a_baro_p (eta, 0);
    foreach_layer() {
      ha.x[] += hf.x[]*ax;
    }
  }
}
#endif
#endif // HOLD_FORCING 
// #endif // FULL_FORCING


int main(int argc, char *argv[])  
{


/** ## Parameters */
  Re = 1000;                     // Reynolds number Re = sqrt(g*lambda**3)/nu
  Pr = 1.0;                      // Prandtl number
  
  // -> Initial conditions
  strat = 0.002 [0,-2];       // [s-2] N^2 stratification
  Ts = 20. [0,0,1];              // [K] Surface temperature (arbitrary)
  P = 0.2 [1, -1];               // energy level (estimated so that kpHs is reasonable)
  coeff_kpL0 = 5 [];             // kpL0 = coeff_kpL0 * pi
  N_mode = 32 [];                // Number of modes in wavenumber space
  N_power = 5 [];                // directional spreading coeff
  F_shape = 0 [];                // shape of the initial spectrum
  /** Forcing */
  ppw = 1. [];
  /** Domain definition */
  N = 256;                       // Horizontal grid number
  nl = 64;                       // Number of layers
  L0 = 200.;                     // Horizontal domain size
  h_ = 20.0 [1];                  // [m] depth of water
  /** Runtime parameters */
  NT0 = 2. [0,1];               // (x Tp) end time of simulation
  t_iniT = 0.1;                  // (x Tp) time where temperature is initialised
  dtout = 0.1 [0,1];             // (x Tp) dt for output in netcdf
  smalltime = 1e-10 [0,1];       // [s] small time increment
  /* physical properties of the fluid */
  betaT = 2e-4 [0,0,-1];         // [K-1] Thermal expansion coeff for water
  Tini = 20. [0,0,1];            // [°C] Reference temperature
  
  
  G = g_;                       // gravity
  theta_H = 0.5;                // implicit parameter

  kp = 2*PI * coeff_kpL0 / L0; // kpL=coeff x 2pi
  double omegap = sqrt(g_*kp);
  double cp = omegap/kp;
  nu = sqrt(g_/kp)*2*PI/kp/Re; // Re = c*lambda/nu
  Tp = 2*PI/sqrt(g_*kp);
  diff_T = nu/Pr;
  
  tend = NT0*Tp;
  t_iniT = t_iniT*Tp;
  dtout = dtout*Tp;

  T_inip = (double*)calloc(nl, sizeof(double));

  // Boundary conditions
  origin (-L0/2., -L0/2.);
  periodic (top);
  periodic (left);

  run();
}


/** ## Initial condition */
event init(i =  0) {

  if (!restore ("restart")) {
    /** We read a spectrum using spectrum.h */
    T_Spectrum spectrum;
    // spectrum = read_spectrum(N_mode);
    spectrum = spectrum_gen_linear(N_mode, N_power, L0, P, kp);  

    /** set eta */
    initial_condition_wave_fft (eta, spectrum, N);

    /** set h*/
    geometric_beta (1./3., true); // if !=0, varying layer thickness
    foreach(cpu) {
      zb[] = -h_;
      double H = eta[] - zb[];
      foreach_layer() {
        h[] = H*beta[point.l];
      } 
    }

    /** set a temporary T_ini field that will be used to initialse T at a later
      time */
    double z = - h_;
    for (int ll=0; ll<nl; ll++){
      z += h_*beta[ll]/2.;
      T_inip[ll] = Tini(z);
      z += h_*beta[ll]/2.;
    }

    /** set currents */
    initial_condition_u_fft (u, spectrum, -h_, N);

    /** add a passive tracer */
    T = new scalar[nl];
    tracers = list_append (tracers, T);
    reset ({T}, 0.);
    fprintf (stderr,"Done initialization!\n");
    free_spectrum(&spectrum);
  }
  else {
    // Restarting from file 'restart'
    // Warning: no check is done to ensure that you restart with the same
    // resolution/parameters ! Basilisk just reload what is in the dump file.
    geometric_beta (1./3., true); // when restarting, remember to specify the grid mapping method
    dtmax = 0.01;
    dt = dtnext (dtmax);
    fprintf (stderr, "Restarting ...\n");
  }
  create_nc({zb, h, u, w, eta, T}, file_out);
}

/** vertical diffusion on T and u,w */
event viscous_term (i++; t<=tend)
{
  foreach()
    vertical_diffusion2 (point, h, T, dt, diff_T, 0., 0.);
  horizontal_diffusion ({T}, diff_T, dt);
  horizontal_diffusion ({u, w}, nu, dt);
}

/** ## Initialise T after wave spinup */
event initT (t = t_iniT){
  foreach(cpu){
    foreach_layer(){
      T[] = T_inip[point.l];
    }
  }
}

/** dump outputs */
event dumpnc (t += dtout){
  fprintf(stderr, "%g\n", t);
  write_nc();
  
  #if 0
  char dname[100];
  sprintf (dname, "dump_t%g", t);
  dump(dname);
  #endif
}

event final_dump (t=end){
  char dname[100];
  sprintf (dname, "dump_t%g", t);
  dump(dname);
}

event cleanup (t = end, last)
{
  delete ({T});
  free(T_inip);
}

/**
## Notes

- for now, we use the acceleration event hook but this is not clean. When
the overloading method of the a_baro macro works I will change that.



*/
