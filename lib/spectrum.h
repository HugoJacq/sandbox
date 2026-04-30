/*
 Code translation from python to C.
 Original code : Jiarong Wu,
 https://github.com/jiarong-wu/multilayer_breaking/blob/main/specgen/specgen.py

g_ has to be defined before #include spectrum.h !

*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "hugoj/lib/interpolate.h"
//#include <gsl/gsl_rng.h>

#define PI 3.14159265358979323846

typedef struct {
  int N_mode;
  //double *kmod;
  double *kx;
  double *ky;
  //double *F_kmod;
  double *F_kxky;
  double *phase;
  double *omega;
} T_Spectrum;

void cart2pol(double x, double y, double *rho, double *phi) {
  *rho = sqrt(x * x + y * y);
  *phi = atan2(y, x);
}

void pol2cart(double rho, double phi, double *x, double *y) {
  *x = rho * cos(phi);
  *y = rho * sin(phi);
}


double randInRange(double min, double max)
{
  return min + (rand() / (RAND_MAX+1.0) * (max - min));
}



// Some common spectra
double spectrum_PM(double P, double kp, double kmod) {
  // Note: P here is in fact P/sqrt(g) in the PM spectrum equation
  return P * pow(kmod, -2.5) * exp(-1.25 * pow(kp / kmod, 2.0));
}
double spectrum_JONSWAP(double alpha, double kp, double kmod) {
  return alpha * pow(kmod, -3.0) * exp(-1.25 * pow(kp / kmod, 2.0));
}
double spectrum_Gaussian(double G, double span, double kp, double kmod) {
  return (G / span) * exp(-0.5 * pow((kmod - kp) / span, 2.0));
}

T_Spectrum spectrum_gen_linear(int N_mode, int N_power, double L, double P,
                               double kp, double thetam=0.) 
{

  /* The function to generate a kx-ky spectrum based on uni-directional
      spectrum and a cos^N (theta) directional spreading.
        Arguments:
            shape: the spectrum shape (a function)
            N_mode: # of modes (Right now it's N_mode for kx and N_mode+1 for
                    ky; has to much what's hard-coded in the spectrum.h
                    headerfile.
            N_power: directional spectrum spreading coefficient
            L: physical domain size
            P: energy of the wave field (dimension=velocity)
            kp: peak wavenumber
            thetam: midline direction (rad, positive anticlockwise, along x = 0.)
  */

  int N_kmod = 128; // Uniform grid in kmod and ktheta, can be finer than N_mode
  int N_theta = 128;
  //double thetam = 0.; // midline direction
  double theta[N_theta];
  double dtheta;
  double Dtheta[N_theta];               // for directional spectrum
  double sum = 0.0;                     // for Dtheta normalization
  double F_kmodtheta[N_kmod * N_theta] ; // directional spectrum
  double kmod[N_kmod];
  double F_kmod[N_kmod];
  int Ntmode = 2*N_mode + 1;

  T_Spectrum spectrum;
  spectrum.N_mode = N_mode;
  spectrum.kx = (double *)calloc(Ntmode, sizeof(double));
  spectrum.ky = (double *)calloc(Ntmode, sizeof(double));  
  spectrum.F_kxky = (double *)calloc(Ntmode*Ntmode, sizeof(double));
  spectrum.phase = (double *)calloc(Ntmode*Ntmode, sizeof(double));
  spectrum.omega = (double *)calloc(Ntmode*Ntmode, sizeof(double));

  // building kmod
  for (int i = 0; i < N_kmod; ++i) {
    kmod[i] =
        2 * PI / L + 1.0 * i / (N_kmod - 1) * (1.41 * 100 * 2 - 2) * PI / L;
  }
  // build Dtheta
  for (int i = 0; i < N_theta; ++i) {
    theta[i] = -PI + 2.0*PI*i / (N_theta - 1);
    Dtheta[i] = fabs(pow(cos(theta[i] - thetam), N_power));  
  }
  // Normalizing Dtheta
  dtheta = theta[1] - theta[0];
  for (int i = 0; i < N_theta - 1; ++i) {
    sum = sum + dtheta * 0.5 * (Dtheta[i] + Dtheta[i+1]); // trapezoid integ
  }
  for (int i = 0; i < N_theta; ++i) {
    Dtheta[i] = Dtheta[i] / sum;
  }

  // Pick the spectrum shape : for now PM only
  // TODO: add more shapes ?
  for (int i = 0; i < N_kmod; ++i) {
    F_kmod[i] = spectrum_PM(P, kp, kmod[i]);
  }
  for (int ik = 0; ik < N_kmod; ++ik) {
    for (int itt = 0; itt < N_theta; ++itt) {
      F_kmodtheta[ik * N_theta + itt] =
          F_kmod[ik] * Dtheta[itt] / kmod[ik];
      // Notice the normalize by kmod !
    }
  }

  // Uniform grid in kx ky
  for (int i = 0; i < Ntmode; ++i) {
    spectrum.kx[i] = -2*PI*N_mode/L + 2*PI/L*i;
  }
  for (int i = 0; i < Ntmode; ++i) {
    //spectrum.ky[i] = 2 * PI / L * (i - N_mode / 2);
    spectrum.ky[i] = -2*PI*N_mode/L + 2*PI/L*i;
  }

  // interp F_kmodtheta on kx,ky grid
  double rho, phi;
  double localkx,localky;
  for (int ix = 0; ix < Ntmode; ++ix) {
    for (int iy = 0; iy < Ntmode; ++iy) {
      // first we get polar coords
      cart2pol(spectrum.kx[ix], spectrum.ky[iy], &rho, &phi);

      // Log out-of-range coordinates to diagnose extrapolation
      // if (rho < kmod[0] || rho > kmod[N_kmod-1])
      //   fprintf(stderr, "rho out of range: %f (kmod: %f to %f)\n",
      //           rho, kmod[0], kmod[N_kmod-1]);
      // if (phi < theta[0] || phi > theta[N_theta-1])
      //   fprintf(stderr, "phi out of range: %f (theta: %f to %f)\n",
      //           phi, theta[0], theta[N_theta-1]);

      // then interp at these coords
      spectrum.F_kxky[ix*Ntmode + iy] = 2*interp_lin(
        kmod, theta, N_kmod, N_theta, rho, phi, F_kmodtheta);

      // we remove the negative kx values
      localkx = spectrum.kx[ix]*cos(thetam) + spectrum.ky[iy]*sin(thetam);
      localky = - spectrum.kx[ix]*sin(thetam) + spectrum.ky[iy]*cos(thetam);
      if (localkx < 0.){
        spectrum.F_kxky[ix*Ntmode + iy] = 0.; // and x2 the half plane
      }
      //
      // we remove the center point too, not defined behavior
      if (ix==N_mode && iy==N_mode){
        //fprintf(stderr, "i'm in ! ix=%d iy=%d\n", ix, iy);
        spectrum.F_kxky[ix*Ntmode + iy] = 0.;
      }
      // Uncomment this if F_kxky is < 0.
      // if (spectrum.F_kxky[ix*Ntmode + iy] < 0.){
      //   fprintf(stderr,"i=%d j=%d %f \n", ix,iy, spectrum.F_kxky[ix*Ntmode + iy]);
      // }
    }
  }

  // Random phase
  int RANDOM=2;
  srand(RANDOM); // We can seed it differently for different runs
  int index = 0;
  double k = 0;
  for (int i=0; i<Ntmode; i++) {
    for (int j=0; j<Ntmode; j++) {
      //index = j*N_mode + i;
      index = i*Ntmode + j;
      k = sqrt(sq(spectrum.kx[i]) + sq(spectrum.ky[j]));
      spectrum.omega[index] = sqrt(g_*k); // we use linear dispersion relation 
      spectrum.phase[index] = randInRange (0, 2.*PI); // random phase in [0,2pi]
    }
  }
  return spectrum;
}

T_Spectrum read_spectrum(int N_mode) {
  /* This function reads a file and extract the spectrum from it.
      The spectra can be generated by spec_gen_linear or another function,
      given that it can be read by 'read_spectrum'

      We look for files like
      F_kxky, kx, ky, omega, phase

      Note: 

      - the length of the array in the files has to match with N_mode ! 

      - if the main program is run using mpi, this read_spectrum function will be
      used by all procs. As its read only, no problem. 

      - we read/write omega and phase. Doing this ensure that the same random number has been used.
      We could give this task to the main proc. 
      
   */

  T_Spectrum spectrum;
  int Ntmode=2*N_mode+1;
  spectrum.N_mode = N_mode;
  spectrum.kx = (double *)calloc(Ntmode, sizeof(double));
  spectrum.ky = (double *)calloc(Ntmode, sizeof(double));  
  spectrum.F_kxky = (double *)calloc(Ntmode*Ntmode, sizeof(double));
  spectrum.phase = (double *)calloc(Ntmode*Ntmode, sizeof(double));
  spectrum.omega = (double *)calloc(Ntmode*Ntmode, sizeof(double));
  int length1D;
  int length2D;
  length1D = Ntmode;
  length2D = Ntmode*Ntmode;
  char filename[100];

  // Read F_kxky
  double * a = (double*) malloc (sizeof(double)*length2D);
  sprintf (filename, "F_kxky");
  FILE * fp = fopen (filename, "rb");
  fread (a, sizeof(double), length2D, fp);
  for (int i=0;i<length2D;i++) {
    spectrum.F_kxky[i] = a[i];
  }
  fclose (fp);

  // Read kx
  double * b1 = (double*) malloc (sizeof(double)*length1D);
  sprintf (filename, "kx");
  FILE * fp2 = fopen (filename, "rb");
  fread (b1, sizeof(double), length1D, fp2);
  for (int i=0;i<length1D;i++) {
    spectrum.kx[i] = b1[i];
  }
  fclose (fp2);

  // Read ky
  double * b2 = (double*) malloc (sizeof(double)*length1D);
  sprintf (filename, "ky");
  FILE * fp3 = fopen (filename, "rb");
  fread (b2, sizeof(double), length1D, fp3);
  for (int i=0;i<length1D;i++) {
    spectrum.ky[i] = b2[i];
  }
  fclose (fp3);

  // Read omega
  double * b3 = (double*) malloc (sizeof(double)*length2D);
  sprintf (filename, "omega");
  FILE * fp4 = fopen (filename, "rb");
  fread (b3, sizeof(double), length2D, fp4);
  for (int i=0;i<length2D;i++) {
    spectrum.omega[i] = b3[i];
  }
  fclose (fp4);

  // Read phase
  double * b4 = (double*) malloc (sizeof(double)*length2D);
  sprintf (filename, "phase");
  FILE * fp5 = fopen (filename, "rb");
  fread (b4, sizeof(double), length2D, fp5);
  for (int i=0;i<length2D;i++) {
    spectrum.phase[i] = b4[i];
  }
  fclose (fp5);

  return spectrum;
}




// Surface elevation following the linear wave theory
// (strictly speaking we look for a solution that is the sum of linear modes, no
// need for the  hypothesis of linear theory)
double wave(double x, double y, int N_grid, T_Spectrum spec, double dir=0.) {
  double eta = 0.0;
  double ampl = 0.0;
  double a = 0.0;
  int index = 0;
  double dkx = spec.kx[1] - spec.kx[0];
  double dky = spec.ky[1] - spec.ky[0];
  int N_mode = spec.N_mode;
  int Ntmode = N_mode*2 + 1;
  //double kx2, ky2;
  for (int i = 0; i < Ntmode; i++) {
    for (int j = 0; j < Ntmode; j++) {
      index = i*Ntmode + j;
      if (spec.F_kxky[index] * dkx * dky < 0.){
        fprintf(stderr,"i=%d j=%d %f \n", i,j, spec.F_kxky[index] * dkx * dky);
      }
      ampl = sqrt(2. * spec.F_kxky[index] * dkx * dky);

      //adding a direction: change of coordinates
      // kx2 = spec.kx[i]*cos(dir) - spec.ky[j]*sin(dir);
      // ky2 = spec.kx[i]*sin(dir) + spec.ky[j]*cos(dir);
      // a = ( kx2*x + ky2*y + spec.phase[index]);



      // a = ( (spec.kx[i]*cos(dir)-spec.ky[j]*sin(dir))*x + 
      //     (spec.kx[i]*sin(dir)+spec.ky[j]*cos(dir))*y +
      //      spec.phase[index]);

      a = (spec.kx[i]*x + spec.ky[j]*y + spec.phase[index]);

      eta += ampl*cos(a);
    }
  }
  //printf("wave(): ampl %f, a %f, eta %f\n", ampl, a, eta);
  return eta;
}

// Velocities following the linear wave theory
double u_x(double x, double y, double z, int N_grid, T_Spectrum spec, double dir=0.) 
{
  int index = 0;
  double u_x = 0.0;
  double ampl = 0.0;
  double a = 0.0;
  double z_actual = 0.0;
  double kmod = 0.0;
  double theta = 0.0;
  double dkx = spec.kx[1] - spec.kx[0];
  double dky = spec.ky[1] - spec.ky[0];
  int N_mode = spec.N_mode;
  //double kx2, ky2;
  int Ntmode = N_mode*2 + 1;
  for (int i = 0; i < Ntmode; i++) {
    for (int j = 0; j < Ntmode; j++) {
      index = i*Ntmode + j;
      ampl = sqrt(2. * spec.F_kxky[index] * dkx * dky);
      z_actual = (z < ampl ? (z) : ampl);

      //adding a direction: change of coordinates
      //kx2 = spec.kx[i]*cos(dir) - spec.ky[j]*sin(dir);
      //ky2 = spec.kx[i]*sin(dir) + spec.ky[j]*cos(dir);

      //kmod = sqrt(sq(kx2) + sq(ky2));
      //theta = atan(ky2 / kx2);
      // a = (spec.kx[i] * x + spec.ky[j] * y + spec.phase[index]);
      //a = ( kx2*x + ky2*y + spec.phase[index]);

      kmod = sqrt(sq(spec.kx[i]) + sq(spec.ky[j]));
      //theta = atan(spec.ky[j] / spec.kx[i]);
      theta = atan2(spec.ky[j], spec.kx[i]);
      a = ( spec.kx[i]*x + spec.ky[j]*y + spec.phase[index]);

      u_x +=
          sqrt(g_ * kmod) * ampl * exp(kmod * z_actual) * cos(a) * cos(theta);
    }
  }
  return u_x;
}
double u_y(double x, double y, double z, int N_grid, T_Spectrum spec, double dir=0.) 
{
  int index = 0;
  double u_y = 0.0;
  double ampl = 0.0;
  double a = 0.0;
  double z_actual = 0.0;
  double kmod = 0.0;
  double theta = 0.0;
  double dkx = spec.kx[1] - spec.kx[0];
  double dky = spec.ky[1] - spec.ky[0];
  int N_mode = spec.N_mode;
  int Ntmode = 2*N_mode + 1;
  //double kx2, ky2;
  for (int i = 0; i < Ntmode; i++) {
    for (int j = 0; j < Ntmode; j++) {
      index = i*Ntmode + j;
      ampl = sqrt(2. * spec.F_kxky[index] * dkx * dky);
      z_actual = (z < ampl ? (z) : ampl);
      // kmod = sqrt(sq(spec.kx[i]) + sq(spec.ky[j]));
      // theta = atan(spec.ky[j] / spec.kx[i]);
      // //a = (spec.kx[i] * x + spec.ky[j] * y + spec.phase[index]);
      // a = ( (spec.kx[i]*cos(dir)-spec.ky[j]*sin(dir))*x + 
      //     (spec.kx[i]*sin(dir)+spec.ky[j]*cos(dir))*y +
      //      spec.phase[index]);


      //adding a direction: change of coordinates
      // kx2 = spec.kx[i]*cos(dir) - spec.ky[j]*sin(dir);
      // ky2 = spec.kx[i]*sin(dir) + spec.ky[j]*cos(dir);
      //
      // kmod = sqrt(sq(kx2) + sq(ky2));
      // theta = atan(ky2 / kx2);
      // a = (spec.kx[i] * x + spec.ky[j] * y + spec.phase[index]);
      //a = ( kx2*x + ky2*y + spec.phase[index]);

      kmod = sqrt(sq(spec.kx[i]) + sq(spec.ky[j]));
      //theta = atan(spec.ky[j] / spec.kx[i]);
      theta = atan2(spec.ky[j], spec.kx[i]);
      a = ( spec.kx[i]*x + spec.ky[j]*y + spec.phase[index]);

      u_y +=
          sqrt(g_ * kmod) * ampl * exp(kmod * z_actual) * cos(a) * sin(theta);
    }
  }
  return u_y;
}

double u_z(double x, double y, double z, int N_grid, T_Spectrum spec, double dir=0.) 
{
  int index = 0;
  double u_z = 0.0;
  double ampl = 0.0;
  double a = 0.0;
  double z_actual = 0.0;
  double kmod = 0.0;
  double dkx = spec.kx[1] - spec.kx[0];
  double dky = spec.ky[1] - spec.ky[0];
  int N_mode = spec.N_mode;
  int Ntmode = 2*N_mode + 1;
  //double kx2, ky2;
  for (int i = 0; i < Ntmode; i++) {
    for (int j = 0; j < Ntmode; j++) {
      index = i*Ntmode + j;
      ampl = sqrt(2. * spec.F_kxky[index] * dkx * dky);
      z_actual = (z < ampl ? (z) : ampl);
      // kmod = sqrt(sq(spec.kx[i]) + sq(spec.ky[j]));
      // //a = (spec.kx[i] * x + spec.ky[j] * y + spec.phase[index]);
      // a = ( (spec.kx[i]*cos(dir)-spec.ky[j]*sin(dir))*x + 
      //     (spec.kx[i]*sin(dir)+spec.ky[j]*cos(dir))*y +
      //      spec.phase[index]);


      //adding a direction: change of coordinates
      // kx2 = spec.kx[i]*cos(dir) - spec.ky[j]*sin(dir);
      // ky2 = spec.kx[i]*sin(dir) + spec.ky[j]*cos(dir);
      //
      // kmod = sqrt(sq(kx2) + sq(ky2));
      // // a = (spec.kx[i] * x + spec.ky[j] * y + spec.phase[index]);
      // a = ( kx2*x + ky2*y + spec.phase[index]);
  

      kmod = sqrt(sq(spec.kx[i]) + sq(spec.ky[j]));
      a = ( spec.kx[i]*x + spec.ky[j]*y + spec.phase[index]);

      u_z += sqrt(g_ * kmod) * ampl * exp(kmod * z_actual) * sin(a);
    }
  }
  return u_z;
}






