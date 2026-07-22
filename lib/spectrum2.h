#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "hugoj/lib/interpolate.h"
#pragma autolink -lgsl -lgslcblas
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#define PI 3.14159265358979323846

/** ## Some definitions and usfull functions */

typedef struct {
  int N_mode;
  double *kx;
  double *ky;
  double *F_kxky;
  double *phase;
  double *omega;
} T_Spectrum;

void free_spectrum(T_Spectrum s) {
    free(s.kx);
    free(s.ky);
    free(s.F_kxky);
    free(s.phase);
    free(s.omega);
}

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

/** ## Some common spectra */
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

/** ## Inverse FFT 
 see the original version from Andrés [here](https://basilisk.fr/sandbox/acastillo/input_fields/initial_conditions_dimonte_fft2.h)
*/
 void fft2D(double *data, int n0, int n1){  
  
  // Inverse FFT along rows 
  for (int i = 0; i < n0; ++i){
    gsl_fft_complex_radix2_backward(data + 2 * i * n1, 1, n1);
  }

  // Inverse FFT along columns
  double *column = malloc(2 * n0 * sizeof(double));
  for (int j = 0; j < n1; ++j){
    for (int i = 0; i < n0; ++i){
      REAL(column,i) = REAL(data, i*n1 + j);
      IMAG(column,i) = IMAG(data, i*n1 + j);
    }
    gsl_fft_complex_radix2_backward(column, 1, n0);
    for (int i = 0; i < n0; ++i)
    {
      REAL(data, i*n1 + j) = REAL(column,i);
      IMAG(data, i*n1 + j) = IMAG(column,i);
    }
  }
  free(column);
}


