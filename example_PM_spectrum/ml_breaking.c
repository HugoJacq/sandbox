/**
  # Wave breaking field decay

No stratification, wave decay (no forcing)
*/
#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#define g_ 9.81

/*

TO DO : 
- add dimensions so that they can be saved in the netcdf

- use correct parameters

*/


/*
 ## PARAMETERS

Dimensions : [Length, Time]
*/
// -> Initial conditions
double P = 0.2 [1, -1];     // energy level (estimated so that kpHs is reasonable)
int coeff_kpL0 = 10 [];     // kpL0 = coeff_kpL0 * pi
int N_mode = 32 [];         // Number of modes in wavenumber space
int N_power = 5 [];         // directional spreading coeff
int F_shape = 0 [];         // shape of the initial spectrum
// -> Domain definition
int N_grid = 5 [];       // 2^N_grid : number of x and y gridpoints
double L = 200.0 [1];       // domain size
int N_layer = 2 [];         // number of layers
double kp = PI*10/200.0 [-1];// peak wave number
double h0 = 1.0 [1];        // depth of water
// -> Runtime parameters
int restart = 0;            // 1: restart, 0: no restart
double tend = 2.0;          // end time of simulation
// -> saving outputs
double dtout = 2.0;        // dt for output in netcdf
int pad = 4;                // number of 0-padding for ouput files
int nout = 1;               // number of the outfile
char fileout[100];          // name of outfile
// -> physical properties
double nu0 = 0.;            // Viscosity for vertical diffusion
int RANDOM = 2;             // For random number generator
double thetaH = 0.5;        // theta_h for dumping fast barotropic modes


// diag
double *u_profile;
double dt_mean = 5;  // s
static FILE * fp;




/*
## Initialisation

Let us define some functions
*/
#define PI 3.141592
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

double randInRange(double min, double max) {
  return min + (rand() / (RAND_MAX+1.0) * (max - min));
}

/*
 The initial wave field spectrum has the Pierson-Moscowitz shape
 */
double spectrum_PM(double P, double kp, double kmod) {
  // Note: P here is in fact P/sqrt(g) in the PM spectrum equation
  return P * pow(kmod, -2.5) * exp(-1.25 * pow(kp / kmod, 2.0));
}
/*
 We use an bilinear interpolation function
 */
int find_index_closest(double arr[], int length, double target) {
  /*
   Find the index of the closest element in 'arr'.
   */
    int res = 0;
    int lo = 0, hi = length - 1;

    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        // Update res if mid is closer to target
        if (fabs(arr[mid] - target) < fabs(arr[res] - target)) {
            res = mid;
        }

        if (arr[mid] == target) {
            return mid;
        }
        else if (arr[mid] < target) {
            lo = mid + 1;
        }
        else {
            hi = mid - 1;
        } 
    }
    return res;
}

void find_ibounds(double arr[], int length, double target, int* low, int* high){
  /* This function finds the bounds in 'arr' that closely match 'target'
   
    we assume that 'target' is inside the range of 'arr'
   
  */

  int closest = find_index_closest(arr, length, target);

  if (closest == length) {
    *low = closest-1;
    *high = closest;
    return;
  }

  if (arr[closest]-target < 0) {
    *low = closest;
    *high = closest+1;
  }
  else {
    *low = closest-1;
    *high = closest;
  }
}


double interp_lin(double x[], double y[], int Nx, int Ny, double xi, double yi, double F[]) {
  /* Interpolate linearly F on grid (x,y) at position xi,yi 
   *
   * Using eq 98 of https://pages.hmc.edu/ruye/MachineLearning/lectures/ch7/node7.html
   * */
  
  int j0, j1, i0, i1;
  double Fi;
  find_ibounds(x, Nx, xi, &i0, &i1);
  find_ibounds(y, Ny, yi, &j0, &j1);
  
  double dx = xi - x[i0];
  double dy = yi - y[j0];
  double deltaX = x[i1] - x[i0];
  double deltaY = y[j1] - y[j0];
  double f11 = F[i1*Nx+j1];
  double f10 = F[i1*Nx+j0];
  double f01 = F[i0*Nx+j1];
  double f00 = F[i0*Nx+j0];

  Fi = ( dx/deltaX * dy/deltaY * f11 +
      dy/deltaY * (1-dx/deltaX) * f01 +
      dx/deltaX * (1-dy/deltaY) * f10 +
      (1-dx/deltaX-dy/deltaY+ dx/deltaX*dy/deltaY) * f00 );
  
  return Fi;
}

/* The function that generate the F(kx,ky) spectrum from an omnidirectional
 spectrum
 */
T_Spectrum spectrum_gen_linear(int N_mode, int N_power, double L, double P,
                               double kp) 
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
  */

  int N_kmod = 64; // Uniform grid in kmod and ktheta, can be finer than N_mode
  int N_theta = 64;
  double thetam = 0.0; // midline direction
  double theta[N_theta];
  double dtheta;
  double Dtheta[N_theta];               // for directional spectrum
  double sum = 0.0;                     // for Dtheta normalization
  double F_kmodtheta[N_kmod * N_theta] ; // directional spectrum
  double kmod[N_kmod];
  double F_kmod[N_kmod];

  T_Spectrum spectrum;
  spectrum.N_mode = N_mode;
  spectrum.kx = (double *)malloc(N_mode * sizeof(double));
  spectrum.ky = (double *)malloc((N_mode + 1) * sizeof(double));  
  spectrum.F_kxky = (double *)malloc(N_mode * (N_mode + 1) * sizeof(double));
  spectrum.phase = (double *)malloc(N_mode * (N_mode + 1) * sizeof(double));
  spectrum.omega = (double *)malloc(N_mode * (N_mode + 1) * sizeof(double));

  // building kmod
  for (int i = 0; i < N_kmod; ++i) {
    //spectrum.kmod[i] =
    kmod[i] =
        2 * PI / L + 1.0 * i / (N_kmod - 1) * (1.41 * 100 * 2 - 2) * PI / L;
  }
  // build Dtheta
  for (int i = 0; i < N_theta; ++i) {
    theta[i] = -0.5 * PI + 1.0 * i / (N_theta - 1) * PI;
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
  // TO DO: add more shapes ?
  for (int i = 0; i < N_kmod; ++i) {
    F_kmod[i] = spectrum_PM(P, kp, kmod[i]);
  }
  for (int ik = 0; ik < N_kmod; ++ik) {
    for (int itt = 0; itt < N_theta; ++itt) {
      F_kmodtheta[ik * N_kmod + itt] =
          F_kmod[ik] * Dtheta[itt] / kmod[ik]; // Notice the normalize by kmod ! 
    }
  }

  // Uniform grid in kx ky
  for (int i = 0; i < N_mode; ++i) {
    spectrum.kx[i] = 2 * PI / L * (i + 1);
  }
  for (int i = 0; i < N_mode + 1; ++i) {
    spectrum.ky[i] = 2 * PI / L * (i - N_mode / 2);
  }

  // interp F_kmodtheta on kx,ky grid
  double rho, phi;
  
  for (int ix = 0; ix < N_mode; ++ix) {
    for (int iy = 0; iy < (N_mode + 1); ++iy) {
      // first we get polar coords
      cart2pol(spectrum.kx[ix], spectrum.ky[iy], &rho, &phi);
      // then interp at these coords
      spectrum.F_kxky[ix * N_mode + iy] = interp_lin(
        kmod, theta, N_kmod, N_theta, rho, phi, F_kmodtheta);

    }
  }

  // Random phase
  int RANDOM=2;
  srand(RANDOM); // We can seed it differently for different runs
  int index = 0;
  double k = 0;
  //double randnum;
  for (int i=0; i<N_mode; i++) {
    for (int j=0; j<N_mode+1; j++) {
      index = j*N_mode + i;
      k = sqrt(sq(spectrum.kx[i]) + sq(spectrum.ky[j]));
      spectrum.omega[index] = sqrt(g_*k);
      spectrum.phase[index] = randInRange (0, 2.*PI);
      //printf(" random number = %f\n", spectrum.phase[index]);
    }
  }
  return spectrum;
}





/* ## Main program */
int main(int argc, char *argv[])  
{
  
  kp = PI * coeff_kpL0 / L; // kpL=coeff x pi peak wavelength
  
  L0 = L;
  nu = nu0;
  N = 1 << N_grid; // 1*2^N_grid
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  CFL_H = 1; 

  // Boundary condition
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  
  // allocate diag
  u_profile = (double*)calloc(nl, sizeof(double));
  fp  = fopen("u_profile.dat","w"); // reset file
  fclose(fp);
   
  fprintf (stderr, "Read in parameters!\n");
  run();
}

/* Initial conditions */
event init(i =  0) {
  geometric_beta (1./3., true); // Varying layer thickness
  if (!restore ("restart")) {
    // Generate linealy spaced kx, ky according to specified # of modes, and
    //  then interpolated on a cartesian grid to get F(kx,ky).
    T_Spectrum spectrum;
    spectrum = read_spectrum(N_mode); // TO DO:

    // step 1: set eta and h
    foreach() {
      zb[] = -h0;
      eta[] = wave(x, y, N_grid, spectrum);
      double H = wave(x, y, N_grid, spectrum) - zb[];
      foreach_layer() {
        h[] = H*beta[point.l]/nl;
      } 
    }
    // step 2: set currents
    foreach() {
      double z = zb[];
      foreach_layer() {
        z += h[]/2.;
        u.x[] = u_x(x, y, z, N_grid, spectrum);
        u.y[] = u_y(x, y, z, N_grid, spectrum);
        w[] = u_z(x, y, z, N_grid, spectrum);
        z += h[]/2.;
      }
    }
  }
  else {
    // We limit the first time step after the restart
    geometric_beta (1./3., true); // when restarting, remember to specify the grid mapping method
    dtmax = 0.01;
    dt = dtnext (dtmax);
  }

  fprintf (stderr,"Done initialization!\n");
}


// event film surface


// event images (t += dtout) {
//   scalar img[];
//   foreach()
//     img[] = u.x[0,0,l];
//   static FILE * fp = fopen ("grid.ppm", "w");
//
//   output_ppm (img, fp, min = -2, max = 2);
// }
//
// event dumpini (t=0.){
//   dump();
// }

/* Diagnostic of the layer averaged zonal velocity */
// This event compute layer average of u.x
event compute_horizontal_avg (i++; t<=tend+1e-10){
 event compute_horizontal_avg (t+=dt_mean; t<=tend+1e-10){
   foreach(reduction(+:u_profile[:nl]))
     foreach_layer(){
       u_profile[point.l] += u.x[] / (N*N) * dt / dt_mean;

     }
 }

 // This event writes to a file the layer average
 event write_diag(t=0., t+=dt_mean){
     // main worker is writing the file
     if (pid()==0) {
       fp  = fopen("u_profile.dat","a");
       if (fp == NULL){
         fprintf(stderr, "Error opening file u_profile.dat");
         return 2;
       }
       for (int i=0; i<nl; ++i) {
         fprintf (fp, "%f %d %g\n", t, i, u_profile[i]);
       }
       fprintf(fp,"\n");
       fclose(fp);
     }
     // Reset the profile for all workers
     for (int i=0; i<nl; ++i) {
       u_profile[i] = 0.0;
     }


// event image (t = end) {
//   clear();
//   static FILE * fp = fopen ("image.ppm", "w");
//   output_ppm(eta, fp, min=-0.1, max=0.1);
//
// }

// event snapshot (t = end)
// {
//   clear();
	
  // 3/4 view
  // view (quat = {0.567, 0.137, 0.196, 0.789},
  //     fov = 30, near = 0.01, far = 1000,
  //     tx = 0.048, ty = -0.001, tz = -2.096,
  //     width = 1398, height = 803);
  //
  //
  // char s[80];
  // sprintf (s, "t = %.2f", t);
  // draw_string (s, size = 80);
  // for (double x = -1; x <= 1; x++)
  //   translate (x) {
  //     squares ("eta", linear = true, z = "eta*10", min = -1.0, max = 1.0 , map=gray);
  //   }
  // box ();
  //
  //
  //
  //
  //
  // colorbar(map=gray, label="eta (m)", min=-1.0,max=1.0, pos={-0.95,-0.5},
  //          levels=10);
  // save ("snap_side.png");
  //
  //
  // // top view
  // clear();
  // view (quat = {0.000, 0.000, 0.000, 1.000},
  //     fov = 30, near = 0.01, far = 1000,
  //     tx = 0.000, ty = 0.000, tz = -2.505,
  //     width = 1398, height = 803);
  //
  //
  // sprintf (s, "t = %.2f", t);
  // draw_string (s, size = 80);
  // for (double x = -1; x <= 1; x++)
  //   translate (x) {
  //     squares ("eta", min = -1.0, max = 1.0 , map=gray);
  //   }
  // box ();
  //
  // colorbar(map=gray, label="eta (m)", min=-1.0,max=1.0, pos={-0.95,-0.5},
  //          levels=10);
  // save ("snap_top.png");


//}

// Clean for my diag of layer avg
event cleanup(t=end){
  free(u_profile);
}

/**
Results: plots
~~~gnuplot 
# Plot the heatmap
set pm3d map
set view map
set xlabel "Time (s)"
set ylabel "Layer"
set cblabel "u.x (m/s)"
set yrange [0:29]
set xrange [0:200]
set terminal pngcairo size 800,600 enhanced font 'Verdana,12'
set output 'u_profile.png'
set size 0.9, 0.9
splot "u_profile.dat" using 1:2:3 with pm3d
unset output
~~~
**/

