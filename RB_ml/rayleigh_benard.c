/*
 # Rayleigh-Benard convection using the multilayer
 
 This code is used to test the multilayer code with a heat flux imposed at
 boundaries.
 
*/
//#include "grid/multigrid.h"
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/diffusion.h"
#include "layered/perfs.h"
#include "bderembl/libs/extra.h"      // parameters from namlist
#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files

#define g_ 9.81


// TODO: add namelist file
/*
DEFAULT PARAMETERS

*/
char namlist[80] = "namelist.toml";    // file name of namlist
char file_out[20] = "out.nc";         // file name of output
// -> Initial conditions
double T0= 20.;             // [K] Initial T
// -> Forcing
double qtop = -100.;           // [W.m-2] Heat flux at top
double qbot = 100.;           // [W.m-2] Heat flux at bottom
// -> Domain definition
int N_grid = 5 [];          // 2^N_grid : number of x and y gridpoints
double L = 200.0 [1];       // domain size
int N_layer = 2 [];         // number of layers
double h0 = 1.0 [1];        // depth of water
// -> Runtime parameters
double tend = 2.0;          // end time of simulation
// -> saving outputs
double dtout = 2.0;         // dt for output in netcdf
// -> physical properties
double nu0 = 0.00025;       // Viscosity for vertical diffusion
double thetaH = 0.5;        // theta_h for dumping fast barotropic modes
// -> stratification related 
double rho0 = 1025.;     // [kg.m-3] reference density
double cp = 4.2e3;       // [J.kg-1.K-1] heat capacity water
double betaT = 2e-4;     // Linear equation of state: drho = betaT*(T0-T) (Vallis 2.4)
double D=1.5e-5;

/*
I should non dimensionalize the problem. 
A=L/H : aspect ratio
Pr=nu/alpha : Prandtl number (momentum viscosity vs thermal viscosity)
Ra=g.BetaT.DeltaT.H^3/(nu.alpha) : Rayleigh number (buoyancy vs viscous forces)
Nu=<w'T'>/ ? : Nusselt number (heat transfer normalized)
*/


#define drho(T) (betaT*(T0-T))
#include "layered/dr.h"

int main(int argc, char *argv[])  
{
 
  // Building a 'params' array with all parameters from the namlist
  params = array_new();
  add_param("N_grid", &N_grid, "int");
  add_param("L", &L, "double");
  add_param("N_layer", &N_layer, "int");
  add_param("h0", &h0, "double");
  add_param("tend", &tend, "double");
  add_param("nu0", &nu0, "double");
  add_param("thetaH", &thetaH, "double");
  add_param("dtout", &dtout, "double");
  add_param("T0", &T0, "double");
  add_param("rho0", &rho0, "double");
  add_param("cp", &cp, "double");
  add_param("betaT", &betaT, "double");
  add_param("qtop", &qtop, "double");
  add_param("qbot", &qbot, "double");

  // Search for the configuration file with a given path or read params.in
  if (argc == 2)
    strcpy(file_param, argv[1]);
  else
    strcpy(file_param, namlist);
  read_params(file_param);
  
  // Settings solver values from namlist values
  L0 = L;
  nu = nu0;
  N = 1 << N_grid; // 1*2^N_grid
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  CFL_H = 1; // Smaller time step


  // Boundary condition
  origin (-L0/2., -L0/2.);
  periodic (right);
  
  fprintf (stderr, "Read in parameters!\n");
  run();
}



event init(i =  0) {
  geometric_beta (0., true); // Varying layer thickness
  // step 1: set eta and h
  foreach() {
    zb[] = -h0;
    eta[] = 0.;
    double H = - zb[];
    foreach_layer() {
      h[] = H/nl;
    } 
  }
  
  // step 2: remap
  vertical_remapping (h, tracers);
  // step 3: set currents
  foreach() {
    foreach_layer() {
      u.x[] = 0.;
      u.y[] = 0.;
      w[] = 0.;
      T[] = T0 -0.01 + (rand() / (RAND_MAX+1.0) * 0.02);
    }
  }
  create_nc((scalar *){zb, h, u, w, eta, T}, file_out);
  fprintf (stderr,"Done initialization!\n");

}
// This implements a heat flux at the surface and bottom
event viscous_term (i++)
{
  foreach() {
    vertical_diffusion (point, h, T, dt, D, qt/(D*rho0*cp), 0., 0.);
    // TODO: add bottom flux (I need to convert Navier slip to neumann)
  }
}

// Writing a 4D netcdf file
event output(t = 0.; t<= tend+1e-10; t+=dtout){
  write_nc();
}


