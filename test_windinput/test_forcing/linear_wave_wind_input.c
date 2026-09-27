/**
# Wave growth: from linear to breaking

## Analytical solution

A linear surface gravity wave will decay in a viscous fluid. The rate of this
decay is $E(t)=E_0 e^{-4\nu k^2 t}$ (Lamb 1932)

## References

Horace Lamb, Hydrodynamics (6th ed., 1932), Chapter XI, Article 348, "Effect of
Viscosity on Water-Waves," pp. 623–625

et Article 349 pour une dérivation plus sérieuse


*/

// Faire le cas test !
// - propre
// - avec dimensions
// - simplifier

#include "grid/multigrid1D.h"
double p0 = 0.00625; // Pa
double rho = 1000;
#if HOLD_FORCING || EXACT_FORCING
// wind pressure param
#define wind_pressure(eta,i)  (p0/rho)*(eta[i+1] - eta[i-1])/(2*Delta)
// total barotropic pressure
#define p_baro(eta,i) (-G*eta[i] - wind_pressure(eta,i))

#define a_baro(eta,i) (gmetric(i)*(p_baro(eta,i)-p_baro(eta,i-1))/Delta) 

#endif // HOLD_FORCING || EXACT_FORCING

#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "bderembl/libs/netcdf_bas.h"

double k_ = 2.*pi, h_ = 1., g_ = 1.0, ak = 0.01; 
double RE = 20000.;


#define T0  (2.*pi/sqrt(g_*k_)) // = sqrt(2PI) = 2.5s
double lam;
double NT0 = 50.;
double etam_i = 0.;
double etavar_i = 0.;
double etavar_previous, etavar_current;
double cp, omega;
double relax_dt;

int main()
{
  L0=1.;
  origin (-L0/2.);
  periodic (right);
  N = 64;
  nl = 64; // 60
  G = g_;
  lam = 2*pi/k_*L0;
  
  omega = sqrt(g_*k_);
  cp = omega/k_;
  nu = cp*lam/RE; 
  theta_H=0.5065; // scheme conserve energy when theta_H = 0.5
  DT=0.02; // fixed DT to study spatial and temporal resolution separately
  NITERMIN=3; // Forces to do more cycles to avoid any influence of the poisson solver
  //TOLERANCE=1e-5;
  #if EXACT_FORCING
  p0 = 4*rho*nu*k_*cp; // Forcing to exactly balance viscous diss
  #endif
  relax_dt = 5*T0; ///(ak); //2*pi/(ak*omega);
  fprintf(stderr, "T0 = %f, lam=%f, g_=%f,omega=%f,nu=%g\n", T0, lam, g_, omega,
          nu);
  run();
}

void eta_stats (double *mean, double *variance)
{
  double sum = 0., var = 0.;
  foreach (reduction(+:sum))
    sum += eta[]; 
  *mean = sum/N;
  foreach (reduction(+:var))
    var += sq(eta[] - *mean);
  *variance = var/N; 
}

event init (i = 0)
{
  //geometric_beta (1./5., true);
  // default beta is 1/nl
  foreach() {
    zb[] = -h_;
    eta[] = ak/k_*cos(k_*x);
    double H = eta[] - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
      #if 1
      z += h_/nl/2; 
      u.x[] = ak/k_*sqrt(g_*k_)*exp(k_*z)*cos(k_*x); // 
      w[] = ak/k_*sqrt(g_*k_)*exp(k_*z)*sin(k_*x);  // 
      z += h_/nl/2; 
      #else
      z +=  h[]/2.;
      u.x[] = ak/k_*sqrt(g_*k_)*exp(k_*z)*cos(k_*x); // 
      w[] = ak/k_*sqrt(g_*k_)*exp(k_*z)*sin(k_*x);  // 
      z += h[]/2.;
      #endif
    }

  }
  // Compute initial wave energy
  eta_stats (&etam_i, &etavar_i);
  fprintf (stderr, "INITIAL eta mean = %.10f, variance = %.10f\n", etam_i, etavar_i);
  fprintf (stderr, "initial p0 = %f\n", p0);
  
  create_nc({zb, eta, h, u.x, w}, "out.nc");
}

/* compute p0 from eta in a 'face_fields' event so that it uses eta from previous timestep,
   before applying the forcing using the macro 'a_baro'
*/

#if HOLD_FORCING
// event face_fields (i++, last)
// {
//   double etam, etavar;
//   eta_stats (&etam, &etavar);
//   p0 = rho*G*(etavar_i-etavar)/(pi*sqrt(etavar_i));
//   fprintf (stderr, "i=%d,eta mean = %g, variance = %g  p0=%g\n", i, etam, etavar, p0);
// }

event update_p0 (i++)
{
  double etam;
  eta_stats (&etam, &etavar_current);
  double dE = G*(etavar_i - etavar_current);
  #if 0
  double eta_rms = sqrt(etavar_i);
  p0 = rho*dE/(pi*eta_rms);
  #else
  double sum = 0.;
  foreach() {
    double integrated_transport=0.;
    double etaxx = (eta[1] + eta[-1]- 2*eta[])/sq(Delta);
    foreach_layer()
      integrated_transport += h[] * u.x[]; // h * detadx * u
    sum += integrated_transport*etaxx*dv();
  }
  p0 = -rho*dE/(relax_dt*sum);
  fprintf(stderr,
        "i=%d dE=%g sum=%g dt=%g p0=%g predicted=%g\n",
        i, dE, sum, relax_dt, p0,
        -p0*sum*relax_dt/rho);
  #endif
//   fprintf(stderr,
//           "i=%d t=%g var=%g dE=%g p0=%g\n",
//           i, t, etavar_current, dE, p0);
// 
}

#endif 


event viscous_term (i++) {
  // vertical diffusion is done in diffusion.h (u.x) and nh.h (w)
  horizontal_diffusion ({u.x, w}, nu, dt);
}


// #if FORCING
// event face_fields (i++, last)
// {
//   /* compute p0 from eta here */
//   ...
// }
// #endif


 // #if FORCING
 // #warning "FORCING accel branch active"
 // event acceleration(i++, last){
 //   foreach_face(x) {
 //     double detadx = 0.;
 //     detadx = (eta[] - eta[-1])/Delta;
 //     ha.x[0,0,nl-1] += hf.x[]*(p0 * detadx * detadx / sqrt(1+detadx*detadx)); //   
 //   }
 // }
 // #endif


event logfile (i++; t <= NT0*T0) // target: at least 300*T0
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    foreach_layer() {
      double norm2 = sq(w[]);
      foreach_dimension()
	norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
    }
    gpe += sq(eta[])*dv();
  }
  fprintf (stdout, "%g %g %g\n", t/T0, ke/2., g_*gpe/2.);
}



event writenc (i+=10; t <= NT0*T0) // target: at least 300*T0
{
  write_nc();
}



/**
~~~pythonplot Wave energy
import numpy as np
import matplotlib.pyplot as plt

g = 1.0
L0 = 1.
ak = 0.05
k = 2*np.pi/L0
lam = 2*np.pi/k
Re = 40000
c = np.sqrt(g*k)/k
NT0 = 100

nu = c*lam/Re
T0 = (2*np.pi/np.sqrt(g*k))

def E_linwave(E0,nu,ak,k,t):
  print('\nWave decay for linear wave (theory)')
  print(f'nu={nu},ak={ak},k={k/np.pi}pi\n')
  return E0*np.exp(-4*nu*k**2*t)

data = np.loadtxt("out",skiprows=1)
#data_ref = np.loadtxt('../no_forcing/out',skiprows=1, usecols=np.arange(8700))
time = data[:,0]
ke  = data[:,1]
gpe = data[:,2]
E = ke + gpe
#timeref = data_ref[:,0]
#Eref = data_ref[:,1]+data_ref[:,2]
E0 = E[0]
Eth = E_linwave(E0,nu,ak,k,time*T0)
#E0 = 1 
fig, ax = plt.subplots(figsize=(5, 5))
#ax.plot(timeref,Eref/E0,c='k',ls='--',label='no forcing')
ax.plot(time,2*ke/E0,color='b', label='2*ke')
ax.plot(time,2*gpe/E0,color='g', label='2*gpe')
ax.plot(time,E/E0, color='k', label='E')
ax.plot(time, Eth/E0, color='r', label=r'$E(t)=E_0 e^{-4 \nu k^2 t}$')
ax.set_xlabel("t/T0")
ax.set_ylabel("E/E0")
#ax.set_xlim([0,NT0])
ax.set_ylim([0,1.5])
ax.legend(loc="lower left")
plt.tight_layout()
plt.savefig("energy.png", dpi=150)
plt.show()

# plot [:2] "out" using 1:($2+$3) w l
~~~
*/
