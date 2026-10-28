#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
# basic spectral tools for data in square domains
# based on https://github.com/adekunleoajayi/powerspec
#
# from https://github.com/bderembl/geostrokit/blob/master/geostrokit/fftlib.py

def azimuthal_integral(spec_2D, Delta, all_kr=False):
    '''
    Compute the azimuthal integral of a 2D spectra

    Parameters
    ----------

    spec_2D: array [(nz,) ny,nx] 2d spectra
    Delta: float, grid step
    all_kr: bool. if False, only get radial wave number
    in the inner circle, if True, gets all wave numbers. 
    Default is False. You should only use True to check Parseval

    Returns
    -------

    kr: array [nkr]  radial wave number
    spec_1D: array [(nz,) nkr] 1d spectra
    '''
    
    #print('import is ok')

    nd = spec_2D.ndim
    if nd == 2:
        spec_2D = spec_2D[None,...]
    nl,N,naux = spec_2D.shape

    k,l,K,kr = get_wavenumber(N,Delta, all_kr)
    dk = kr[1] - kr[0]
    spec_1D = np.zeros((nl,len(kr)))
    for i in range(kr.size):
        kfilt =  (K>=kr[i] - 0.5*dk) & (K<kr[i] + 0.5*dk)
        #    the azimuthal integral is the average value*2*pi*k
        # but to get the same value of the integral for the 1d spetrum
        # and the 2d spectrum, it is better to just sum the cells*dk

        spec_1D[:,i] = (spec_2D[:,kfilt].sum(axis=-1))*dk #/(2*np.pi) #*kr[i]*2*np.pi/Nbin

        # Nbin = kfilt.sum()
        # spec_1D[:,i] = (spec_2D[:,kfilt].sum(axis=-1))*kr[i]/(Nbin)*2*np.pi
        
        # the loop is missing the value at K=0:
        # add it with all_kr option
    return kr, spec_1D.squeeze()


def get_wavenumber(N, Delta, all_kr=False):
    '''
    Compute wavenumber and radial wavenumber 
    wave number units are in cycle per unit length

    Parameters
    ----------

    N : int
    Delta: float, grid step
    all_kr: bool. if False, only get radial wave number
        in the inner circle, if True, gets all wave numbers. 
        Default is False. You should only use True to check Parseval

    Returns
    -------

    k: array [ny,nx] zonal wave numbers
    l: array [ny,nx] meridional wave numbers
    K: array [ny,nx]  (k^2 + l^2)
    kr: array [nkr]  radial wave number
    '''

    kx = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(N,Delta)) # two sided  
    k,l = np.meshgrid(kx,kx)
    K = np.sqrt(k**2 + l**2)
    if all_kr == False:
        kr = k[0,int(N/2)+1:]
    elif all_kr == True:
        kmax = K.max()
        dk = np.abs(kx[2]-kx[1])
        kr = dk*np.arange(0,int(kmax/dk)+2)
    return k,l,K,kr


def get_spec_2D(psi1, psi2, Delta, window=None):
  ''' 
  Compute the 2D cross power spectrum of psi1 and psi2
  If psi1=psi2, then this is simply the power spectrum

   normalization such that parseval is ok: 
   E = 1/V*(np.sum(psi**2)*Delta**2) = np.sum(psi**2)/N**2
     = np.sum(spec_2D)*dk**2
     ~ np.sum(spec_1D)*dk
  (the last equality is not exactly true because the 1d spectra only contains 
  the inner circle to the square k,l. Set all_kr=True in get_spec_1D)
  
  with V = (N*Delta)**2

  Parameters
  ----------

  psi1 : array [(nz,) ny,nx]
  psi2 : array [(nz,) ny,nx]
  Delta: grid step
  window: None or "hanning"

  Returns
  -------

  k: array [ny,nx] zonal wave numbers
  l: array [ny,nx] meridional wave numbers
  spec_2D: array [(nz,) ny,nx] 2d spectra
  '''
  
  N = psi1.shape[-1]
  k,l,K,kr = get_wavenumber(N,Delta)

  if window == 'hanning':
    win1d = np.hanning(N)
    window = win1d[None,:]*win1d[:,None]
  else:
    window = np.ones((N,N))

  psi1 = window*psi1
  psi2 = window*psi2

  psi1_hat = np.fft.fft2(psi1)
  psi2_hat = np.fft.fft2(psi2)
  spec_2D = (psi1_hat*psi2_hat.conj()).real*( Delta / (2*np.pi) )**2/N**2
  spec_2D = np.fft.fftshift(spec_2D, axes=(-1,-2))
  return k, l, spec_2D


def get_spec_1D(psi1, psi2, Delta, window=None, all_kr= False, averaging='radial'):
  '''
  Compute the 1D power spectrum of the data by azimuthal integral
  of the 2D spectra

  Parameters
  ----------

  psi1 : array [(nz,) ny,nx]
  psi2 : array [(nz,) ny,nx]
  Delta: grid step
  window: None or "hanning"
  all_kr: bool. if False, only get radial wave number
    in the inner circle, if True, gets all wave numbers. 
    Default is False. You should only use True to check Parseval
  averaging: 'radial' (default), or 'xy' 

  Returns
  -------

  if averaging == 'radial'

  kr: array [N] radial wave number
  spec_1D: array [(nz,) N] 1d spectra

  elif averaging == 'xy'

  k: array [N-1] zonal wave numbers
  l: array [N-1] meridional wave numbers
  spec_x: array [(nz,) N-1] 1d spectra in zonal direction
  spec_y: array [(nz,) N-1] 1d spectra in meridional direction

  '''

  k, l, spec_2D = get_spec_2D(psi1, psi2, Delta, window)

  if averaging == 'radial':
    kr, spec_1D = azimuthal_integral(spec_2D, Delta, all_kr)
    return kr, spec_1D

  elif averaging == 'xy':
    si = psi1.shape
    N = si[-1]
    if all_kr == False:
      kx = k[0,int(N/2)+1:]
      ky = l[int(N/2)+1:,0]
    else:
      kx = -k[0,:int(N/2)+1][::-1]
      ky = -l[:int(N/2)+1,0][::-1]
      


    dk = kx[1] - kx[0]
    dl = ky[1] - ky[0]
    spec_x = np.sum(spec_2D,axis=-2)*dl
    spec_y = np.sum(spec_2D,axis=-1)*dk

    if all_kr == False:
      spec_x = spec_x[...,int(N/2)+1:] + spec_x[...,1:int(N/2)][::-1]
      spec_y = spec_y[...,int(N/2)+1:] + spec_y[...,1:int(N/2)][::-1]
    else:
      spec_x[...,1:int(N/2)] += spec_x[...,int(N/2)+1:][...,::-1]
      spec_y[...,1:int(N/2)] += spec_y[...,int(N/2)+1:][...,::-1]
      spec_x = spec_x[...,:int(N/2)+1][...,::-1]
      spec_y = spec_y[...,:int(N/2)+1][...,::-1]


    return kx, ky, spec_x, spec_y




def get_spec_flux(psi1, psi2, Delta, window=None):
  '''
  Compute spectral flux

  Parameters
  ----------

  psi1 : array [(nz,) ny,nx]
  psi2 : array [(nz,) ny,nx]
  Delta: grid step
  window: None or "hanning"

  Returns
  -------

  kr: array [N] radial wave number
  flux: array [(nz,) N] spectral flux

  '''
  k, l, spec_2D = get_spec_2D(psi1, psi2, Delta, window)

  # kr,spec_1D = azimuthal_integral(spec_2D,Delta)
  # flux = -np.cumsum(spec_1D)*dk # integrate from low wavenumbers
  # flux = np.cumsum(spec_1D[::-1])[::-1]*dk # integrate from high wavenumbers

  nd = spec_2D.ndim
  if nd == 2:
    spec_2D = spec_2D[None,...]
  nl,N,naux = spec_2D.shape

  k,l,K,kr = get_wavenumber(N,Delta)
  dk = kr[1] - kr[0]

  flux = np.zeros((nl,len(kr)))
  for i in range(kr.size):
    kfilt =  (kr[i] <= K ) 
    flux[:,i] = (spec_2D[:,kfilt]).sum(axis=-1)*dk*dk
    # next line is probably faster. check it is the same result
    #flux[:,i] = (spec_2D*kfilt).sum(axis=(1,2))*dk*dk

  return kr, flux.squeeze()


def convolve2D(psi1, psi2, norm=False, psi_s=1, kr=0, k_target=0, dk=1, K=0):
  '''
  Convolution of psi1 and psi2 (assume these two fields are in spectral space)

  Parameters
  ----------

  psi1 : array [ny,nx] Shiffted Fourier coefs of psi1
  psi2 : array [ny,nx] Shiffted Fourier coefs of psi1
  norm: Bool if True, normalize the convolution by 1/N**2 (optional)
  psi_s : array [ny,nx] Shiffted Fourier coefs of psi_s (for scalar product)(optional)
  kr: array [N] radial wave number (optional)
  k_target: float (optional)
  dk: float (optional)
  K: wave number magnitude (optional)

  Returns
  -------

  out: array [ny,nx], convolution of psi1*psi2
  kmag: array [ny,nx], magnitude of k (1st dimension is K of psi2 second dimension is K of psi1)

  '''

  # only convol odd arrays
  si = psi1.shape
  N = si[-1]
  Nc = int((N-1)/2)

  if norm:
    cnorm = 1/N**2  # N or N-1??
  else:
    cnorm = 1

  flag_reshape = 0
  if N % 2 == 0:
    flag_reshape = 1
    psi1 = psi1[1:,1:]
    psi2 = psi2[1:,1:]
    if isinstance(K,np.ndarray):
      K = K[1:,1:]
    if isinstance(psi_s,np.ndarray):
      psi_s = psi_s[1:,1:]

  si = psi1.shape
  N = si[-1]
  Nc = int((N-1)/2)

  if not isinstance(psi_s,np.ndarray):
    psi_s = np.ones((N,N))

  out = np.zeros_like(psi1)

  if k_target:
    nk = len(kr)
    kmag = np.zeros((nk,nk))
    kfilt = (K>=k_target - 0.5*dk) & (K<k_target + 0.5*dk)
  else:
    nk = 1
    kmag = np.zeros((nk,nk))
    kfilt = np.zeros((N,N))

  psif = np.flip(np.flip(psi2,-1),-2)

  for j,i in np.ndindex((N,N)):

    # indices of the flipped matrix
    i1 = None if i >= Nc else Nc-i
    i2 = None if i <= Nc else Nc-i
    j1 = None if j >= Nc else Nc-j
    j2 = None if j <= Nc else Nc-j

    # indices of the non-flipped matrix
    k1 = None if i <= Nc else i-Nc
    k2 = None if i >= Nc else i-Nc
    l1 = None if j <= Nc else j-Nc
    l2 = None if j >= Nc else j-Nc
    #print("(",i,",",j,")",i1,i2,j1,j2, "--", k1,k2,l1,l2)
    out[j,i] = psi_s[j,i].conj()*np.sum(psif[j1:j2,i1:i2]*psi1[l1:l2,k1:k2])*cnorm

    if kfilt[j,i]:
      k_ind1 = np.round((K[l1:l2,k1:k2]/dk).flatten()).astype(int)
      k_ind2 = np.round((K[j1:j2,i1:i2]/dk).flatten()).astype(int)
      kmag[k_ind2,k_ind1] += (psi_s[j,i].conj()*
                              (psif[j1:j2,i1:i2]*psi1[l1:l2,k1:k2]*cnorm).flatten()).real

  if flag_reshape:
    out = np.pad(out,((1,0),(1,0)))

  if k_target == 0:
    return out
  else:
    return out, kmag


def directional_wave_spectrum(
    eta,
    u,
    v,
    dx,
    dy,
    g=9.81,
    z=None,
    w=None,
):
    """
    Estimate directional deep-water wave components from one snapshot.

    The space-time Fourier transform of elevation is, 
    $$\eta_t({\bm x_h}) = \int \int \hat{\eta}(\bm{k},\omega) e^{i(-\omega t + \bm{k}\cdot\bm{x}_h )} d\bm{k} d\omega$$

    where $\bm{x}_h = (x,y)$

    Then assuming linear deep-water waves, $\omega$ can take only two values: $\omega = \pm \sqrt{gk}$ where $k=\| \bm{k}\|$.
    Associated complex amplitudes are $a^{\pm}_t (\bm{k}) = \hat{\eta}(\bm{k},\omega = \pm \sqrt{gk}) e^{i(\mp\omega t)}$, and elevation becomes

    $$\eta_t({\bm x_h}) = \int (a^+_t(\bm k) + a^-_t(\bm k)) e^{i \bm{k}\cdot\bm{x}_h} d\bm{k} $$

    Now if we consider a snapshot at fixed time $t=t_0$, then the \textit{spatial} Fourier transform is 
    $$\eta_{t_0}({\bm x_h}) = \int \hat{\eta}_{t_0}(\bm k) e^{i \bm{k}\cdot\bm{x}_h} d\bm{k} $$

    and we have $\hat{\eta}_{t_0}(\bm k) = a^+_{t_0}(\bm k) + a^-_{t_0}(\bm k)$. 
    But we cannot compute $a^+$ and $a^-$ separately. The spectrum from one snapshot is then said to be "ambiguous". 
    For example, from $\eta_{t=t_0} = \cos (k x )$ it is not possible to know if this wave elevation corresponds to forward or backward propagating wave. 

    Additional information to disambiguate the spectrum can be provided by the horizontal velocity  spatial spectrum,

    $$\bm{u}_{t_0}({\bm x_h},z) = \int (\hat{u}_{t_0}(\bm k,z),\hat{v}_{t_0}(\bm k,z)) e^{i \bm{k}\cdot\bm{x}_h} d\bm{k} $$

    If we assume again linear deep-water waves, then the Fourier amplitude of velocity satisfies

    $$  (\hat{u},\hat{v}) = a^+ \sqrt{gk} e^{kz} \frac{\bm k}{k} + a^- (-\sqrt{gk}) e^{kz} \frac{\bm k}{k} $$

    Let us denote the projection of Fourier amplitudes on wave-vector by $\hat{u}_{\parallel} := (\hat{u},\hat{v})\cdot \frac{\bm k}{k} = \frac{\hat{u} k_x + \hat{v} k_y}{k} $ 

    Combining this to $\hat{\eta} = a^+ + a^-$ leads to 

    $$ a^+ = \frac{1}{2} (\hat{\eta} + \frac{\hat{u}_{\parallel}}{\sqrt{gk} e^{kz} }), \qquad a^- = \frac{1}{2} (\hat{\eta} - \frac{\hat{u}_{\parallel}}{\sqrt{gk} e^{kz} })$$


    Normalization is such that Parseval equality reads: 
    $$ \frac{1}{\cal A} \int_{\cal A} \eta^2(x,y) dx dy = \int | \hat{\eta}(k_x,k_y) |^2 d k_x d k_y$$

    Parameters
    ----------
    eta : (ny, nx) array
        Surface elevation [m].
    u, v : (ny, nx) arrays
        Horizontal velocities at z=0 [m/s].
    dx, dy : float
        Horizontal grid spacing [m].
    g : float
        Gravitational acceleration [m/s^2].
    z : (nz,) array, optional
        Vertical coordinates for a subsurface velocity field,
        with z=0 at the mean surface and z<0 underwater.
        Not used by this function; surface velocities are expected.
    w : (ny, nx) array, optional
        Surface vertical velocity [m/s], for diagnostics.

    Returns
    -------
    dict
        kx, ky, k, theta, omega:
            Wavenumber, direction and frequency grids.
        eta_hat:
            Fourier coefficients of elevation.
        a_plus, a_minus:
            Complex amplitudes of the two propagation branches.
        E_plus, E_minus:
            Modal elevation powers [m^2].
        diagnostics:
            Useful linear-wave consistency metrics.
    """
    eta = np.asarray(eta, dtype=float)
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)

    if eta.shape != u.shape or eta.shape != v.shape:
        raise ValueError("eta, u, and v must have identical shapes.")

    ny, nx = eta.shape

    # Remove means to eliminate the zero-wavenumber component.
    eta = eta - np.mean(eta)
    u = u - np.mean(u)
    v = v - np.mean(v)

    # Physical angular wavenumbers [rad/m].
    kx_1d = 2 * np.pi * np.fft.fftfreq(nx, d=dx)
    ky_1d = 2 * np.pi * np.fft.fftfreq(ny, d=dy)
    kx, ky = np.meshgrid(kx_1d, ky_1d)

    k = np.sqrt(kx**2 + ky**2)
    omega = np.sqrt(g * k)

    # Fourier-series coefficients: 
    # normalization such that 1/A int_A eta^2 dxdy = int abs(eta_hat)^2 dkx dky 
    
    norm = 1/(2*np.pi)*np.sqrt(dx*dy/(nx*ny))
    eta_hat = np.fft.fft2(eta) * norm
    u_hat = np.fft.fft2(u) *norm
    v_hat = np.fft.fft2(v) *norm

    # Longitudinal velocity in Fourier space.
    u_long_hat = np.zeros_like(eta_hat)
    nonzero = k > 0

    u_long_hat[nonzero] = (
        kx[nonzero] * u_hat[nonzero]
        + ky[nonzero] * v_hat[nonzero]
    ) / k[nonzero]

    # Branch decomposition.
    factor = np.zeros_like(k)
    factor[nonzero] = omega[nonzero] / (
        g * k[nonzero]
    )

    a_plus = 0.5 * (
        eta_hat + factor * u_long_hat
    )
    a_minus = 0.5 * (
        eta_hat - factor * u_long_hat
    )

    # Elevation variance contributions per Fourier mode.
    E_plus = np.abs(a_plus)**2
    E_minus = np.abs(a_minus)**2

    # Direction of the wavevector, measured from +x.
    theta = np.arctan2(ky, kx)

    # Optional vertical-velocity consistency diagnostic.
    diagnostics = {}
    # if w is not None:
    #     w = np.asarray(w, dtype=float)
    #     if w.shape != eta.shape:
    #         raise ValueError("w must have the same shape as eta.")

    #     w_hat = np.fft.fft2(w - np.mean(w)) / (nx * ny)

    #     mask = nonzero
    #     denominator = np.sqrt(np.sum(np.abs(u_long_hat[mask])**2))
    #     numerator = np.sqrt(
    #         np.sum(np.abs(w_hat[mask] + 1j*u_long_hat[mask])**2)
    #     )
    #     diagnostics["relative_w_residual"] = (
    #         numerator / denominator if denominator > 0 else np.nan
    #     )

    return {
        "kx": kx,
        "ky": ky,
        "k": k,
        "theta": theta,
        "omega": omega,
        "eta_hat": eta_hat,
        "u_long_hat": u_long_hat,
        "a_plus": a_plus,
        "a_minus": a_minus,
        "E_plus": E_plus,
        "E_minus": E_minus,
        "diagnostics": diagnostics,
    }



def bin_directional_spectrum(
    result,
    nk=40,
    ntheta=72,
    kmin=None,
    kmax=None,
):
    """
    Bin the branch amplitudes into a directional wavenumber spectrum.

    Returns
    -------
    k_edges : array
        Radial wavenumber bin edges [rad/m].
    theta_edges : array
        Direction bin edges [rad].
    spectrum : (nk, ntheta) array
        Modal elevation variance per (k, theta) bin [m^2].
        Divide by bin area dk*dtheta to obtain a density.
    """
    k = result["k"]
    kx = result["kx"]
    ky = result["ky"]
    Eplus = result["E_plus"]
    Eminus = result["E_minus"]
    theta = result["theta"]

    # One representative of each conjugate pair.
    half_plane = (kx > 0) | ((kx == 0) & (ky > 0))
    valid = half_plane & (k > 0)

    if kmin is None:
        kmin = np.min(k[valid])
    if kmax is None:
        kmax = np.max(k[valid])

    k_edges = np.linspace(kmin, kmax, nk + 1)
    theta_edges = np.linspace(-np.pi, np.pi, ntheta + 1)

    kp = k[valid]
    thp = theta[valid]

    # Forward and reverse branches, each counted once.
    k_all = np.concatenate([kp, kp])
    theta_all = np.concatenate([
        thp,
        (thp + np.pi + np.pi) % (2*np.pi) - np.pi,
    ])
    energy_all = np.concatenate([
        Eplus[valid],
        Eminus[valid],
    ])

    spectrum, _, _ = np.histogram2d(
        k_all,
        theta_all,
        bins=[k_edges, theta_edges],
        weights=energy_all,
    )

    dk = np.diff(k_edges)[:, None]
    dtheta = np.diff(theta_edges)[None, :]
    spectrum_density = spectrum / (dk * dtheta)

    return k_edges, theta_edges, spectrum_density


def plot_directional_spectrum(
    k_edges,
    theta_edges,
    spectrum_density,
    cmap="viridis",
):
    """
    Plot directional wavenumber spectrum in polar coordinates.

    Parameters
    ----------
    k_edges : array (nk + 1,)
        Wavenumber bin edges [rad/m].
    theta_edges : array (ntheta + 1,)
        Direction bin edges [rad].
    spectrum_density : array (nk, ntheta)
        Directional spectrum density.
    """
    fig, ax = plt.subplots(
        subplot_kw={"projection": "polar"},
        constrained_layout=True,
    )

    mesh = ax.pcolormesh(
        theta_edges,
        k_edges,
        spectrum_density,
        shading="auto",
        cmap=cmap,
    )

    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)
    ax.set_title("Directional wavenumber spectrum")

    cbar = fig.colorbar(mesh, ax=ax, pad=0.1)
    cbar.set_label(r"$E(k,\theta)$")

    return fig, ax