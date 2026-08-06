"""
# Breaking front detection for a multilayer simulation

This script tests the breaking detection method from Wu et al. 2023

"""

import xarray as xr
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# add libpy
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../libpy/")
sys.path.append(filename)
from breaking import detect_ridges, detect_slope, detect_slope, simple_mapping

# ==================================
# Defining some know surfaces
# ==================================


ny, nx = 64, 64
A = 5.0
w = 4.0


def nan_eta():
    arr = np.zeros((32, 32), dtype=float)
    arr[10, 10] = np.nan
    return arr


def flat_eta():
    """No ridge anywhere -> curvature should be (numerically) zero everywhere,
    and simple_mapping should return an all-zero histogram."""
    return np.zeros((ny, nx), dtype=float)


def gaussian_ridge_eta():
    """
    A ridge elongated along the row (y) axis, centred at column x0=32,
    with amplitude A and width w. eta[i, j] = A * exp(-(j - x0)^2 / (2 w^2)).

    Analytic expectation at the crest (j == x0), away from edges:
      - curvature along the ridge (row direction)  ~ 0
      - curvature across the ridge (column direction) = -A / w^2  (peak of a
        Gaussian: d2/dx2 [A*exp(-x^2/2w^2)] at x=0 is -A/w^2)
    Since hessian_matrix_eigvals returns eigenvalues sorted descending,
    we expect maxima_ridges[i, x0] ~ 0 and minima_ridges[i, x0] ~ -A/w^2
    (up to the smoothing introduced by the Gaussian derivative sigma).
    """
    x0 = nx // 2
    x = np.arange(nx)
    profile = A * np.exp(-((x - x0) ** 2) / (2 * w**2))
    eta = np.tile(profile, (ny, 1))
    return {"eta": eta, "A": A, "w": w, "x0": x0, "ny": ny, "nx": nx}


def radial_gaussian_bump_eta():
    """
    A radially symmetric bump: eta = A * exp(-r^2 / (2 w^2)).
    At the centre, both principal curvatures are equal and negative
    (concave down in every direction): -A / w^2. This is a stricter,
    axis-independent check than the ridge fixture above.
    """
    cy, cx = 32, 32
    y, x = np.mgrid[0:ny, 0:nx]
    r2 = (x - cx) ** 2 + (y - cy) ** 2
    eta = A * np.exp(-r2 / (2 * w**2))
    return {"eta": eta, "A": A, "w": w, "cy": cy, "cx": cx}


def uniform_velocity_fields():
    """Constant ux, uy fields with a known speed, for simple_mapping
    histogram checks (so ux_a**2 + uy_a**2 has one predictable value
    wherever a_ == 1)."""
    ux = np.full((ny, nx), 3.0)
    uy = np.full((ny, nx), 4.0)  # speed = 5.0 wherever a_ == 1 (3-4-5 triangle)
    return {"ux": ux, "uy": uy, "expected_speed": 5.0}


nan_eta = nan_eta()
flat_eta = flat_eta()
gaussian_ridge_eta = gaussian_ridge_eta()
radial_gaussian_bump_eta = radial_gaussian_bump_eta()
uniform_velocity_fields = uniform_velocity_fields()


x = np.arange(nx)
fig, axes = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True, dpi=100)
ax = axes.flatten()
ax[0].pcolormesh(x, x, gaussian_ridge_eta["eta"], vmin=0, vmax=A)
ax[0].set_title("gaussian eta (bug visuel contour)")
s = ax[1].pcolormesh(x, x, radial_gaussian_bump_eta["eta"], vmin=0, vmax=A)
ax[1].set_title("radial gaussian eta")
plt.colorbar(s, ax=ax)

jslice = 0  # ny // 2
fig2, axes2 = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True, dpi=100)
ax2 = axes2.flatten()
ax2[0].plot(x, gaussian_ridge_eta["eta"][jslice, :])
ax2[0].set_ylim([0, A + 1])
ax2[0].set_title(f"gaussian 1D (j={jslice})")
ax2[1].plot(x, radial_gaussian_bump_eta["eta"][ny // 2, :], c="k", label="mid")
ax2[1].set_ylim([0, A + 1])
ax2[1].plot(
    x, radial_gaussian_bump_eta["eta"][ny // 2 + 2, :], c="darkblue", label="ny//2 + 2"
)
ax2[1].set_ylim([0, A + 1])
ax2[1].set_title("gaussian 2D")
# =================================
# testing detect_ridges()
# =================================


# 1 ) Defining some test
def test_output_shapes_and_dtype(flat_eta):
    maxima, minima = detect_ridges(flat_eta, sigma=1.0)
    assert maxima.shape == flat_eta.shape
    assert minima.shape == flat_eta.shape
    assert np.issubdtype(maxima.dtype, np.floating)
    assert np.issubdtype(minima.dtype, np.floating)


def test_flat_field_has_zero_curvature_everywhere(flat_eta):
    maxima, minima = detect_ridges(flat_eta, sigma=1.0)
    np.testing.assert_allclose(maxima, 0.0, atol=1e-10)
    np.testing.assert_allclose(minima, 0.0, atol=1e-10)


def test_maxima_ridges_are_never_smaller_than_minima_ridges(gaussian_ridge_eta):
    """
    Pins down the contract that `hessian_matrix_eigvals` returns eigenvalues
    sorted descending, i.e. maxima_ridges[i, j] >= minima_ridges[i, j]
    everywhere. This is the assumption BUG #3 depends on -- if this ever
    stops holding (e.g. a skimage version change), simple_mapping's use of
    `b` as "the strongly negative curvature marking a ridge" breaks silently.
    """
    eta = gaussian_ridge_eta["eta"]
    maxima, minima = detect_ridges(eta, sigma=1.0)
    assert np.all(maxima >= minima - 1e-12)


def test_ridge_crest_has_near_zero_along_ridge_curvature(gaussian_ridge_eta):
    """At the crest of a ridge elongated along rows, the along-ridge
    (row-direction) curvature should be ~0, which -- since eigenvalues are
    sorted descending -- should show up as maxima_ridges ~ 0 at the crest."""
    eta = gaussian_ridge_eta["eta"]
    x0 = gaussian_ridge_eta["x0"]
    ny = gaussian_ridge_eta["ny"]
    maxima, minima = detect_ridges(eta, sigma=1.0)
    interior = slice(10, ny - 10)  # avoid edge artefacts from Gaussian smoothing
    np.testing.assert_allclose(maxima[interior, x0], 0.0, atol=0.05)


def test_ridge_crest_has_strongly_negative_cross_ridge_curvature(gaussian_ridge_eta):
    """Cross-ridge curvature at the crest should be close to the analytic
    second derivative of a Gaussian at its peak: -A / w^2."""
    eta = gaussian_ridge_eta["eta"]
    A, w, x0, ny = (
        gaussian_ridge_eta["A"],
        gaussian_ridge_eta["w"],
        gaussian_ridge_eta["x0"],
        gaussian_ridge_eta["ny"],
    )
    maxima, minima = detect_ridges(eta, sigma=1.0)
    expected_curvature = -A / w**2
    interior = slice(10, ny - 10)
    # Gaussian-derivative smoothing (sigma=1) broadens the effective width
    # a little, so we allow a generous relative tolerance here rather than
    # pin an exact value -- tighten this once the analytic relationship
    # between `sigma` and effective smoothing is characterised separately.
    np.testing.assert_allclose(minima[interior, x0], expected_curvature, rtol=0.3)


def test_far_from_ridge_curvature_is_near_zero(gaussian_ridge_eta):
    eta = gaussian_ridge_eta["eta"]
    ny = gaussian_ridge_eta["ny"]
    maxima, minima = detect_ridges(eta, sigma=1.0)
    interior = slice(10, ny - 10)
    np.testing.assert_allclose(maxima[interior, 0], 0.0, atol=0.05)
    np.testing.assert_allclose(minima[interior, 0], 0.0, atol=0.05)


def test_radial_bump_is_isotropic_at_centre(radial_gaussian_bump_eta: dict):
    """At the exact centre of a radially symmetric bump, both principal
    curvatures are equal (isotropic), so maxima_ridges ~ minima_ridges
    there -- a stricter, axis-independent sanity check than the ridge test."""
    eta = radial_gaussian_bump_eta["eta"]
    cy, cx = radial_gaussian_bump_eta["cy"], radial_gaussian_bump_eta["cx"]
    A, w = radial_gaussian_bump_eta["A"], radial_gaussian_bump_eta["w"]
    maxima, minima = detect_ridges(eta, sigma=1.0)
    expected_curvature = -A / w**2
    np.testing.assert_allclose(maxima[cy, cx], minima[cy, cx], atol=0.05)
    np.testing.assert_allclose(minima[cy, cx], expected_curvature, rtol=0.3)


def test_runs_for_various_sigma(flat_eta, sigma):
    # Smoke test across the sigma range we expect to use in production.
    maxima, minima = detect_ridges(flat_eta, sigma=sigma)
    assert maxima.shape == flat_eta.shape


def test_nan_input_propagates_rather_than_crashing_silently(nan_eta):
    """Document current NaN behaviour explicitly -- undecided whether this
    should be the long-term contract, but the pipeline should never silently
    treat a NaN as a valid ridge without us deciding that on purpose."""
    maxima, minima = detect_ridges(nan_eta, sigma=1.0)
    assert np.isnan(maxima).any() or np.isnan(minima).any()


# 2) Run the tests

test_output_shapes_and_dtype(flat_eta)
test_flat_field_has_zero_curvature_everywhere(flat_eta)
test_maxima_ridges_are_never_smaller_than_minima_ridges(gaussian_ridge_eta)
test_ridge_crest_has_near_zero_along_ridge_curvature(gaussian_ridge_eta)
test_ridge_crest_has_strongly_negative_cross_ridge_curvature(gaussian_ridge_eta)
test_far_from_ridge_curvature_is_near_zero(gaussian_ridge_eta)
test_radial_bump_is_isotropic_at_centre(radial_gaussian_bump_eta)
for sigma in [0.5, 1.0, 2.0, 3.0]:
    test_runs_for_various_sigma(flat_eta, sigma)
test_nan_input_propagates_rather_than_crashing_silently(nan_eta)


gmax, gmin = detect_ridges(gaussian_ridge_eta["eta"])

rgmax, rgmin = detect_ridges(radial_gaussian_bump_eta["eta"])


ax[0].contour(x, x, gmin, colors="b")
ax[0].contour(x, x, gmax, colors="r")
ax2[0].scatter(
    x,
    np.where(gmin[ny // 2, :] < 0, gaussian_ridge_eta["eta"][jslice, :], np.nan),
    marker="x",
    c="cyan",
    label="curv < 0",
)
ax2[0].scatter(
    x,
    np.where(gmax[ny // 2, :] > 0, gaussian_ridge_eta["eta"][jslice, :], np.nan),
    marker="x",
    c="r",
    label="curv > 0",
)
ax2[0].scatter(
    x,
    np.where(gmin[ny // 2, :] < -0.05, gaussian_ridge_eta["eta"][jslice, :], np.nan),
    marker="x",
    c="k",
    label="curv < thresh = -0.05 ",
)

"""
Contours are the result of detect_ridges()
"""
ax[1].contour(x, x, rgmin, colors="b")
ax[1].contour(x, x, rgmax, colors="r")
ax[1].set_aspect(1)
ax[0].set_aspect(1)
ax2[0].legend()
ax2[1].legend()


print("Tests for detect_ridges() passed !")

# =================================
# testing simple_mapping()
# =================================


class TestFlatField:
    def test_no_ridge_gives_all_zero_histogram(self, flat_eta, uniform_velocity_fields):
        """
        Flat field should compute a 0. curvature which is > to threshold
        so no ridge should be detected !
        """
        hist = simple_mapping(
            flat_eta,
            uniform_velocity_fields["ux"],
            uniform_velocity_fields["uy"],
            delta=1.0,
            bins=BINS,
            method=0,
            threshold=-0.01,
            EXTRA_FILTER=False,
        )
        target = np.zeros(len(BINS) - 1)
        np.testing.assert_array_equal(hist, target)

    def test_extra_filter_toggle_does_not_crash_on_flat_field(
        self, flat_eta, uniform_velocity_fields
    ):
        for extra in (True, False):
            hist = simple_mapping(
                flat_eta,
                uniform_velocity_fields["ux"],
                uniform_velocity_fields["uy"],
                delta=1.0,
                bins=BINS,
                method=0,
                threshold=-0.01,
                EXTRA_FILTER=extra,
            )
            assert hist.shape == (len(BINS) - 1,)


class TestRidgeCase:
    """
    Cross-ridge crest detected via a strongly negative curvature, with a
    known uniform velocity field (speed = 5.0 wherever a front is flagged).
    All detected mass must land in the bin containing 5.0, nothing else.
    """

    def test_detected_front_velocity_lands_in_expected_bin(
        self, gaussian_ridge_eta, uniform_velocity_fields
    ):
        eta = gaussian_ridge_eta["eta"]
        hist, mask = simple_mapping(
            eta,
            uniform_velocity_fields["ux"],
            uniform_velocity_fields["uy"],
            delta=1.0,
            bins=BINS,
            method=0,
            threshold=-0.05,
            EXTRA_FILTER=False,
            return_mask=True,
        )
        assert hist.sum() > 0, "expected at least one detected front pixel"
        bin_5 = np.digitize([5.0], BINS)[0] - 1
        other_bins = [i for i in range(len(hist)) if i != bin_5]
        assert hist[bin_5] == hist.sum(), (
            "all detected-front pixels have speed exactly 5.0 (uniform ux=3, "
            "uy=4), so every count should fall in the bin containing 5.0"
        )
        for i in other_bins:
            assert hist[i] == 0

        return mask

    def test_threshold_that_excludes_ridge_gives_zero_histogram(
        self, gaussian_ridge_eta, uniform_velocity_fields
    ):
        eta = gaussian_ridge_eta["eta"]
        hist = simple_mapping(
            eta,
            uniform_velocity_fields["ux"],
            uniform_velocity_fields["uy"],
            delta=1.0,
            bins=BINS,
            method=0,
            threshold=-1000.0,  # nothing has curvature this negative
            EXTRA_FILTER=False,
        )
        np.testing.assert_array_equal(hist, np.zeros(len(BINS) - 1))


# run the tests
BINS = np.array([0.0, 2.0, 4.0, 6.0, 8.0])  # speed 5.0 falls in bin [4, 6)

flatfield = TestFlatField()
flatfield.test_no_ridge_gives_all_zero_histogram(flat_eta, uniform_velocity_fields)
flatfield.test_extra_filter_toggle_does_not_crash_on_flat_field(
    flat_eta, uniform_velocity_fields
)

ridgefield = TestRidgeCase()
gmask = ridgefield.test_detected_front_velocity_lands_in_expected_bin(
    gaussian_ridge_eta, uniform_velocity_fields
)
"""
This test shows a bug when using a 1D gradient (along x) on a 2D field.
the detected edge should be a straite line, but on the side spurious cells are detected.
This behavior is not seen in the test below
"""
rgmask = ridgefield.test_detected_front_velocity_lands_in_expected_bin(
    radial_gaussian_bump_eta, uniform_velocity_fields
)


ridgefield.test_threshold_that_excludes_ridge_gives_zero_histogram(
    gaussian_ridge_eta, uniform_velocity_fields
)
print("Tests simple_mappting() passed!")


"""
Mask of detected ridge is plotted on the figure in pink
"""

ax[0].pcolormesh(
    x,
    x,
    np.ma.masked_where(gmask, np.ones(gmask.shape)),
    cmap="cool",
    vmax=1,
    vmin=0,
)
ax[1].pcolormesh(
    x,
    x,
    np.ma.masked_where(rgmask, np.ones(gmask.shape)),
    cmap="cool",
    vmax=1,
    vmin=0,
)


fig.savefig("ridge_gaussian.png")
fig2.savefig("ridge_gaussian_slice.png")

print("check produced figure that everything is ok")

plt.show()
