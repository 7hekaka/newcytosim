from __future__ import annotations
from typing import Tuple
import numpy as np


from dataclasses import dataclass

@dataclass
class TopBottomBands:
    top_bias: float
    top_fraction_band: float
    n_top: int
    n_bot: int

def compute_top_bottom_bands(z: np.ndarray, top_band_frac: float = 1/3) -> TopBottomBands:
    """Compute top/bottom-band bias using top & bottom bands of thickness f.
       Returns both [-1,1] bias and [0,1] top-fraction among the two bands."""
    if z.size == 0:
        return TopBottomBands(np.nan, np.nan, 0, 0)

    zmin = float(np.min(z))
    zmax = float(np.max(z))
    span = max(zmax - zmin, 1e-12)  # avoid div-by-zero

    # band thresholds
    z_bot_max = zmin + top_band_frac * span
    z_top_min = zmin + (1.0 - top_band_frac) * span

    n_top = int(np.sum(z >= z_top_min))
    n_bot = int(np.sum(z <= z_bot_max))

    denom = (n_top + n_bot)
    if denom == 0:
        return TopBottomBands(np.nan, np.nan, n_top, n_bot)

    bias = (n_top - n_bot) / denom
    frac = n_top / denom
    return TopBottomBands(bias, frac, n_top, n_bot)

# ---------- helpers ----------
def unwrap_angles(a: np.ndarray) -> np.ndarray:
    return np.unwrap(a)


def theta_from_xy(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Azimuth angle in (−π, π]."""
    return np.arctan2(y, x)


def circular_mean_theta(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """
    Mean azimuth (−π..π) via atan2(<sin>, <cos>) and resultant length R (0..1).
    Returns (mean_theta, R).
    """
    if x.size == 0:
        return 0.0, 0.0
    ang = np.arctan2(y, x)
    s = np.mean(np.sin(ang))
    c = np.mean(np.cos(ang))
    return float(np.arctan2(s, c)), float(np.hypot(s, c))


# ---------- single-frame metrics ----------
def top_bias_index(z: np.ndarray, z0: float = 0.0) -> float:
    """
    Fractional excess of 'top' over 'bottom': (N_top - N_bot)/N_tot, in [-1,1].
    z0 is the midplane (default 0).
    """
    if z.size == 0:
        return 0.0
    top = np.count_nonzero(z > z0)
    bot = np.count_nonzero(z < z0)
    tot = top + bot
    if tot == 0:
        return 0.0
    return float((top - bot) / tot)


def nematic_order_xy(fidxy: np.ndarray) -> float:
    """
    2D nematic order S2 from segment tangents in xy within each filament.
    fidxy columns: 0=fid, 1=x, 2=y
    """
    if fidxy.size == 0:
        return 0.0

    order = np.lexsort((np.arange(fidxy.shape[0]), fidxy[:, 0]))
    data = fidxy[order]

    c2_sum = 0.0
    s2_sum = 0.0
    n = 0
    for k in range(1, data.shape[0]):
        if data[k, 0] != data[k - 1, 0]:
            continue
        dx = data[k, 1] - data[k - 1, 1]
        dy = data[k, 2] - data[k - 1, 2]
        norm = np.hypot(dx, dy)
        if norm <= 0:
            continue
        th = np.arctan2(dy, dx)
        c2_sum += np.cos(2 * th)
        s2_sum += np.sin(2 * th)
        n += 1
    if n == 0:
        return 0.0
    return float(np.hypot(c2_sum, s2_sum) / n)


def annulus_envelope_zspan(r: np.ndarray, z: np.ndarray, rbins: np.ndarray) -> float:
    """
    Point-weighted mean z-span across radial bins. Returns 0 if no points.
    """
    if r.size == 0:
        return 0.0
    idx = np.digitize(r, rbins) - 1
    zspan, w = [], []
    for b in range(len(rbins) - 1):
        mask = (idx == b)
        if np.any(mask):
            zspan.append(np.ptp(z[mask]))
            w.append(mask.sum())
    if not zspan:
        return 0.0
    zspan = np.array(zspan, dtype=float)
    w = np.array(w, dtype=float)
    return float(np.average(zspan, weights=w))


def voxelized_union_fraction(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    nr: int = 24,
    nth: int = 36,
    nz: int = 24,
) -> float:
    """
    Cylindrical voxel occupancy fraction within the frame's bounding box.
    Returns the fraction of (r,theta,z) voxels occupied by at least one point.
    """
    if x.size == 0:
        return 0.0

    r = np.hypot(x, y)
    rmax = float(np.max(r)) if r.size else 1.0
    zmin, zmax = (float(np.min(z)), float(np.max(z))) if z.size else (0.0, 1.0)

    rbins = np.linspace(0.0, rmax + 1e-9, nr + 1)
    zbins = np.linspace(zmin, zmax + 1e-9, nz + 1)
    thetas = np.arctan2(y, x)
    tbins = np.linspace(-np.pi, np.pi + 1e-9, nth + 1)

    ir = np.clip(np.digitize(r, rbins) - 1, 0, nr - 1)
    iz = np.clip(np.digitize(z, zbins) - 1, 0, nz - 1)
    it = np.clip(np.digitize(thetas, tbins) - 1, 0, nth - 1)

    occ = np.zeros((nr, nth, nz), dtype=bool)
    occ[ir, it, iz] = True
    return float(occ.mean())


# ---------- frame-to-frame metrics (throughput/flux) ----------
def azimuthal_throughput(
    prev_xy: Tuple[np.ndarray, np.ndarray],
    curr_xy: Tuple[np.ndarray, np.ndarray],
    dt: float,
) -> Tuple[float, float]:
    """
    Azimuthal 'throughput' between two frames sharing point identity/order.
    Returns (net, abs_mean) angular velocity in rad/s:
      - net      = mean(dtheta_unwrapped)/dt  (signed)
      - abs_mean = mean(|dtheta_unwrapped|)/dt
    """
    (x0, y0), (x1, y1) = prev_xy, curr_xy
    if x0.size == 0 or x1.size == 0 or dt <= 0:
        return 0.0, 0.0

    th0 = theta_from_xy(x0, y0)
    th1 = theta_from_xy(x1, y1)
    dth = unwrap_angles(th1) - unwrap_angles(th0)

    net = float(np.mean(dth) / dt)
    abs_mean = float(np.mean(np.abs(dth)) / dt)
    return net, abs_mean


def radial_flux(
    prev_r: np.ndarray,
    curr_r: np.ndarray,
    dt: float,
) -> Tuple[float, float]:
    """
    Radial flux between two frames sharing point identity/order.
    Returns (net, abs_mean) radial speed in length/s:
      - net      = mean(dr)/dt  (outward positive)
      - abs_mean = mean(|dr|)/dt
    """
    if prev_r.size == 0 or curr_r.size == 0 or dt <= 0:
        return 0.0, 0.0

    dr = curr_r - prev_r
    net = float(np.mean(dr) / dt)
    abs_mean = float(np.mean(np.abs(dr)) / dt)
    return net, abs_mean

def voxelized_union_per_filament(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    nfil: int,
    nr: int = 24,
    nth: int = 36,
    nz: int = 24,
) -> float:
    """
    Same as voxelized_union_fraction but normalized by number of filaments.
    Avoids scale differences when comparing different filament densities.
    """
    if nfil <= 0:
        return np.nan
    base = voxelized_union_fraction(x, y, z, nr=nr, nth=nth, nz=nz)
    return base / float(nfil)
