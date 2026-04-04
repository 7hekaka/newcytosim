from __future__ import annotations
from typing import Dict, List, Tuple
import numpy as np
from .metrics import (
    top_bias_index, azimuthal_throughput, radial_flux,
    nematic_order_xy, annulus_envelope_zspan,
    voxelized_union_fraction, voxelized_union_per_filament,
    compute_top_bottom_bands,
)


def _unpack(item):
    """Normalize different frame payloads into vectors."""
    if isinstance(item, (list, tuple)) and item and isinstance(item[0], (list, tuple, np.ndarray)):
        item = np.asarray(item)

    if isinstance(item, np.ndarray) and item.ndim == 2:
        cols = item.shape[1]
        if cols == 6:  # frame, fid, x, y, z, w
            fid = item[:, 1].astype(np.int64, copy=False)
            x, y, z, w = item[:, 2], item[:, 3], item[:, 4], item[:, 5]
            return fid, x, y, z, w
        elif cols == 5:
            # fid, x, y, z, w  or frame, fid, x, y, z
            if np.issubdtype(item[:, 0].dtype, np.number) and np.all(item[:, 0] > 1e3):
                # probably frame numbers
                fid = item[:, 1].astype(np.int64, copy=False)
                x, y, z = item[:, 2], item[:, 3], item[:, 4]
                return fid, x, y, z, None
            fid = item[:, 0].astype(np.int64, copy=False)
            x, y, z, w = item[:, 1], item[:, 2], item[:, 3], item[:, 4]
            return fid, x, y, z, w
        elif cols == 4:  # fid, x, y, z
            fid = item[:, 0].astype(np.int64, copy=False)
            x, y, z = item[:, 1], item[:, 2], item[:, 3]
            return fid, x, y, z, None
        else:
            raise ValueError(f"Unsupported 2D frame with {cols} columns")

    if isinstance(item, (list, tuple)):
        if len(item) == 6:
            _frame, fid, x, y, z, w = item
            return np.array([fid]), np.array([x]), np.array([y]), np.array([z]), np.array([w])
        if len(item) == 5:
            fid, x, y, z, w = item
            return np.array([fid]), np.array([x]), np.array([y]), np.array([z]), np.array([w])
        if len(item) == 4:
            fid, x, y, z = item
            return np.array([fid]), np.array([x]), np.array([y]), np.array([z]), None

    raise ValueError("Unsupported frame payload.")


def analyze_from_frames(
    frames: List[tuple],
    dt: float,
    nfil: int = 0,
    drop_first: bool = False,
    rbins: np.ndarray | None = None,
) -> Dict[str, np.ndarray]:
    """Compute all frame-based and frame-to-frame metrics."""
    if not frames:
        return {k: np.array([]) for k in [
            "top_bias","top_fraction_band","nematic_xy","annulus_zspan",
            "voxel_union","voxel_union_per_fil","azimuthal_net","azimuthal_abs",
            "radial_net","radial_abs","t_perframe","t_flux",
        ]}

    # ---- per-frame metrics
    top_bias, top_fraction_band = [], []
    nematic_xy, annulus_zspan_vals = [], []
    voxel_union, voxel_union_per_fil = [], []

    for item in frames:
        fid, x, y, z, w = _unpack(item)
        r = np.hypot(x, y)

        # ---- top/bottom 1/3 bands
        bands = compute_top_bottom_bands(z, top_band_frac=1/3)
        top_bias.append(bands.top_bias)
        top_fraction_band.append(bands.top_fraction_band)

        # ---- nematic order in xy plane
        nematic_xy.append(nematic_order_xy(np.stack([fid, x, y], axis=1)))

        # ---- annulus envelope
        if rbins is None:
            rmax = float(np.max(r)) if r.size else 1.0
            rb = np.linspace(0.0, rmax + 1e-9, 16)
        annulus_zspan_vals.append(annulus_envelope_zspan(r, z, rb))

        # ---- voxelized occupancy
        vu = voxelized_union_fraction(x, y, z)
        voxel_union.append(vu)
        voxel_union_per_fil.append(voxelized_union_per_filament(x, y, z, nfil=nfil))

    # convert lists → arrays
    top_bias = np.asarray(top_bias)
    top_fraction_band = np.asarray(top_fraction_band)
    nematic_xy = np.asarray(nematic_xy)
    annulus_zspan_vals = np.asarray(annulus_zspan_vals)
    voxel_union = np.asarray(voxel_union)
    voxel_union_per_fil = np.asarray(voxel_union_per_fil)

    # optionally drop frame 0
    if drop_first and top_bias.size:
        top_bias = top_bias[1:]
        top_fraction_band = top_fraction_band[1:]
        nematic_xy = nematic_xy[1:]
        annulus_zspan_vals = annulus_zspan_vals[1:]
        voxel_union = voxel_union[1:]
        voxel_union_per_fil = voxel_union_per_fil[1:]
        perframe_offset = 1
    else:
        perframe_offset = 0

    # ---- frame-to-frame metrics
    azi_net, azi_abs, r_net, r_abs = [], [], [], []
    for k in range(1, len(frames)):
        _fid0, x0, y0, z0, _w0 = _unpack(frames[k-1])
        _fid1, x1, y1, z1, _w1 = _unpack(frames[k])
        net, ab = azimuthal_throughput((x0, y0), (x1, y1), dt)
        azi_net.append(net)
        azi_abs.append(ab)
        r0 = np.hypot(x0, y0)
        r1 = np.hypot(x1, y1)
        rn, ra = radial_flux(r0, r1, dt)
        r_net.append(rn)
        r_abs.append(ra)

    azi_net = np.asarray(azi_net)
    azi_abs = np.asarray(azi_abs)
    r_net = np.asarray(r_net)
    r_abs = np.asarray(r_abs)

    # ---- time vectors
    t_perframe = np.arange(perframe_offset, perframe_offset + top_bias.size, dtype=float) * dt
    t_flux = np.arange(1, 1 + azi_net.size, dtype=float) * dt

    return {
        "top_bias": top_bias,
        "top_fraction_band": top_fraction_band,
        "nematic_xy": nematic_xy,
        "annulus_zspan": annulus_zspan_vals,
        "voxel_union": voxel_union,
        "voxel_union_per_fil": voxel_union_per_fil,
        "azimuthal_net": azi_net,
        "azimuthal_abs": azi_abs,
        "radial_net": r_net,
        "radial_abs": r_abs,
        "t_perframe": t_perframe,
        "t_flux": t_flux,
    }
