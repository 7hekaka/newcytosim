# analysis/plotting.py
import os
import matplotlib.pyplot as plt

def _maybe_minutes(t, minutes: bool):
    return (t / 60.0) if minutes else t

def make_plots(res, minutes: bool = False, save_prefix: str | None = None):
    # Top bias (per-frame)
    plt.figure()
    t = _maybe_minutes(res["t_perframe"], minutes)
    plt.plot(t, res["top_bias"])
    plt.xlabel("time (min)" if minutes else "time (s)")
    plt.ylabel("Top-bias index")
    plt.title("Top-bias vs time")
    if save_prefix:
        plt.savefig(f"{save_prefix}_top_bias.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Nematic order (per-frame)
    plt.figure()
    t = _maybe_minutes(res["t_perframe"], minutes)
    plt.plot(t, res["nematic_xy"])
    plt.xlabel("time (min)" if minutes else "time (s)")
    plt.ylabel("Nematic order (XY)")
    plt.title("Nematic order vs time")
    if save_prefix:
        plt.savefig(f"{save_prefix}_nematic_xy.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Radial / azimuthal flux (frame-to-frame)
    for key, label in [
        ("azimuthal_net", "Azimuthal throughput (net)"),
        ("azimuthal_abs", "Azimuthal throughput (abs)"),
        ("radial_net", "Radial flux (net)"),
        ("radial_abs", "Radial flux (abs)"),
    ]:
        plt.figure()
        t = _maybe_minutes(res["t_flux"], minutes)
        plt.plot(t, res[key])
        plt.xlabel("time (min)" if minutes else "time (s)")
        plt.ylabel(label)
        plt.title(f"{label} vs time")
        if save_prefix:
            safe_key = key.replace(" ", "_")
            plt.savefig(f"{save_prefix}_{safe_key}.png", dpi=200, bbox_inches="tight")
        plt.close()

    # Envelope / voxel union (per-frame)
    plt.figure()
    t = _maybe_minutes(res["t_perframe"], minutes)
    plt.plot(t, res["annulus_zspan"])
    plt.xlabel("time (min)" if minutes else "time (s)")
    plt.ylabel("Annulus envelope Z-span")
    plt.title("Annulus envelope Z-span vs time")
    if save_prefix:
        plt.savefig(f"{save_prefix}_annulus_zspan.png", dpi=200, bbox_inches="tight")
    plt.close()

    plt.figure()
    t = _maybe_minutes(res["t_perframe"], minutes)
    plt.plot(t, res["voxel_union"])
    plt.xlabel("time (min)" if minutes else "time (s)")
    plt.ylabel("Voxelized union fraction")
    plt.title("Voxelized union vs time")
    if save_prefix:
        plt.savefig(f"{save_prefix}_voxel_union.png", dpi=200, bbox_inches="tight")
    plt.close()

        # Top fraction (banded, 0–1 scale)
    if "top_fraction_band" in res:
        plt.figure()
        t = _maybe_minutes(res["t_perframe"], minutes)
        plt.plot(t, res["top_fraction_band"])
        plt.xlabel("time (min)" if minutes else "time (s)")
        plt.ylabel("Top fraction (banded)")
        plt.title("Top fraction (top vs bottom 1/3 bands)")
        if save_prefix:
            plt.savefig(f"{save_prefix}_top_fraction_band.png", dpi=200, bbox_inches="tight")
        plt.close()

    # Voxelized union per filament (size‐normalized occupancy)
    if "voxel_union_per_fil" in res:
        plt.figure()
        t = _maybe_minutes(res["t_perframe"], minutes)
        plt.plot(t, res["voxel_union_per_fil"])
        plt.xlabel("time (min)" if minutes else "time (s)")
        plt.ylabel("Voxelized occupancy / filament")
        plt.title("Voxelized union per filament vs time")
        if save_prefix:
            plt.savefig(f"{save_prefix}_voxel_union_per_fil.png", dpi=200, bbox_inches="tight")
        plt.close()

