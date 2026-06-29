from __future__ import annotations

import math
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "turnover_continuous_supply_smoke" / "top_pulsed_rotatable_fast"

SEED = 20260616

INNER = 2.625
OUTER = 2.875
TOP = 7.5
BOTTOM = -7.5
CHEWER_TOP = -3.0
FILAMENT_LEN = 2.5

INITIAL_FILAMENTS = 12
PULSE_FILAMENTS = 3
N_PULSES = 7
INITIAL_CROSSLINKERS = 48
CROSSLINKERS_PER_PULSE = 12
N_CHEWERS = 48
N_CLUSTERS = 4
MOTORS_PER_CLUSTER = 12

RUN_SECONDS = 60.0
TIME_STEP = 0.004
RUN_STEPS = int(RUN_SECONDS / TIME_STEP)
FRAMES_PER_BLOCK = 20


def annulus_point(rng: random.Random, radius: float | None = None) -> tuple[float, float]:
    theta = rng.uniform(0.0, 2.0 * math.pi)
    if radius is None:
        radius = rng.uniform(INNER + 0.05, OUTER - 0.05)
    return radius * math.cos(theta), radius * math.sin(theta)


def actin_line(rng: random.Random, z_minus_low: float, z_plus_high: float) -> str:
    z_plus = rng.uniform(z_minus_low + FILAMENT_LEN, z_plus_high)
    z_minus = z_plus - FILAMENT_LEN
    x, y = annulus_point(rng)
    return (
        f"new 1 actin {{ position_ends = {x:.3f} {y:.3f} {z_minus:.3f}, "
        f"{x:.3f} {y:.3f} {z_plus:.3f}; end_state = green, white; }}"
    )


def motor_lines(rng: random.Random) -> list[str]:
    lines = ["% --- BEGIN rotatable membrane-bound motor clusters ---"]
    cluster_z = [-5.6, -1.9, 1.9, 5.6]
    for cluster, z0 in enumerate(cluster_z[:N_CLUSTERS]):
        theta = rng.uniform(0.0, 2.0 * math.pi)
        lines.append(f"% cluster {cluster:02d}: theta={theta:.4f} z={z0:.3f}")
        for _ in range(MOTORS_PER_CLUSTER):
            dtheta = rng.gauss(0.0, 0.045)
            z = max(BOTTOM + 0.2, min(TOP - 0.2, rng.gauss(z0, 0.18)))
            radius = INNER + 0.01
            x = radius * math.cos(theta + dtheta)
            y = radius * math.sin(theta + dtheta)
            lines.append(f"new 1 myosin1 {{ position = {x:.3f} {y:.3f} {z:.3f} }}")
    lines.append("% --- END rotatable membrane-bound motor clusters ---")
    return lines


def run_block(label: str) -> list[str]:
    return [
        f"% {label}",
        f"run {RUN_STEPS} system",
        "{",
        f"    nb_frames = {FRAMES_PER_BLOCK}",
        "}",
        "",
    ]


def build_config() -> str:
    rng = random.Random(SEED)
    lines: list[str] = [
        "% Continuous top-supply smoke test: pulsed aligned actin birth + bottom minus-end chewing.",
        "% Reduced annulus: half linear dimensions of the 20 um mini system.",
        "% New top filaments are inserted every 60 s; crosslinkers are replenished with each pulse.",
        "% Fast-local variant: steric interactions are disabled for quick mechanism testing.",
        "set simul system",
        "{",
        "    steric = 0",
        f"    time_step = {TIME_STEP}",
        "    kT = 0.0042",
        "    viscosity = 0.1",
        "}",
        "",
        "set system display",
        "{",
        "    back_color = white",
        "}",
        "",
        "set space stripbox",
        "{",
        "    shape = annulus",
        "    display = ( color = blue; )",
        "}",
        "",
        "new stripbox",
        "{",
        f"    outer = {OUTER}",
        f"    inner = {INNER}",
        f"    top = {TOP}",
        f"    bottom = {BOTTOM}",
        "}",
        "",
        "set space chewer_zone",
        "{",
        "    shape = annulus",
        "    display = ( visible = 0; color = red; )",
        "}",
        "",
        "new chewer_zone",
        "{",
        f"    outer = {OUTER}",
        f"    inner = {INNER}",
        f"    top = {CHEWER_TOP}",
        f"    bottom = {BOTTOM}",
        "}",
        "",
        "set hand actin_binder",
        "{",
        "    binding_rate = 10",
        "    binding_range = 0.150",
        "    unbinding_rate = 0.08",
        "    unbinding_force = 5",
        "    display = ( color = cyan; )",
        "}",
        "",
        "set couple crosslinker",
        "{",
        "    hand1 = actin_binder",
        "    hand2 = actin_binder",
        "    stiffness = 10",
        "    diffusion = 10",
        "    specificity = parallel",
        "}",
        "",
        "set fiber actin",
        "{",
        "    rigidity = 0.075",
        "    segmentation = 0.18",
        "    steric = 0",
        "    confine = inside, 200",
        "",
        "    activity = grow",
        "    growing_speed = 0.15, 0",
        "    growing_force = inf, inf",
        "    total_polymer = 192.0",
        "    min_length = 0.18",
        "    max_chewing_speed = 1.5",
        "",
        "    display = ( color = black; )",
        "}",
        "",
        "set hand chewer1",
        "{",
        "    binding_rate = 50.0",
        "    binding_range = 0.15",
        "    unbinding_rate = 0.2",
        "    activity = chew",
        "    bind_also_end = minus_end",
        "    bind_only_end = minus_end, 0.5",
        "    hold_growing_end = 0, 1",
        "    chewing_speed = 0.8",
        "    diffusion = 0.2",
        "    display = ( color = red; size = 4; )",
        "}",
        "",
        "set single cutter1",
        "{",
        "    hand = chewer1",
        "    stiffness = 1",
        "    activity = diffuse",
        "    diffusion = 10.0",
        "    confine = inside, 0, chewer_zone",
        "}",
        "",
        "set hand motor1",
        "{",
        "    binding_rate = 10",
        "    binding_range = 0.02",
        "    unbinding_rate = 0.3",
        "    activity = move",
        "    unloaded_speed = 2.0",
        "    stall_force = 4",
        "    unbinding_force = 2.5",
        "",
        "    display = ( color = blue; size = 4; )",
        "}",
        "",
        "set single myosin1",
        "{",
        "    hand = motor1",
        "    stiffness = 100",
        "    diffusion = 10",
        "    activity = fixed",
        "}",
        "",
        "% Initial top-biased aligned filaments: minus end lower, plus/growing end higher.",
    ]

    for _ in range(INITIAL_FILAMENTS):
        lines.append(actin_line(rng, BOTTOM + 0.4, TOP - 0.2))

    lines.extend(
        [
            "",
            f"new {N_CHEWERS} cutter1 {{ position = inside, chewer_zone; }}",
            f"new {INITIAL_CROSSLINKERS} crosslinker {{ position = inside; }}",
            "",
        ]
    )
    lines.extend(motor_lines(rng))
    lines.append("")
    lines.extend(run_block("initial 60 s with motors, crosslinkers, and bottom chewers"))

    for pulse in range(1, N_PULSES + 1):
        lines.append(f"% --- top actin supply pulse {pulse:02d}: +{PULSE_FILAMENTS} filaments ---")
        for _ in range(PULSE_FILAMENTS):
            lines.append(actin_line(rng, TOP - FILAMENT_LEN - 0.25, TOP - 0.05))
        lines.append(f"new {CROSSLINKERS_PER_PULSE} crosslinker {{ position = inside; }}")
        lines.append("")
        lines.extend(run_block(f"post-pulse {pulse:02d}: 60 s"))

    return "\n".join(lines)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    config = build_config()
    (OUTDIR / "config.cym").write_text(config + "\n", encoding="utf-8")
    (OUTDIR / "README.md").write_text(
        "\n".join(
            [
                "# Continuous top-supply smoke test",
                "",
                "Fast-local reduced annulus with rotatable membrane-bound motors, crosslinkers,",
                "bottom minus-end chewers, and pulsed top-biased actin insertion every 60 seconds.",
                "Steric interactions are disabled in this variant so the mechanism can be tested quickly.",
                "",
                f"- annulus: inner `{INNER}` um, outer `{OUTER}` um, z `{BOTTOM}` to `{TOP}` um",
                f"- chewer zone: z `{BOTTOM}` to `{CHEWER_TOP}` um",
                f"- initial filaments: `{INITIAL_FILAMENTS}`",
                f"- pulse supply: `{PULSE_FILAMENTS}` filaments/min for `{N_PULSES}` pulses",
                f"- crosslinkers: `{INITIAL_CROSSLINKERS}` initial, `{CROSSLINKERS_PER_PULSE}` per pulse",
                f"- chewers: `{N_CHEWERS}`",
                f"- motors: `{N_CLUSTERS}` clusters x `{MOTORS_PER_CLUSTER}` motors",
                "",
                "Run locally with:",
                "",
                "```bash",
                "../../build_mini_turnover/bin/sim config.cym",
                "```",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(OUTDIR / "config.cym")


if __name__ == "__main__":
    main()
