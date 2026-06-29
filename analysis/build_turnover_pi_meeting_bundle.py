#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "analysis" / "results"
FAIR = RESULTS / "turnover_fair_comparisons"
OUT = RESULTS / "pi_meeting_turnover_2026-04-23"
PLAY = ROOT / "turnover_growth_pool_bias_bottom_chewers" / "play"

CAMPAIGNS = [
    "turnover_growth_sweep",
    "turnover_growth_pool_bias_sweep",
    "turnover_growth_pool_bias_bottom_chewers",
    "turnover_growth_pool_bias_bottom_chewers_diffuse",
    "turnover_growth_pool_bias_motor_escalation",
]

CAMPAIGN_LABELS = {
    "turnover_growth_sweep": "baseline growth sweep",
    "turnover_growth_pool_bias_sweep": "pool-biased growth sweep",
    "turnover_growth_pool_bias_bottom_chewers": "fixed bottom chewers",
    "turnover_growth_pool_bias_bottom_chewers_diffuse": "diffuse bottom chewers",
    "turnover_growth_pool_bias_motor_escalation": "motor escalation",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def f(row: dict[str, str], key: str, default: float = float("nan")) -> float:
    try:
        value = row.get(key, "")
        if value == "":
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def fmt(value: float, digits: int = 3) -> str:
    if not math.isfinite(value):
        return "n/a"
    return f"{value:.{digits}f}"


def load_run_rows() -> dict[str, list[dict[str, str]]]:
    out = {}
    for campaign in CAMPAIGNS:
        out[campaign] = read_csv(RESULTS / campaign / "run_metrics.csv")
    return out


def find_run(
    run_rows: dict[str, list[dict[str, str]]],
    campaign: str,
    scenario: str,
    condition: str,
    *,
    total_motors: int | None = None,
) -> dict[str, str]:
    matches = []
    for row in run_rows[campaign]:
        if row["scenario"] != scenario or row["condition"] != condition:
            continue
        if total_motors is not None and int(float(row.get("total_motors", "0") or 0)) != total_motors:
            continue
        matches.append(row)
    if not matches:
        raise KeyError((campaign, scenario, condition, total_motors))
    return matches[0]


def top_rows(rows: list[dict[str, str]], key: str, n: int = 8, **filters: str) -> list[dict[str, str]]:
    filtered = []
    for row in rows:
        if any(row.get(field) != value for field, value in filters.items()):
            continue
        value = f(row, key)
        if math.isfinite(value):
            filtered.append(row)
    return sorted(filtered, key=lambda row: f(row, key), reverse=True)[:n]


def copy_asset(src: Path, dest_rel: str, missing: list[str]) -> str:
    dest = OUT / dest_rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.copy2(src, dest)
    else:
        missing.append(str(src))
    return dest_rel


def make_movie_row(
    run_rows: dict[str, list[dict[str, str]]],
    priority: int,
    story: str,
    campaign: str,
    scenario: str,
    condition: str,
    why_show: str,
    *,
    total_motors: int | None = None,
) -> dict[str, object]:
    row = find_run(run_rows, campaign, scenario, condition, total_motors=total_motors)
    run_path = Path(row["run_path"])
    play_cmd = f"{PLAY} {run_path} objects.cmo"
    render_cmd = f"{PLAY} {run_path} objects.cmo movie period=5"
    return {
        "priority": priority,
        "story": story,
        "campaign": campaign,
        "campaign_label": CAMPAIGN_LABELS[campaign],
        "scenario": scenario,
        "scenario_label": row["scenario_label"],
        "condition": condition,
        "condition_label": row["condition_label"],
        "run_dir": row["run_dir"],
        "run_path": str(run_path),
        "auc_vz_abs": fmt(f(row, "auc_vz_abs")),
        "excess_auc_vz_abs": fmt(f(row, "excess_auc_vz_abs")),
        "matched_active_transport_fraction": fmt(f(row, "matched_active_transport_fraction")),
        "total_motors": row.get("total_motors", ""),
        "why_show": why_show,
        "play_command": play_cmd,
        "render_command": render_cmd,
    }


def make_selected_movies(run_rows: dict[str, list[dict[str, str]]]) -> list[dict[str, object]]:
    return [
        make_movie_row(
            run_rows,
            1,
            "Main positive result",
            "turnover_growth_pool_bias_bottom_chewers_diffuse",
            "sc02_top_two_thirds_aligned",
            "c03_rotatable_xlink",
            "Strongest fair rotatable + xlink motor gain; best single trajectory to open the result.",
        ),
        make_movie_row(
            run_rows,
            2,
            "Matched passive control",
            "turnover_growth_pool_bias_bottom_chewers_diffuse",
            "sc02_top_two_thirds_aligned",
            "c01_nomotor_xlink",
            "Use immediately after the main movie to show the matched no-motor baseline.",
        ),
        make_movie_row(
            run_rows,
            3,
            "Rotatable versus fixed-global control",
            "turnover_growth_pool_bias_bottom_chewers_diffuse",
            "sc02_top_two_thirds_aligned",
            "c05_fixed_global_xlink",
            "Shows that the high transport is not just crosslinking; motor orientation/adaptability matters.",
        ),
        make_movie_row(
            run_rows,
            4,
            "Diffuse chewers help localized growth",
            "turnover_growth_pool_bias_bottom_chewers_diffuse",
            "sc04_top_cap_aligned",
            "c03_rotatable_xlink",
            "Largest diffuse-minus-fixed bottom-chewer AUC gain among matched full motor conditions.",
        ),
        make_movie_row(
            run_rows,
            5,
            "Fixed chewer comparator",
            "turnover_growth_pool_bias_bottom_chewers",
            "sc04_top_cap_aligned",
            "c03_rotatable_xlink",
            "Compare against priority 4 to show what changes when bottom chewers diffuse.",
        ),
        make_movie_row(
            run_rows,
            6,
            "Motor escalation upper bound",
            "turnover_growth_pool_bias_motor_escalation",
            "sc02_top_two_thirds_aligned",
            "c03_rotatable_xlink",
            "Highest-AUC top-two-thirds motor escalation condition; useful as an upper-bound phenotype.",
            total_motors=1920,
        ),
        make_movie_row(
            run_rows,
            7,
            "Pool-bias baseline for escalation",
            "turnover_growth_pool_bias_sweep",
            "sc02_top_two_thirds_aligned",
            "c03_rotatable_xlink",
            "Compare against priority 6 to show the effect of adding more motor capacity.",
        ),
        make_movie_row(
            run_rows,
            8,
            "Why raw comparisons needed controls",
            "turnover_growth_sweep",
            "sc02_top_two_thirds_aligned",
            "c03_rotatable_xlink",
            "Baseline growth-only version of the same scenario; use to motivate fair matched comparisons.",
        ),
    ]


def make_headline_findings(effect_rows: list[dict[str, str]], reference_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    motor_auc = top_rows(effect_rows, "delta_auc_vz_abs", effect="motor_gain_rotatable_xlink")[0]
    motor_excess = top_rows(effect_rows, "delta_excess_auc_vz_abs", effect="motor_gain_rotatable_xlink")[0]
    orientation = top_rows(effect_rows, "delta_auc_vz_abs", effect="orientation_gain_xlink")[0]
    xlink = top_rows(effect_rows, "delta_auc_vz_abs", effect="xlink_gain_rotatable")[0]
    diffuse = top_rows(
        reference_rows,
        "delta_auc_vz_abs",
        target_campaign="turnover_growth_pool_bias_bottom_chewers_diffuse",
        reference_campaign="turnover_growth_pool_bias_bottom_chewers",
        condition="c03_rotatable_xlink",
    )[0]
    escalation = top_rows(
        reference_rows,
        "delta_auc_vz_abs",
        target_campaign="turnover_growth_pool_bias_motor_escalation",
        reference_campaign="turnover_growth_pool_bias_sweep",
        condition="c03_rotatable_xlink",
    )[0]
    pool_bias = top_rows(
        reference_rows,
        "delta_auc_vz_abs",
        target_campaign="turnover_growth_pool_bias_sweep",
        reference_campaign="turnover_growth_sweep",
        condition="c03_rotatable_xlink",
    )[0]

    return [
        {
            "finding": "Best fair motor effect",
            "claim": "Rotatable motors with crosslinkers are strongest in top-two-thirds growth with diffuse bottom chewers.",
            "support": f"Delta AUC = {fmt(f(motor_auc, 'delta_auc_vz_abs'))} um vs matched no-motor + xlink.",
            "campaign": motor_auc["campaign"],
            "scenario": motor_auc["scenario"],
            "condition": motor_auc["numerator_condition"],
        },
        {
            "finding": "Best excess transport",
            "claim": "The same condition remains strongest after subtracting matched passive AUC.",
            "support": f"Delta excess AUC = {fmt(f(motor_excess, 'delta_excess_auc_vz_abs'))} um.",
            "campaign": motor_excess["campaign"],
            "scenario": motor_excess["scenario"],
            "condition": motor_excess["numerator_condition"],
        },
        {
            "finding": "Rotatable beats fixed-global",
            "claim": "Motor orientation/adaptability is a major determinant; fixed-global motors underperform in matched xlinked cases.",
            "support": f"Largest rotatable-minus-fixed AUC delta = {fmt(f(orientation, 'delta_auc_vz_abs'))} um.",
            "campaign": orientation["campaign"],
            "scenario": orientation["scenario"],
            "condition": orientation["numerator_condition"],
        },
        {
            "finding": "Crosslinkers matter",
            "claim": "Crosslinkers amplify rotatable motor transport rather than acting as a passive add-on.",
            "support": f"Largest rotatable xlink gain = {fmt(f(xlink, 'delta_auc_vz_abs'))} um.",
            "campaign": xlink["campaign"],
            "scenario": xlink["scenario"],
            "condition": xlink["numerator_condition"],
        },
        {
            "finding": "Diffuse chewers improve localized phenotypes",
            "claim": "Diffusing bottom chewers improve the top-cap rotatable+xlink trajectory relative to fixed chewers.",
            "support": (
                f"Diffuse-minus-fixed AUC delta = {fmt(f(diffuse, 'delta_auc_vz_abs'))} um; "
                f"matched active-fraction delta = {fmt(f(diffuse, 'delta_matched_active_transport_fraction'))}."
            ),
            "campaign": diffuse["target_campaign"],
            "scenario": diffuse["scenario"],
            "condition": diffuse["condition"],
        },
        {
            "finding": "Motor escalation is a strong upper bound",
            "claim": "Increasing motor capacity boosts all rotatable+xlink scenarios, but this campaign lacks matched no-motor controls.",
            "support": f"Top escalation-minus-pool-bias AUC delta = {fmt(f(escalation, 'delta_auc_vz_abs'))} um.",
            "campaign": escalation["target_campaign"],
            "scenario": escalation["scenario"],
            "condition": escalation["condition"],
        },
        {
            "finding": "Pool-biased growth explains part of the raw movement",
            "claim": "Top-biased growth can inflate raw movement, so matched controls are required for fair interpretation.",
            "support": f"Top pool-bias-minus-growth AUC delta = {fmt(f(pool_bias, 'delta_auc_vz_abs'))} um in the full motor condition.",
            "campaign": pool_bias["target_campaign"],
            "scenario": pool_bias["scenario"],
            "condition": pool_bias["condition"],
        },
    ]


def copy_figures() -> list[str]:
    missing: list[str] = []
    assets = [
        (
            FAIR / "within_campaign" / "auc_vz_abs" / "motor_gain_rotatable_xlink.png",
            "figures/01_within_motor_gain_rotatable_xlink_auc.png",
        ),
        (
            FAIR / "within_campaign" / "excess_auc_vz_abs" / "motor_gain_rotatable_xlink.png",
            "figures/02_within_motor_gain_rotatable_xlink_excess_auc.png",
        ),
        (
            FAIR / "within_campaign" / "matched_active_transport_fraction" / "motor_gain_rotatable_xlink.png",
            "figures/03_within_motor_gain_rotatable_xlink_matched_active_fraction.png",
        ),
        (
            FAIR / "within_campaign" / "auc_vz_abs" / "orientation_gain_xlink.png",
            "figures/04_within_orientation_gain_xlink_auc.png",
        ),
        (
            FAIR / "within_campaign" / "auc_vz_abs" / "xlink_gain_rotatable.png",
            "figures/05_within_xlink_gain_rotatable_auc.png",
        ),
        (
            FAIR / "cross_campaign" / "auc_vz_abs" / "turnover_growth_pool_bias_sweep__minus__turnover_growth_sweep.png",
            "figures/06_cross_pool_bias_minus_growth_auc.png",
        ),
        (
            FAIR / "cross_campaign" / "auc_vz_abs" / "turnover_growth_pool_bias_bottom_chewers_diffuse__minus__turnover_growth_pool_bias_bottom_chewers.png",
            "figures/07_cross_diffuse_minus_fixed_chewers_auc.png",
        ),
        (
            FAIR / "cross_campaign" / "matched_active_transport_fraction" / "turnover_growth_pool_bias_bottom_chewers_diffuse__minus__turnover_growth_pool_bias_bottom_chewers.png",
            "figures/08_cross_diffuse_minus_fixed_chewers_matched_active_fraction.png",
        ),
        (
            FAIR / "cross_campaign" / "bottom_half_mass_fraction_final" / "turnover_growth_pool_bias_bottom_chewers_diffuse__minus__turnover_growth_pool_bias_bottom_chewers.png",
            "figures/09_cross_diffuse_minus_fixed_chewers_bottom_half_mass.png",
        ),
        (
            FAIR / "cross_campaign" / "auc_vz_abs" / "turnover_growth_pool_bias_motor_escalation__minus__turnover_growth_pool_bias_sweep.png",
            "figures/10_cross_motor_escalation_minus_pool_bias_auc.png",
        ),
        (
            RESULTS / "turnover_growth_pool_bias_bottom_chewers_diffuse" / "kymographs" / "full_timeline_kymograph_grid.png",
            "figures/11_diffuse_bottom_chewers_full_timeline_kymograph_grid.png",
        ),
        (
            RESULTS / "turnover_growth_pool_bias_bottom_chewers_diffuse" / "timecourses" / "by_scenario" / "sc02_top_two_thirds_aligned_timecourses.png",
            "figures/12_diffuse_sc02_top_two_thirds_timecourses.png",
        ),
        (
            RESULTS / "turnover_growth_pool_bias_bottom_chewers_diffuse" / "timecourses" / "by_scenario" / "sc04_top_cap_aligned_timecourses.png",
            "figures/13_diffuse_sc04_top_cap_timecourses.png",
        ),
        (
            RESULTS / "turnover_growth_pool_bias_motor_escalation" / "timecourses" / "by_scenario" / "sc02_top_two_thirds_aligned_timecourses.png",
            "figures/14_motor_escalation_sc02_top_two_thirds_timecourses.png",
        ),
        (
            RESULTS / "turnover_growth_pool_bias_bottom_chewers_diffuse" / "unwrapped_annulus" / "heatmaps" / "t120_heatmaps.png",
            "figures/15_diffuse_unwrapped_heatmaps_t120.png",
        ),
        (
            RESULTS / "turnover_growth_pool_bias_bottom_chewers_diffuse" / "unwrapped_annulus" / "heatmaps" / "t720_heatmaps.png",
            "figures/16_diffuse_unwrapped_heatmaps_t720.png",
        ),
    ]
    copied = [copy_asset(src, dest, missing) for src, dest in assets]
    if missing:
        (OUT / "missing_assets.txt").write_text("\n".join(missing) + "\n", encoding="utf-8")
    return copied


def write_markdown(
    findings: list[dict[str, object]],
    movie_rows: list[dict[str, object]],
    copied_assets: list[str],
) -> None:
    main = movie_rows[0]
    control = movie_rows[1]
    fixed = movie_rows[2]
    topcap = movie_rows[3]
    fixed_topcap = movie_rows[4]
    escalation = movie_rows[5]
    pool_bias = movie_rows[6]

    summary = [
        "# Turnover One-Trajectory Sweep: PI Meeting Summary",
        "",
        "## Short version",
        (
            "The cleanest story is that top-biased actin renewal plus rotatable motors and crosslinkers can produce strong axial transport, "
            "but raw movement metrics are not fair by themselves because growth, turnover, and polymer diffusion already move and redistribute material. "
            "We therefore switched to matched controls: motorized cases are compared against no-motor controls from the same campaign, scenario, and xlink state; "
            "cross-campaign comparisons are only made for matched scenario-condition pairs."
        ),
        "",
        "## What we did",
        "- Ran single-trajectory turnover campaigns spanning baseline growth, pool-biased growth, fixed bottom chewers, diffuse bottom chewers, and motor escalation.",
        "- Used five spatial growth scenarios: full-height mixed polarity, full-height aligned, top two-thirds aligned, top one-third aligned, and top-cap aligned.",
        "- Compared no-motor, rotatable motor, fixed-global motor, no-crosslinker, and crosslinked conditions where available.",
        "- Added fair metrics: matched active fraction, excess AUC above matched controls, mass-normalized AUC, mass retention, swirl penalty, seed-zone retention, chewer-band mass, bottom-half mass, and axial spread entropy.",
        "",
        "## Main findings to say out loud",
    ]
    for row in findings:
        summary.append(f"- {row['finding']}: {row['claim']} {row['support']}")

    summary.extend(
        [
            "",
            "## Recommended show order",
            f"- Start with `figures/01_within_motor_gain_rotatable_xlink_auc.png`, then show movie `{main['run_dir']}` from `{main['campaign']}`.",
            f"- Immediately show matched no-motor control `{control['run_dir']}` and fixed-global comparator `{fixed['run_dir']}` for the same scenario.",
            (
                f"- Use `figures/07_cross_diffuse_minus_fixed_chewers_auc.png` and compare "
                f"`{topcap['campaign']}/{topcap['run_dir']}` versus `{fixed_topcap['campaign']}/{fixed_topcap['run_dir']}` "
                "to explain why diffuse bottom chewers are worth following up."
            ),
            (
                f"- Use `figures/10_cross_motor_escalation_minus_pool_bias_auc.png` and compare "
                f"`{escalation['campaign']}/{escalation['run_dir']}` versus `{pool_bias['campaign']}/{pool_bias['run_dir']}` "
                "only as an upper-bound result, because escalation lacks matched no-motor controls."
            ),
            "",
            "## Caveats",
            "- These are one-trajectory sweeps, so treat them as screening data rather than final statistics.",
            "- Top-biased scenarios can look good partly because passive growth and diffusion fill space; matched controls are therefore central to the interpretation.",
            "- Motor escalation shows large AUC increases, but active-fraction deltas are not computed there because that campaign does not include matched no-motor controls.",
            "- The bottom-chewer mechanism appears to change bottom-half mass and localized transport, but the exact causal mechanism needs the next diffuse-chewer-focused analysis and replicates.",
            "",
            "## Next steps",
            "- For the next analysis, focus on diffuse bottom chewers versus fixed bottom chewers in the top-two-thirds and top-cap scenarios.",
            "- Prioritize replicated runs for `sc02_top_two_thirds_aligned` and `sc04_top_cap_aligned` with `c01_nomotor_xlink`, `c03_rotatable_xlink`, and `c05_fixed_global_xlink`.",
            "- Add movies or image exports for the selected runs so the meeting story does not depend on opening every raw Cytosim trajectory live.",
        ]
    )

    movies_md = [
        "# Suggested Turnover Movies",
        "",
        "No pre-rendered movie files were present in the turnover folders. Use these `play` commands for live replay, or the listed render command to export a movie from Cytosim.",
        "",
    ]
    for row in movie_rows:
        movies_md.extend(
            [
                f"## {row['priority']}. {row['story']}",
                f"- Run: `{row['campaign']}` / `{row['run_dir']}` / `{row['scenario']}` / `{row['condition']}`.",
                f"- Metrics: AUC `{row['auc_vz_abs']}`, excess AUC `{row['excess_auc_vz_abs']}`, matched active fraction `{row['matched_active_transport_fraction']}`.",
                f"- Why: {row['why_show']}",
                f"- Play: `{row['play_command']}`",
                f"- Render: `{row['render_command']}`",
                "",
            ]
        )

    readme = [
        "# PI Meeting Turnover Bundle",
        "",
        "Curated package for the PI meeting on 2026-04-23.",
        "",
        "## Read First",
        "- `ONE_PAGE_SUMMARY.md`: short narrative, core findings, caveats, and next steps.",
        "- `MOVIES_TO_SHOW.md`: exact run IDs and `play` commands for the recommended trajectories.",
        "- `tables/headline_findings.csv`: machine-readable version of the headline claims.",
        "- `tables/selected_movie_runs.csv`: run paths and metrics for the movie shortlist.",
        "",
        "## Figures",
    ]
    readme.extend(f"- `{asset}`" for asset in copied_assets)
    readme.extend(
        [
            "",
            "## Source Data",
            "- Main fair-comparison source: `analysis/results/turnover_fair_comparisons`.",
            "- Per-campaign metrics: `analysis/results/turnover_growth_*`.",
            "- Raw trajectories remain in the corresponding root-level `turnover_growth_*` folders.",
        ]
    )

    (OUT / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    (OUT / "ONE_PAGE_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    (OUT / "MOVIES_TO_SHOW.md").write_text("\n".join(movies_md) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)
    (OUT / "figures").mkdir(exist_ok=True)

    effect_rows = read_csv(FAIR / "effect_deltas.csv")
    reference_rows = read_csv(FAIR / "reference_deltas.csv")
    run_rows = load_run_rows()

    findings = make_headline_findings(effect_rows, reference_rows)
    movie_rows = make_selected_movies(run_rows)
    copied_assets = copy_figures()

    write_csv(OUT / "tables" / "headline_findings.csv", findings)
    write_csv(OUT / "tables" / "selected_movie_runs.csv", movie_rows)
    write_csv(
        OUT / "tables" / "top_motor_gain_rotatable_xlink_by_auc.csv",
        top_rows(effect_rows, "delta_auc_vz_abs", effect="motor_gain_rotatable_xlink"),
    )
    write_csv(
        OUT / "tables" / "top_rotatable_minus_fixed_xlink_by_auc.csv",
        top_rows(effect_rows, "delta_auc_vz_abs", effect="orientation_gain_xlink"),
    )
    write_csv(
        OUT / "tables" / "top_diffuse_minus_fixed_bottom_chewers.csv",
        top_rows(
            reference_rows,
            "delta_auc_vz_abs",
            target_campaign="turnover_growth_pool_bias_bottom_chewers_diffuse",
            reference_campaign="turnover_growth_pool_bias_bottom_chewers",
        ),
    )
    write_csv(
        OUT / "tables" / "top_motor_escalation_minus_pool_bias.csv",
        top_rows(
            reference_rows,
            "delta_auc_vz_abs",
            target_campaign="turnover_growth_pool_bias_motor_escalation",
            reference_campaign="turnover_growth_pool_bias_sweep",
        ),
    )
    shutil.copy2(FAIR / "README.md", OUT / "tables" / "fair_comparisons_README.md")
    shutil.copy2(FAIR / "scenario_matched_thresholds.csv", OUT / "tables" / "scenario_matched_thresholds.csv")

    write_markdown(findings, movie_rows, copied_assets)
    print(OUT)


if __name__ == "__main__":
    main()
