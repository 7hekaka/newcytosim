# Codex Sync Protocol

Last updated: 2026-06-28

This file is the cross-machine handoff note for Codex sessions. A new Codex instance should read this file and `PROJECT_PLANNER.md` before doing project work.

## Canonical Workspace

- Active repo: `/home/thekaka/project/cytosim`
- Windows view of the same repo: `\\wsl.localhost\Ubuntu-22.04\home\thekaka\project\cytosim`
- Avoid `C:\home\thekaka\project\cytosim`; it has been a stale/bad mirror.
- GitHub remote: `git@github.com:7hekaka/newcytosim.git`
- Upstream Cytosim remote: `https://gitlab.com/f-nedelec/cytosim`

## What Git Should Sync

Sync these through Git:

- Cytosim source edits.
- Generator scripts.
- Analysis scripts.
- Campaign manifests, README files, and selected `config.cym` files.
- Planner, protocol, and manuscript/figure source text.
- Small final figures only when they are intentionally curated.

Keep these out of Git unless deliberately curated:

- Raw simulation job folders.
- `objects.cmo`, `properties.cmp`, `messages.cmo`, logs, frame dumps, movies, and large returned analysis folders.
- Machine-specific build folders and compiled binaries.

Large returned job folders should stay on the desktop/work PC, cluster scratch, or external storage. The Mac/laptop should normally hold only code, configs, scripts, manifests, and selected small outputs.

## Start Of Session Checklist

Run this before starting real work:

```bash
cd /home/thekaka/project/cytosim
git fetch --all --prune
git status --short --branch
sed -n '1,180p' CODEX_SYNC.md
sed -n '1,220p' PROJECT_PLANNER.md
```

If the branch is dirty, divergent, or behind remote, do not run broad `git pull` or `git add .`. First identify whether the local changes are part of the current task.

## End Of Session Checklist

Before ending a work session:

1. Update `PROJECT_PLANNER.md` when priorities, completed tasks, data locations, or next actions changed.
2. Update the current handoff snapshot below.
3. Stage only intentional files.
4. Commit with a specific message.
5. Push the branch.

Useful checks:

```bash
git status --short --branch
git diff --name-only
git diff --cached --name-only
```

Avoid `git add .` in this repo because many output folders are large and easy to stage accidentally.

## Current Handoff Snapshot

Date: 2026-06-28

Current major project state:

- Paper/manuscript branch: metric-focused Results/Discussion draft and figures passed PI/user review. User is moving remaining edits to Overleaf, so local paper work is mostly paused.
- Central-cell F-actin turnover branch: full comparable continuous-supply simulations are currently running on the cluster. This is the full-size central-cell-like campaign with aggressive bottom minus-end depolymerization and top actin pulse supply.
- Current running campaign: `turnover_continuous_supply_full_comparable`.
- Expected cluster location: `/lustre/isaac24/scratch/kacheamp/newcytosim/turnover_continuous_supply_full_comparable`.
- Campaign comparison: motor-free xlink controls versus rotatable motor-cluster cases, with 10 replicates each.
- Key model ingredients: full annulus scale, top actin supply every 60 s, aggressive bottom chewers, 80 motor clusters x 12 motors per cluster for motor cases.
- Heavy returned outputs should not be expected on laptop/Mac unless explicitly copied there.

Current Git caution:

- As of this snapshot, local `main` was ahead of `origin/main` by 2 commits and behind by 6 commits, with many staged/untracked files. Clean sync should happen through deliberate curated commits or a clean sync branch, not broad pull/push from the dirty tree.

Next likely actions:

- After full comparable turnover jobs return, analyze motor versus no-motor cases using post-motor/intended analysis windows.
- Recreate collaborator-style cyan snapshots and smoother movies for the full comparable campaign.
- Update presentation slides with full-size results after the analysis is stable.

## Laptop Or New Work Computer Setup

On the Mac laptop:

```bash
mkdir -p ~/project
cd ~/project
git clone git@github.com:7hekaka/newcytosim.git cytosim
cd cytosim
git fetch --all --prune
git status --short --branch
```

On a Windows work computer, use WSL2 Ubuntu and clone inside WSL:

```bash
mkdir -p ~/project
cd ~/project
git clone git@github.com:7hekaka/newcytosim.git cytosim
cd cytosim
git fetch --all --prune
git status --short --branch
```

Then tell the local Codex:

```text
Read CODEX_SYNC.md and PROJECT_PLANNER.md first. Use the cloned repo path on this machine as the active workspace. Do not assume heavy job folders are present locally.
```

## Logging Rule

Every meaningful work block should leave enough trace for the next machine to continue:

- What changed.
- Why it changed.
- Which files/folders matter.
- Which heavy outputs are local-only.
- What the next concrete action is.

Use `PROJECT_PLANNER.md` for project tasks and this file for cross-machine handoff/state.
