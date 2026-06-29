#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../../" && pwd)"
"$root/build_mini_turnover/bin/sim" config.cym
