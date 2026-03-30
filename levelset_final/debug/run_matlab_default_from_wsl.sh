#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MATLAB_WIN="/mnt/d/Program Files/matlab/R2024a/bin/matlab.exe"
MATLAB_CD_WIN="$(wslpath -w "$PROJECT_ROOT")"

"$MATLAB_WIN" -batch "cd('$MATLAB_CD_WIN'); addpath(genpath(pwd), '-begin'); results = fiber_levelset('default');"
