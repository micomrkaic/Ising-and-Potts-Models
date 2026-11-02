#!/bin/sh
# push_to_github.sh — put this directory on GitHub as a new repository
#
# Usage:  ./push_to_github.sh [repo-name] [public|private]
#         defaults: ising_parallel, public
#
# Needs git, and either the GitHub CLI (gh, logged in: `gh auth login`) or an
# SSH key registered with GitHub (then the repo is created via the API-free
# path: you create it empty on github.com first and this script pushes to it).
set -eu

REPO=${1:-ising_parallel}
VIS=${2:-public}
USER=${GITHUB_USER:-micomrkaic}
BRANCH=main

cd "$(dirname "$0")"

if [ ! -f ising_parallel.c ] || [ ! -f Makefile ]; then
    echo "run this from the directory containing ising_parallel.c" >&2; exit 1
fi

# ── .gitignore ─────────────────────────────────────────────────────────────
[ -f .gitignore ] || cat > .gitignore <<'EOF'
# build products
ising_parallel
ising_parallel_debug
ising_parallel_tsan
ising_parallel_prof
gmon.out
*.o
# run products
binder_L*.csv
binder.gp
ising_*.bmp
ising_*.mp4
*.tar.gz
EOF

# ── local repository ───────────────────────────────────────────────────────
if [ ! -d .git ]; then
    git init -q -b "$BRANCH" 2>/dev/null || { git init -q; git checkout -q -b "$BRANCH"; }
fi
git add -A
if git diff --cached --quiet; then
    echo "nothing new to commit"
else
    git commit -q -F - <<'EOF'
ising_parallel 2.1.1: parallel 2-D Ising/Potts Monte Carlo with SDL2

Two-colour checkerboard Metropolis over row bands (any thread count),
serial Wolff with calibrated fixed-length sweeps, q-state Potts model,
headless batch mode with jackknife errors and Binder cumulant, Onsager
exact reference values, man page, regression tests (make check) and a
Binder-crossing T_c locator (binder.sh).

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015UYsLDf32GcbVmYok4dVCm
EOF
    echo "committed"
fi

# ── remote ─────────────────────────────────────────────────────────────────
if git remote get-url origin >/dev/null 2>&1; then
    echo "remote 'origin' already set: $(git remote get-url origin)"
elif command -v gh >/dev/null 2>&1 && gh auth status >/dev/null 2>&1; then
    gh repo create "$USER/$REPO" "--$VIS" --source=. --remote=origin \
        --description "Multithreaded 2-D Ising and Potts Monte Carlo with SDL2 visualisation"
else
    echo "gh not available; using SSH remote. Create an EMPTY repo named '$REPO' at"
    echo "  https://github.com/new   (no README, no .gitignore, no licence)"
    printf "then press Enter to push, or Ctrl-C to abort: "
    read -r _
    git remote add origin "git@github.com:$USER/$REPO.git"
fi

git push -u origin "$BRANCH"
echo "pushed to https://github.com/$USER/$REPO"
