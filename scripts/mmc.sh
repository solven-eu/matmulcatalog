#!/usr/bin/env bash
#
# mmc.sh — Matrix-Multiplication-Catalog operations runner.
#
# One entry point for the recurring operations on this repo, grouped into the
# three categories we actually think in:
#
#   compute    Slow optimisation work, focused on a band / single shape:
#                ranks  → search for new / better schemes   (SchemeSweep)
#                slp    → minimise additive complexity (SLP) (MaterialiseAdditionsSlp)
#              These can run for minutes→hours; launched with a big heap and
#              heap-dump-on-OOM, per CLAUDE.md.
#
#   index      Regenerate the docs manifests the GitHub-Pages browser reads:
#                derived/cited bounds, then catalog.json. Fast.
#
#   sanitize   Cheap, deterministic, idempotent canonicalisation of scheme JSONs:
#                format → MatrixJsonFormatter (rows inline, matrices vertical)
#                zt     → compute & stamp the `zt` ternary-integer flag (MaterialiseZT)
#                all    → zt then format (zt already reformats what it touches)
#
# Usage:
#   scripts/mmc.sh                      # interactive menu
#   scripts/mmc.sh sanitize             # interactive sanitize sub-menu
#   scripts/mmc.sh sanitize all [path]
#   scripts/mmc.sh sanitize zt [path]
#   scripts/mmc.sh sanitize format [path]
#   scripts/mmc.sh index
#   scripts/mmc.sh compute ranks  [SchemeSweep args…]
#   scripts/mmc.sh compute slp    [MaterialiseAdditionsSlp args…]
#
set -euo pipefail

# Repo root = parent of this script's dir.
cd "$(dirname "$0")/.."
ROOT="$(pwd)"

CP_CACHE="target/mmc-classpath.txt"

# ── classpath ───────────────────────────────────────────────────────────────
# Build (if needed) and cache the runtime classpath. Rebuilt when pom.xml is
# newer than the cache, or target/classes is missing.
ensure_classpath() {
  if [[ ! -d target/classes || ! -f target/test-classes/.compiled-marker || pom.xml -nt "$CP_CACHE" ]]; then
    echo "▸ compiling (mvn test-compile)…" >&2
    mvn -q -o -ntp test-compile
    : > target/test-classes/.compiled-marker
  fi
  if [[ ! -f "$CP_CACHE" || pom.xml -nt "$CP_CACHE" ]]; then
    echo "▸ resolving classpath…" >&2
    mvn -q -o dependency:build-classpath -Dmdep.outputFile="$CP_CACHE" >/dev/null
  fi
  CP="target/classes:target/test-classes:$(cat "$CP_CACHE")"
}

# run_java <heap> <mainclass> [args…]
run_java() {
  local heap="$1"; shift
  local main="$1"; shift
  ensure_classpath
  mkdir -p target/oom-dumps
  echo "▸ java -Xmx${heap} ${main} $*" >&2
  MAVEN_OPTS="" java -Xmx"${heap}" \
      -XX:+HeapDumpOnOutOfMemoryError -XX:HeapDumpPath=target/oom-dumps/ \
      -cp "$CP" "$main" "$@"
}

# ── operations ────────────────────────────────────────────────────────────────
op_sanitize() {
  local sub="${1:-}"; [[ $# -gt 0 ]] && shift || true
  case "$sub" in
    zt)      run_java 2g eu.solven.matmul.docs.migrate.MaterialiseZT "$@" ;;
    format)  run_java 2g eu.solven.matmul.catalog.ReformatSchemes "$@" ;;
    all|"")
      # zt first (it reformats the Z files it stamps), then format the rest.
      run_java 2g eu.solven.matmul.docs.migrate.MaterialiseZT "$@"
      run_java 2g eu.solven.matmul.catalog.ReformatSchemes "$@"
      ;;
    *) echo "unknown sanitize op: $sub (expected: zt | format | all)" >&2; exit 1 ;;
  esac
}

op_index() {
  # Bounds first, then the consumer-facing manifest.
  run_java 2g eu.solven.matmul.docs.generate.GenerateDerivedBounds
  run_java 2g eu.solven.matmul.docs.generate.GenerateCitedBounds
  run_java 2g eu.solven.matmul.docs.generate.GenerateCatalogManifest
}

op_compute() {
  local sub="${1:-}"; [[ $# -gt 0 ]] && shift || true
  case "$sub" in
    ranks) run_java 4g eu.solven.matmul.docs.SchemeSweep "$@" ;;
    slp)   run_java 4g eu.solven.matmul.docs.migrate.MaterialiseAdditionsSlp "$@" ;;
    *) echo "unknown compute op: $sub (expected: ranks | slp [args…])" >&2; exit 1 ;;
  esac
}

# ── interactive menus ─────────────────────────────────────────────────────────
menu_sanitize() {
  echo "Sanitize — cheap, idempotent canonicalisation of scheme JSONs:"
  local choice
  select choice in "all (zt + format)" "zt only" "format only" "cancel"; do
    case "$choice" in
      "all (zt + format)") op_sanitize all; break ;;
      "zt only")           op_sanitize zt; break ;;
      "format only")       op_sanitize format; break ;;
      cancel|"")           echo "cancelled"; break ;;
    esac
  done
}

menu_compute() {
  echo "Compute — SLOW; typically scoped to a band / single shape (pass args after choosing):"
  local choice
  select choice in "evaluate new ranks (SchemeSweep)" "evaluate SLP / additions (MaterialiseAdditionsSlp)" "cancel"; do
    case "$choice" in
      "evaluate new ranks (SchemeSweep)")
        echo "Re-run non-interactively with args, e.g.:  scripts/mmc.sh compute ranks --mode=materialize --cubic=2-32" >&2
        op_compute ranks; break ;;
      "evaluate SLP / additions (MaterialiseAdditionsSlp)")
        echo "Re-run non-interactively with args, e.g.:  scripts/mmc.sh compute slp <band/shape args>" >&2
        op_compute slp; break ;;
      cancel|"") echo "cancelled"; break ;;
    esac
  done
}

menu_top() {
  echo "What do you want to run?"
  local choice
  select choice in \
      "compute  — evaluate new ranks / SLP (slow, band-focused)" \
      "index    — regenerate docs manifests (catalog.json + bounds)" \
      "sanitize — format JSON + compute ZT (fast, idempotent)" \
      "quit"; do
    case "$choice" in
      compute*)  menu_compute; break ;;
      index*)    op_index; break ;;
      sanitize*) menu_sanitize; break ;;
      quit|"")   echo "bye"; break ;;
    esac
  done
}

# ── dispatch ──────────────────────────────────────────────────────────────────
case "${1:-}" in
  "")        menu_top ;;
  compute)   shift; [[ $# -gt 0 ]] && op_compute "$@" || menu_compute ;;
  index)     op_index ;;
  sanitize)  shift; [[ $# -gt 0 ]] && op_sanitize "$@" || menu_sanitize ;;
  -h|--help|help)
    sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//' ;;
  *) echo "unknown operation: $1 (expected: compute | index | sanitize)" >&2; exit 1 ;;
esac
