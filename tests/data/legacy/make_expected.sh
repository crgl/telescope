#!/bin/bash
# Regenerate expected/ from the reference implementation: Python Telescope 1.0.4.1
# (hanalysis fork, this repository's telescope/ directory) for the py310 and py37 cases,
# and this program itself for the cases marked "rust" (which have no Python equivalent).
#
#   TESTBED=~/herv_working_data/telescope_ccle bash make_expected.sh
#
# TESTBED must hold two environments with Telescope installed: env/ (Python 3.10,
# numpy 1.26.4, scipy 1.15.2) and env_py37/ (Python 3.7, numpy 1.21.6, scipy 1.7.3).
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd); repo=$(cd "$here/../../.." && pwd)
T=${TESTBED:-$HOME/herv_working_data/telescope_ccle}
work=$(mktemp -d); mkdir -p "$here/expected"

# One line per record of an updated BAM: the fields --updated_sam is responsible for.
updated_table() {
  "$T/env/bin/samtools" view "$1" | awk -F'\t' 'BEGIN{OFS="\t"} {
    zf=zt=zb=xp=yc=".";
    for (i=12; i<=NF; i++) { t=substr($i,1,2); v=substr($i,6);
      if (t=="ZF") zf=v; else if (t=="ZT") zt=v; else if (t=="ZB") zb=v; else if (t=="XP") xp=v; else if (t=="YC") yc=v }
    print $1,$2,$3,$4,$5,zf,zt,zb,xp,yc }'
}

grep -v '^#' "$here/cases.tsv" | while IFS=$'\t' read -r name ref aln gtf args upd; do
  [ -z "$name" ] && continue
  [ "$args" = "-" ] && args=""
  out="$work/$name"; mkdir -p "$out"
  flag=""; [ "$upd" = yes ] && flag="--updated_sam"
  case $ref in
    py310) (cd "$out" && "$T/env/bin/telescope" assign $flag --outdir "$out" --exp_tag t $args "$here/$aln" "$here/$gtf" 2> log) ;;
    py37)  # the Python run is plain Telescope; --tie_hash only tells the Rust side which Python to match
           (cd "$out" && "$T/env_py37/bin/telescope" assign $flag --outdir "$out" --exp_tag t ${args/--tie_hash python37/} "$here/$aln" "$here/$gtf" 2> log) ;;
    rust)  "$repo/target/release/telescope_rs" assign --quiet $flag --outdir "$out" --exp_tag t $args "$here/$aln" "$here/$gtf" 2> "$out/log" ;;
  esac
  if [ ! -s "$out/t-telescope_report.tsv" ]; then echo "FAILED $name"; tail -2 "$out/log"; continue; fi
  cp "$out/t-telescope_report.tsv" "$here/expected/$name.report.tsv"
  [ "$upd" = yes ] && updated_table "$out/t-updated.bam" > "$here/expected/$name.updated.tsv"
  echo "ok $name ($ref): $(($(wc -l < "$here/expected/$name.report.tsv") - 2)) features"
done
