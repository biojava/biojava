#!/usr/bin/env bash
# One command: fetch TornadoVM (if needed) and the test data, build, run the parity tests on the GPU and the benchmarks.
#   ./bench.sh            everything (~10 min, most of it BioJava's own CPU baselines)
#   ./bench.sh asa|sw     one benchmark
set -euo pipefail
cd "$(dirname "$0")"
source scripts/env.sh
tornado --devices | grep -E "device=|--" | head -4

mkdir -p data/pdb
for id in 1smt 4hhb 1aon 4v6x 3j3q; do
  [[ -s data/pdb/$id.cif.gz ]] || curl -sfL -o data/pdb/$id.cif.gz "https://files.rcsb.org/download/$id.cif.gz"
done
cp -n src/test/resources/PF00104_small.fasta data/ 2>/dev/null || true
(cd .. && mvn -q -B -Ptornadovm -pl biojava-tornadovm -am -DskipTests install)

what="${1:-all}"
if [[ "$what" == all ]]; then
  echo "== parity tests on the GPU"; scripts/test-gpu.sh 2>/dev/null | grep -E "tests (successful|failed)"
fi
if [[ "$what" == all || "$what" == asa ]]; then
  echo "== ASA (Shrake-Rupley, 1000 sphere points)"
  scripts/run.sh org.biojava.nbio.tornadovm.bench.AsaBench 1SMT 4HHB 1AON 4V6X 3J3Q 2>/dev/null | grep -vE "cold|DEBUG|INFO|WARN"
fi
if [[ "$what" == all || "$what" == sw ]]; then
  echo "== All-pairs alignment scores (BLOSUM62, gaps 10/1)"
  scripts/run.sh org.biojava.nbio.tornadovm.bench.SwBench random:40:100:300 random:200:100:500 2>/dev/null | grep -vE "DEBUG|INFO|WARN"
  scripts/run.sh -Dbench.skipcpu=true org.biojava.nbio.tornadovm.bench.SwBench fasta:data/PF00104_small.fasta random:1000:200:500 2>/dev/null | grep -vE "DEBUG|INFO|WARN|threads in"
fi
