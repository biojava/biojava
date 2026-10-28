# biojava-tornadovm

GPU paths for two BioJava hot spots, with [TornadoVM](https://github.com/beehive-lab/TornadoVM). The
classes are drop-in replacements and return **exactly** the same results as the BioJava code they
replace:

| BioJava | Drop-in replacement | Result |
|---|---|---|
| `org.biojava.nbio.structure.asa.AsaCalculator` | `org.biojava.nbio.tornadovm.asa.TornadoAsaCalculator` | bit-identical per-atom and per-residue ASA |
| `Alignments.getAllPairsScores` / `getAllPairsScorers` (`GLOBAL`, `LOCAL`) | `org.biojava.nbio.tornadovm.align.TornadoAlignments` | identical scores and scorers |

Without a TornadoVM runtime, below a size threshold, or after a GPU failure, the classes run an exact
multi-threaded CPU path instead.

This module needs **JDK 21** and is not part of the default build: it builds only with the `tornadovm`
profile, so the rest of BioJava keeps its JDK 11 baseline.

## Run it

Linux with an NVIDIA GPU, using [SDKMAN](https://sdkman.io). For AMD or Intel GPUs, install
`tornadovm 8.0.0-jdk21-opencl` and set `TORNADO_BACKEND=opencl`.

```bash
# 1. JDK 21 and TornadoVM 8.0.0
curl -s "https://get.sdkman.io" | bash && source ~/.sdkman/bin/sdkman-init.sh   # if you have no SDKMAN yet
sdk install java 21.0.2-open
sdk install tornadovm 8.0.0-jdk21-cuda      # answering "n" to "set as default" is fine
export JAVA_HOME=$HOME/.sdkman/candidates/java/21.0.2-open

# 2. Build BioJava with this module
git clone https://github.com/biojava/biojava.git && cd biojava
mvn -Ptornadovm -pl biojava-tornadovm -am install -DskipTests

# 3. The tests on the GPU, then the benchmarks
cd biojava-tornadovm
scripts/test-gpu.sh          # JUnit under TornadoVM: GPU results == BioJava results (23 tests)
./bench.sh                   # parity tests + ASA and all-pairs benchmarks (downloads 5 PDB entries)
./bench.sh asa               # or one benchmark: asa | sw
```

The scripts find `tornadovm 8.0.0-jdk21-<backend>` in SDKMAN by its version, whatever your default
is. If it is missing, they install it with `sdk install`, without making it the default. Without
SDKMAN, they download the SDK from the
[TornadoVM release](https://github.com/beehive-lab/TornadoVM/releases/tag/v8.0.0) into `.tornado/`.
`TORNADO_SDK=/path/to/sdk` uses a given SDK.

`mvn -Ptornadovm -pl biojava-tornadovm test` runs the same tests on a plain JVM, which exercises the
CPU paths.

## Use it in an application

```java
// before: new AsaCalculator(structure, AsaCalculator.DEFAULT_PROBE_SIZE, 1000, nThreads, false)
TornadoAsaCalculator asa = new TornadoAsaCalculator(structure, AsaCalculator.DEFAULT_PROBE_SIZE,
        AsaCalculator.DEFAULT_N_SPHERE_POINTS, nThreads, false);
double[] perAtom = asa.calculateAsas();

// before: Alignments.getAllPairsScores(seqs, PairwiseSequenceScorerType.GLOBAL, gaps, matrix)
double[] scores = TornadoAlignments.getAllPairsScores(seqs, PairwiseSequenceScorerType.GLOBAL, gaps, matrix);
```

Start the application on JDK 21 with the TornadoVM SDK's JVM options. Without them, it runs on the
CPU paths:

```bash
export TORNADOVM_HOME=$HOME/.sdkman/candidates/tornadovm/8.0.0-jdk21-cuda
sed "s|\${TORNADOVM_HOME}|$TORNADOVM_HOME|g" $TORNADOVM_HOME/tornado-argfile.template > tornado.args
$JAVA_HOME/bin/java @tornado.args -cp my-app.jar:<dependencies> com.example.Main
```

`tornado-api` is a `provided` dependency: at run time it comes from the SDK, so do not bundle it.

- `-Dbiojava.tornado=off` forces the CPU paths; `-Dbiojava.tornado=force` ignores the size thresholds
  (GPU from 1,000 atoms for ASA, and from 1e7 DP cells in total for alignment).
- `-Dasa.fused.device=<driver>:<device>` and `-Dsw.scores.device=<driver>:<device>` choose the device;
  `tornado --devices` lists them.
- GPUs without FP64 (e.g. Intel iGPUs) are detected. ASA then uses a single-precision kernel and
  recomputes near-boundary cases exactly on the CPU, so results stay identical.
- `*_IDENTITIES`, `*_SIMILARITIES`, `KMERS` and `WU_MANBER` scorers are delegated to BioJava.

## Results

TornadoVM 8.0.0, CUDA, RTX 4090; i9-13900K with 32 threads. End to end, including transfers. All
results are identical to BioJava's. "Lean Java" is the GPU kernel's algorithm in plain Java on all
cores, which is also this module's CPU path.

| Workload | BioJava | lean Java (32 thr) | GPU | vs BioJava | vs lean Java |
|---|---:|---:|---:|---:|---:|
| ASA 4V6X, 238k atoms | 558 ms (32 thr) | 109 ms | 35.7 ms | 15.6x | 3.0x |
| ASA 3J3Q, 2.4M atoms | 6.92 s (32 thr) | 1.16 s | 373 ms | 18.5x | 3.1x |
| All-pairs global, 200 proteins, 20k pairs | 9.97 s | 417 ms | 17.0 ms | 586x | 24.5x |
| All-pairs local, 200 proteins, 20k pairs | 10.5 s | 293 ms | 17.0 ms | 615x | 17.2x |
| All-pairs global, 1,000 proteins, 500k pairs | (skipped) | 15.0 s | 265 ms | | 56.6x |

Tested on an RTX 4090 with the CUDA and OpenCL backends, and on an Intel UHD 770 with OpenCL. The 23
GPU parity tests pass on all three.
