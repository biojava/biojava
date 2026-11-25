#!/usr/bin/env bash
# Runs a main class of this project under the TornadoVM runtime.
#   scripts/run.sh [-Dprop=value ...] <main-class> [args...]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/scripts/env.sh"
if [[ ! -f "$ROOT/target/cp.txt" || "$ROOT/pom.xml" -nt "$ROOT/target/cp.txt" ]]; then
  (cd "$ROOT" && mvn -q -B dependency:build-classpath -Dmdep.outputFile=target/cp.txt -Dmdep.includeScope=test \
     -DexcludeArtifactIds=tornado-api \
     -DexcludeGroupIds=org.junit.jupiter,org.junit.platform,org.junit.vintage,org.opentest4j,org.apiguardian >/dev/null)
fi
# JUnit is left out: the SDK puts its own JUnit on the module path, and two versions do not mix.
# The SDK's JVM options (module path, exports, native library path), with $TORNADOVM_HOME filled in.
ARGS="$ROOT/target/tornado.args"
if [[ ! -f "$ARGS" || "$TORNADOVM_HOME/tornado-argfile.template" -nt "$ARGS" ]] || ! grep -qF "$TORNADOVM_HOME/lib" "$ARGS"; then
  sed "s|\${TORNADOVM_HOME}|$TORNADOVM_HOME|g" "$TORNADOVM_HOME/tornado-argfile.template" > "$ARGS"
fi
JVM=()
while [[ $# -gt 0 && "$1" == -* ]]; do JVM+=("$1"); shift; done
CP="$ROOT/target/classes:$ROOT/target/test-classes:$(cat "$ROOT/target/cp.txt")"
MAIN="$1"; shift
exec "$JAVA_HOME/bin/java" "@$ARGS" -Xmx24g -Dtornado.device.memory=20GB ${JVM[@]+"${JVM[@]}"} -cp "$CP" "$MAIN" "$@"
