#!/usr/bin/env bash
# Runs the JUnit tests under the TornadoVM runtime, so that the GPU paths are exercised
# (a plain `mvn test` has no TornadoVM runtime and only checks the CPU fallbacks).
# JUnit's console launcher comes from the TornadoVM SDK (JUnit 6, `execute` subcommand).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/scripts/env.sh"
(cd "$ROOT/.." && mvn -q -B -Ptornadovm -pl biojava-tornadovm -am -DskipTests install && cd "$ROOT" && mvn -q -B test-compile)
exec "$ROOT/scripts/run.sh" --add-modules=org.junit.jupiter -Dbiojava.tornado.expectGpu=true org.junit.platform.console.ConsoleLauncher execute \
  --disable-banner --details=tree --scan-classpath "$ROOT/target/test-classes" "$@"
