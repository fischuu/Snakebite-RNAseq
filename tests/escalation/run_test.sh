#!/usr/bin/env bash
# Test that a failing rule escalates to a bigger SLURM job. Needs Docker, and takes a few minutes:
# it builds the mini cluster, and runs the real pipeline code on it.
#
# Usage: run_test.sh
# IMAGE=<image> run_test.sh uses that image of the mini cluster instead of building one, for example
# one with another version of Snakemake.
set -euo pipefail
repo="$(cd "$(dirname "$0")/../.." && pwd)"
image=${IMAGE:-snakebite-rnaseq-minicluster}
[ -n "${IMAGE:-}" ] || docker build -q -t "$image" "$repo/tests/minicluster" > /dev/null
name=rnaseq-escalation-test
docker container rm -f "$name" > /dev/null 2>&1 || true
docker run -d --name "$name" --hostname localhost --privileged -v "$repo":/pipeline "$image" sleep infinity > /dev/null
trap 'docker container rm -f "$name" > /dev/null 2>&1' EXIT
docker exec "$name" start-cluster > /dev/null
docker exec "$name" bash /pipeline/tests/escalation/inside_cluster.sh
