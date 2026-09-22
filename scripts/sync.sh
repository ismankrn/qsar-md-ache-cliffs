#!/usr/bin/env bash
# Sync the project between this laptop and the Hetzner server.
#
# Run this FROM THE LAPTOP. It uses the `qsar` SSH alias, which carries the
# port (443, because 22 is blocked on campus) and agent forwarding.
#
#   ./sync.sh pull     server -> laptop   (fetch results)
#   ./sync.sh push     laptop -> server   (send code/data changes)
#   ./sync.sh status   show what differs, change nothing
#
# Git handles scripts/ and results/ with history; this handles the bulk
# directories that are awkward or wasteful to version.
set -euo pipefail

REMOTE=qsar
RPATH=/root/qsar-md
LPATH="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

EXCLUDES=(
  --exclude '.git'
  --exclude '__pycache__'
  --exclude '*.pyc'
  --exclude '.DS_Store'
  --exclude 'env_create.log'
)

case "${1:-}" in
  pull)
    echo "server -> laptop"
    rsync -avh --progress "${EXCLUDES[@]}" \
      "$REMOTE:$RPATH/" "$LPATH/"
    ;;
  push)
    echo "laptop -> server"
    rsync -avh --progress "${EXCLUDES[@]}" \
      "$LPATH/" "$REMOTE:$RPATH/"
    ;;
  status)
    echo "differences (dry run, nothing is changed)"
    echo "--- would be pulled from the server ---"
    rsync -an --itemize-changes "${EXCLUDES[@]}" \
      "$REMOTE:$RPATH/" "$LPATH/" | grep -E '^[<>ch]' || echo "  (none)"
    echo "--- would be pushed to the server ---"
    rsync -an --itemize-changes "${EXCLUDES[@]}" \
      "$LPATH/" "$REMOTE:$RPATH/" | grep -E '^[<>ch]' || echo "  (none)"
    ;;
  *)
    echo "usage: $0 {pull|push|status}" >&2; exit 1
    ;;
esac
