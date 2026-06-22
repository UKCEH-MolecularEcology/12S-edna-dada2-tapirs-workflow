#!/usr/bin/env bash
# Install usearch v12.0-beta1 (Linux x86_64)
# Run once before executing the workflow.
# Usage: bash bin/install_usearch.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN_DIR="$REPO_ROOT/bin"
DEST="$BIN_DIR/usearch"

if [[ -x "$DEST" ]]; then
  echo "usearch already installed at $DEST"
  "$DEST" 2>&1 | head -1
  exit 0
fi

OS="$(uname -s)"
ARCH="$(uname -m)"

case "$OS-$ARCH" in
  Linux-x86_64)  ASSET="usearch_linux_x86_12.0-beta" ;;
  Linux-aarch64) ASSET="usearch_linux_arch64_12.0-beta" ;;
  Darwin-arm64)  ASSET="usearch_osx_m_12.0-beta" ;;
  Darwin-x86_64) ASSET="usearch_osx_x86_12.0-beta" ;;
  *)
    echo "Unsupported platform: $OS-$ARCH" >&2
    exit 1 ;;
esac

URL="https://github.com/rcedgar/usearch12/releases/download/v12.0-beta1/$ASSET"

echo "Downloading usearch from $URL ..."
curl -fL -o "$DEST" "$URL"
chmod +x "$DEST"

echo "Installed: $DEST"
"$DEST" 2>&1 | head -1
