#!/bin/bash
# Build the macOS app, fix code signing (needed because this repo lives in
# iCloud Drive, which stamps extended attributes that break `codesign --deep`),
# and install it to /Applications (outside iCloud, so it won't be re-tainted).
#
# Usage: scripts/build_mac.sh

set -euo pipefail

cd "$(dirname "$0")/.."

APP_NAME="MPXV Luminex QC.app"
DIST_APP="dist/$APP_NAME"
DEST_APP="/Applications/$APP_NAME"

echo "==> Running PyInstaller..."
.venv/bin/python -m PyInstaller mpox-luminex-qc.spec --noconfirm || true
# PyInstaller's own codesign step fails here because dist/ lives in iCloud
# Drive (iCloud stamps com.apple.FinderInfo on files, which codesign --deep
# rejects). That failure is expected and harmless — we re-sign properly below
# after copying out of iCloud.

echo "==> Installing to /Applications (outside iCloud, so xattrs won't be re-added)..."
rm -rf "$DEST_APP"
cp -R "$DIST_APP" "/Applications/"
xattr -cr "$DEST_APP"
codesign --force --deep --sign - "$DEST_APP"

echo "==> Verifying signature..."
codesign --verify --deep --strict "$DEST_APP"

echo "==> Done. Installed: $DEST_APP"
