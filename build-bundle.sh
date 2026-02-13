#!/bin/bash
#
# Build the self-contained skills bundle for ClawHub publishing.
# Run this before `clawhub publish bioskills-installer/`.
#

set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUNDLE_DIR="$REPO_ROOT/bioskills-installer/data"
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

echo "Building bioSkills bundle..."
echo ""

rsync -a \
    --filter=':- .gitignore' \
    --exclude='.git/' \
    --exclude='.claude/' \
    --exclude='.gitignore' \
    --exclude='bioskills-installer/' \
    --exclude='install-claude.sh' \
    --exclude='install-codex.sh' \
    --exclude='install-gemini.sh' \
    --exclude='build-bundle.sh' \
    "$REPO_ROOT/" "$TMPDIR/"

mkdir -p "$BUNDLE_DIR"
tar czf "$BUNDLE_DIR/skills-bundle.tar.gz" -C "$TMPDIR" .
shasum -a 256 "$BUNDLE_DIR/skills-bundle.tar.gz" | awk '{print $1}' > "$BUNDLE_DIR/skills-bundle.sha256"

skill_count=$(find "$TMPDIR" -mindepth 3 -name "SKILL.md" -type f | wc -l | tr -d ' ')
bundle_size=$(du -h "$BUNDLE_DIR/skills-bundle.tar.gz" | awk '{print $1}')

echo "Bundle created: $BUNDLE_DIR/skills-bundle.tar.gz"
echo "  Skills: $skill_count"
echo "  Size:   $bundle_size"
echo "  SHA256: $(cat "$BUNDLE_DIR/skills-bundle.sha256")"
