#!/bin/bash
#
# Install all bioSkills to OpenClaw from the bundled skills archive.
# No network access required - all skills are included in this package.
#
# Usage:
#   bash install-bioskills.sh                                        # Install all 425 skills
#   bash install-bioskills.sh --categories "single-cell,variant-calling"  # Selective
#   bash install-bioskills.sh --update                               # Only update changed skills
#   bash install-bioskills.sh --uninstall                            # Remove all bio-* skills

set -e

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUNDLE_PATH="$SELF_DIR/../data/skills-bundle.tar.gz"
CHECKSUM_PATH="$SELF_DIR/../data/skills-bundle.sha256"
SKILLS_DIR="$HOME/.openclaw/skills"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

for arg in "$@"; do
    if [ "$arg" = "--uninstall" ]; then
        if [ ! -d "$SKILLS_DIR" ]; then
            echo "No skills directory found at: $SKILLS_DIR"
            exit 0
        fi
        echo "Removing bioSkills from: $SKILLS_DIR"
        removed=0
        for skill_dir in "$SKILLS_DIR"/bio-*; do
            if [ -d "$skill_dir" ]; then
                rm -rf "$skill_dir"
                removed=$((removed + 1))
            fi
        done
        echo "Removed $removed skills."
        if [ -d "$HOME/.openclaw/bioskills-repo" ]; then
            rm -rf "$HOME/.openclaw/bioskills-repo"
            echo "Removed legacy cached repository."
        fi
        exit 0
    fi
done

if [ ! -f "$BUNDLE_PATH" ]; then
    echo -e "${RED}Error: Skills bundle not found at: $BUNDLE_PATH${NC}"
    exit 1
fi

if [ -f "$CHECKSUM_PATH" ]; then
    if command -v shasum &>/dev/null; then
        expected=$(cat "$CHECKSUM_PATH")
        actual=$(shasum -a 256 "$BUNDLE_PATH" | awk '{print $1}')
        if [ "$expected" != "$actual" ]; then
            echo -e "${RED}Error: Bundle checksum mismatch. File may be corrupted or tampered with.${NC}"
            echo "Expected: $expected"
            echo "Actual:   $actual"
            exit 1
        fi
    else
        echo -e "${YELLOW}Warning: shasum not found, skipping checksum verification${NC}"
    fi
fi

echo "bioSkills installer"
echo "==================="
echo ""

WORK_DIR=$(mktemp -d)
trap 'rm -rf "$WORK_DIR"' EXIT

echo "Extracting skills bundle..."
tar xzf "$BUNDLE_PATH" -C "$WORK_DIR"
echo ""

if [ -f "$WORK_DIR/install-openclaw.sh" ]; then
    bash "$WORK_DIR/install-openclaw.sh" "$@"
else
    echo -e "${RED}Error: install-openclaw.sh not found in bundle${NC}"
    exit 1
fi
