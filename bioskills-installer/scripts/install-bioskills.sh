#!/bin/bash
#
# Install all bioSkills to OpenClaw
# Bundled with the bioskills meta-skill for ClawHub distribution
#
# Usage:
#   bash install-bioskills.sh                              # Install all 412 skills
#   bash install-bioskills.sh --categories "single-cell,variant-calling"  # Selective
#   bash install-bioskills.sh --update                     # Only update changed skills
#   bash install-bioskills.sh --uninstall                  # Remove all bio-* skills

set -e

REPO_URL="https://github.com/GPTomics/bioSkills.git"
INSTALL_DIR="$HOME/.openclaw/bioskills-repo"
SKILLS_DIR="$HOME/.openclaw/skills"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

CATEGORY_FILTER=""
UPDATE_MODE=false
UNINSTALL_MODE=false
VERBOSE=false

print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Install all bioSkills to OpenClaw"
    echo ""
    echo "Options:"
    echo "  --categories CATS     Install only specified categories (comma-separated)"
    echo "  --update              Only update skills that have changed"
    echo "  --uninstall           Remove all bio-* prefixed skills"
    echo "  --verbose             Show detailed output"
    echo "  --help                Show this help message"
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --categories)
            if [[ -n "$2" && ! "$2" =~ ^-- ]]; then
                CATEGORY_FILTER=$(echo "$2" | tr -d ' ')
                shift
            else
                echo "Error: --categories requires a comma-separated list"
                exit 1
            fi
            shift
            ;;
        --update) UPDATE_MODE=true; shift ;;
        --uninstall) UNINSTALL_MODE=true; shift ;;
        --verbose|-v) VERBOSE=true; shift ;;
        --help|-h) print_usage; exit 0 ;;
        *) echo "Unknown option: $1"; print_usage; exit 1 ;;
    esac
done

if [ "$UNINSTALL_MODE" = true ]; then
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
    [ -d "$INSTALL_DIR" ] && rm -rf "$INSTALL_DIR" && echo "Removed cached repository."
    exit 0
fi

if ! command -v git &>/dev/null; then
    echo -e "${RED}Error: git is required but not found${NC}"
    exit 1
fi

echo "bioSkills installer"
echo "==================="
echo ""

if [ -d "$INSTALL_DIR/.git" ]; then
    echo "Updating cached repository..."
    git -C "$INSTALL_DIR" pull --quiet 2>/dev/null || {
        echo "Cache stale, re-cloning..."
        tmpdir=$(mktemp -d)
        if git clone --quiet --branch main "$REPO_URL" "$tmpdir"; then
            rm -rf "$INSTALL_DIR"
            mv "$tmpdir" "$INSTALL_DIR"
        else
            rm -rf "$tmpdir"
            echo -e "${YELLOW}Warning: Could not update. Using cached version.${NC}"
        fi
    }
else
    echo "Cloning bioSkills repository..."
    mkdir -p "$(dirname "$INSTALL_DIR")"
    git clone --quiet --branch main "$REPO_URL" "$INSTALL_DIR"
fi

echo ""

EXTRA_ARGS=""
[ -n "$CATEGORY_FILTER" ] && EXTRA_ARGS="$EXTRA_ARGS --categories $CATEGORY_FILTER"
[ "$UPDATE_MODE" = true ] && EXTRA_ARGS="$EXTRA_ARGS --update"
[ "$VERBOSE" = true ] && EXTRA_ARGS="$EXTRA_ARGS --verbose"

if [ -f "$INSTALL_DIR/install-openclaw.sh" ]; then
    bash "$INSTALL_DIR/install-openclaw.sh" $EXTRA_ARGS
else
    echo -e "${RED}Error: install-openclaw.sh not found in repository${NC}"
    exit 1
fi
