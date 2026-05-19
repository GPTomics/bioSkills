#!/bin/bash
#
# Install bioSkills to Antigravity CLI
#
# Antigravity CLI is Google's replacement for Gemini CLI (Gemini CLI sunset
# 2026-06-18). It adopts Anthropic's open Agent Skills standard.
#
# Skills layout:
#   Global:    ~/.gemini/antigravity/skills/<skill-name>/SKILL.md
#   Workspace: <project>/.agents/skills/<skill-name>/SKILL.md
#
# Usage:
#   ./install-antigravity.sh              # Install to ~/.gemini/antigravity/skills/
#   ./install-antigravity.sh --project    # Install to current project's .agents/skills/
#   ./install-antigravity.sh --categories "single-cell,variant-calling"  # Selective install
#   ./install-antigravity.sh --validate   # Validate all skills before installing
#   ./install-antigravity.sh --update     # Only update changed skills
#   ./install-antigravity.sh --uninstall  # Remove all bio-* skills

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/install-common.sh"

TOOL_NAME="Antigravity CLI"
DEFAULT_TARGET_DIR="$HOME/.gemini/antigravity/skills"
PROJECT_SUBDIR=".agents/skills"

print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Install bioSkills to Antigravity CLI"
    echo ""
    echo "Options:"
    print_common_options
}

copy_skill_files() {
    local src_dir="$1" target_dir="$2"
    cp "$src_dir/SKILL.md" "$target_dir/SKILL.md" || return 1
    [ -f "$src_dir/usage-guide.md" ] && \
        cp "$src_dir/usage-guide.md" "$target_dir/"
}

run_installer "$@"
