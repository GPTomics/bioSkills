#!/bin/bash
#
# Install bioSkills to Claude Code
#
# Usage:
#   ./install-claude.sh              # Install globally to ~/.claude/skills/
#   ./install-claude.sh --project    # Install to current project's .claude/skills/
#   ./install-claude.sh --project /path/to/project  # Install to specific project
#   ./install-claude.sh --categories "single-cell,variant-calling"  # Selective install
#   ./install-claude.sh --validate   # Validate all skills before installing
#   ./install-claude.sh --update     # Only update changed skills
#   ./install-claude.sh --uninstall  # Remove all bio-* skills

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/install-common.sh"

TOOL_NAME="Claude Code"
DEFAULT_TARGET_DIR="$HOME/.claude/skills"
PROJECT_SUBDIR=".claude/skills"

print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Install bioSkills to Claude Code"
    echo ""
    echo "Options:"
    print_common_options
}

copy_skill_files() {
    local src_dir="$1" target_dir="$2"
    cp "$src_dir/SKILL.md" "$target_dir/SKILL.md" 2>/dev/null || return 1
    if [ -f "$src_dir/usage-guide.md" ]; then
        cp "$src_dir/usage-guide.md" "$target_dir/" 2>/dev/null || true
    fi
}

run_installer "$@"
