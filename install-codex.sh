#!/bin/bash
#
# Install bioSkills to Codex CLI
#
# Usage:
#   ./install-codex.sh              # Install to ~/.agents/skills/
#   ./install-codex.sh --project    # Install to current project's .agents/skills/
#   ./install-codex.sh --categories "single-cell,variant-calling"  # Selective install
#   ./install-codex.sh --validate   # Validate all skills before installing
#   ./install-codex.sh --update     # Only update changed skills
#   ./install-codex.sh --uninstall  # Remove all bio-* skills

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/install-common.sh"

TOOL_NAME="Codex CLI"
DEFAULT_TARGET_DIR="$HOME/.agents/skills"
PROJECT_SUBDIR=".agents/skills"

print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Install bioSkills to Codex CLI"
    echo ""
    echo "Options:"
    print_common_options
}

copy_skill_files() {
    local src_dir="$1" target_dir="$2"
    awk '
    BEGIN { in_fm=0; fm_done=0; tool_type=""; primary_tool="" }
    /^---$/ && !in_fm { in_fm=1; print; next }
    /^---$/ && in_fm {
        if (tool_type != "" || primary_tool != "") {
            print "metadata:"
            if (tool_type != "") print "  tool_type: " tool_type
            if (primary_tool != "") print "  primary_tool: " primary_tool
        }
        in_fm=0; fm_done=1; print; next
    }
    in_fm && /^tool_type:/ { tool_type=$2; next }
    in_fm && /^primary_tool:/ { sub(/^primary_tool: */, ""); primary_tool=$0; next }
    { print }
    ' "$src_dir/SKILL.md" > "$target_dir/SKILL.md" || return 1
    if [ -d "$src_dir/examples" ]; then
        rsync -a --exclude='*.pyc' --exclude='__pycache__' \
            "$src_dir/examples/" "$target_dir/scripts/"
    fi
    if [ -f "$src_dir/usage-guide.md" ]; then
        mkdir -p "$target_dir/references"
        cp "$src_dir/usage-guide.md" "$target_dir/references/"
    fi
}

run_installer "$@"
