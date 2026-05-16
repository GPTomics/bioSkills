#!/bin/bash
#
# Install bioSkills to OpenCode
#
# OpenCode auto-discovers Agent Skills from .opencode/skills/, .claude/skills/,
# and .agents/skills/ (project and global). This installer targets the native
# ~/.config/opencode/skills/ path; users who already ran install-claude.sh or
# install-codex.sh do not need to run this installer.
#
# Usage:
#   ./install-opencode.sh              # Install to ~/.config/opencode/skills/
#   ./install-opencode.sh --project    # Install to current project's .opencode/skills/
#   ./install-opencode.sh --categories "single-cell,variant-calling"  # Selective install
#   ./install-opencode.sh --validate   # Validate all skills before installing
#   ./install-opencode.sh --update     # Only update changed skills
#   ./install-opencode.sh --uninstall  # Remove all bio-* skills

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/install-common.sh"

TOOL_NAME="OpenCode"
DEFAULT_TARGET_DIR="$HOME/.config/opencode/skills"
PROJECT_SUBDIR=".opencode/skills"

print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Install bioSkills to OpenCode"
    echo ""
    echo "Options:"
    print_common_options
}

copy_skill_files() {
    local src_dir="$1" target_dir="$2"
    local full_name=$(basename "$target_dir")
    awk -v new_name="$full_name" '
    BEGIN { in_fm=0; fm_done=0; tool_type=""; primary_tool=""; name_written=0 }
    /^---$/ && !in_fm { in_fm=1; print; next }
    /^---$/ && in_fm {
        if (!name_written) { print "name: " new_name; name_written=1 }
        if (tool_type != "" || primary_tool != "") {
            print "metadata:"
            if (tool_type != "") print "  tool_type: " tool_type
            if (primary_tool != "") print "  primary_tool: " primary_tool
        }
        in_fm=0; fm_done=1; print; next
    }
    in_fm && /^name:/ { print "name: " new_name; name_written=1; next }
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
