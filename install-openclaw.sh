#!/bin/bash
#
# Install bioSkills to OpenClaw
#
# Usage:
#   ./install-openclaw.sh                                          # Install all skills globally
#   ./install-openclaw.sh --categories "single-cell,variant-calling"  # Selective install
#   ./install-openclaw.sh --project /path/to/workspace             # Install to workspace
#   ./install-openclaw.sh --tool-type-metadata                     # Add OpenClaw dependency metadata
#   ./install-openclaw.sh --dry-run                                # Preview install + token estimate
#   ./install-openclaw.sh --validate                               # Validate all skills
#   ./install-openclaw.sh --update                                 # Only update changed skills
#   ./install-openclaw.sh --uninstall                              # Remove all bio-* skills

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/install-common.sh"

TOOL_NAME="OpenClaw"
DEFAULT_TARGET_DIR="$HOME/.openclaw/skills"
PROJECT_SUBDIR="skills"

TOOL_TYPE_METADATA=false

print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Install bioSkills to OpenClaw"
    echo ""
    echo "Options:"
    print_common_options
    echo ""
    echo "OpenClaw-specific options:"
    echo "  --tool-type-metadata  Add OpenClaw dependency metadata to installed skills"
}

parse_extra_arg() {
    case "$1" in
        --tool-type-metadata) TOOL_TYPE_METADATA=true; EXTRA_SHIFT=1; return 0 ;;
        *) return 1 ;;
    esac
}

inject_metadata() {
    local src="$1" dst="$2" tool_type="$3"
    local meta=""
    case "$tool_type" in
        python) meta='metadata: {"openclaw":{"requires":{"bins":["python3"]},"os":["darwin","linux"]}}' ;;
        r)      meta='metadata: {"openclaw":{"requires":{"bins":["Rscript"]},"os":["darwin","linux"]}}' ;;
        mixed)  meta='metadata: {"openclaw":{"requires":{"anyBins":["python3","Rscript"]},"os":["darwin","linux"]}}' ;;
        *)      cp "$src" "$dst"; return ;;
    esac
    awk -v meta="$meta" 'BEGIN{n=0} /^---$/{n++; if(n==2){print meta}} {print}' "$src" > "$dst"
}

copy_skill_files() {
    local src_dir="$1" target_dir="$2"
    if [ "$TOOL_TYPE_METADATA" = true ]; then
        local tool_type=$(grep "^tool_type:" "$src_dir/SKILL.md" | head -1 | awk -F': ' '{print $2}')
        inject_metadata "$src_dir/SKILL.md" "$target_dir/SKILL.md" "$tool_type" || return 1
    else
        cp "$src_dir/SKILL.md" "$target_dir/SKILL.md" || return 1
    fi
    if [ -d "$src_dir/examples" ]; then
        rsync -a --exclude='*.pyc' --exclude='__pycache__' \
            "$src_dir/examples/" "$target_dir/examples/"
    fi
    if [ -f "$src_dir/usage-guide.md" ]; then
        cp "$src_dir/usage-guide.md" "$target_dir/"
    fi
}

check_project_conflicts() {
    local target_dir="$1"
    [ "$INSTALL_MODE" != "project" ] && return 0
    if [ -d "$target_dir" ]; then
        local non_bio=$(find "$target_dir" -maxdepth 1 -mindepth 1 -type d ! -name "bio-*" | head -1)
        if [ -n "$non_bio" ] && [ "$FORCE" != true ]; then
            echo -e "${YELLOW}Warning: $target_dir already contains non-bioSkills content.${NC}"
            echo "Use --force to install anyway."
            exit 1
        fi
    fi
    return 0
}

post_install_message() {
    echo "Skills are now available in OpenClaw."
    echo "Start a new OpenClaw session for skills to take effect."
}

run_installer "$@"
