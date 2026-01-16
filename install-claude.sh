#!/bin/bash
#
# Install bioSkills to Claude Code
#
# Usage:
#   ./install.sh              # Install globally to ~/.claude/skills/
#   ./install.sh --project    # Install to current project's .claude/skills/
#   ./install.sh --project /path/to/project  # Install to specific project

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALL_MODE="global"
PROJECT_PATH=""

print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --global              Install to ~/.claude/skills/ (default)"
    echo "  --project [PATH]      Install to .claude/skills/ in current or specified directory"
    echo "  --list                List available skills"
    echo "  --help                Show this help message"
}

list_skills() {
    echo "Available bioSkills:"
    echo ""
    find "$SCRIPT_DIR" -name "SKILL.md" -type f | while read -r skill_file; do
        skill_dir=$(dirname "$skill_file")
        skill_name=$(basename "$skill_dir")
        category=$(basename "$(dirname "$skill_dir")")
        description=$(grep -A1 "^description:" "$skill_file" | head -1 | sed 's/description: //')
        echo "  $category/$skill_name"
        echo "    $description" | cut -c1-80
        echo ""
    done
}

install_skills() {
    local target_dir="$1"

    echo "Installing bioSkills to: $target_dir"
    echo ""

    mkdir -p "$target_dir"

    find "$SCRIPT_DIR" -name "SKILL.md" -type f | while read -r skill_file; do
        skill_dir=$(dirname "$skill_file")
        skill_name=$(basename "$skill_dir")
        category=$(basename "$(dirname "$skill_dir")")
        full_skill_name="bio-${category}-${skill_name}"

        target_skill_dir="$target_dir/$full_skill_name"
        mkdir -p "$target_skill_dir"

        cp "$skill_file" "$target_skill_dir/SKILL.md"

        if [ -f "$skill_dir/usage-guide.md" ]; then
            cp "$skill_dir/usage-guide.md" "$target_skill_dir/"
        fi

        echo "  Installed: $full_skill_name"
    done

    echo ""
    echo "Installation complete."
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --global)
            INSTALL_MODE="global"
            shift
            ;;
        --project)
            INSTALL_MODE="project"
            if [[ -n "$2" && ! "$2" =~ ^-- ]]; then
                PROJECT_PATH="$2"
                shift
            fi
            shift
            ;;
        --list)
            list_skills
            exit 0
            ;;
        --help|-h)
            print_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            print_usage
            exit 1
            ;;
    esac
done

if [ "$INSTALL_MODE" = "global" ]; then
    TARGET_DIR="$HOME/.claude/skills"
else
    if [ -n "$PROJECT_PATH" ]; then
        TARGET_DIR="$PROJECT_PATH/.claude/skills"
    else
        TARGET_DIR="$(pwd)/.claude/skills"
    fi
fi

install_skills "$TARGET_DIR"

echo ""
echo "Skills are now available in Claude Code."
echo "They will be auto-invoked when relevant to your bioinformatics tasks."
