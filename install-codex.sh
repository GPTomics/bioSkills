#!/bin/bash
#
# Install bioSkills to Codex CLI
#
# Usage:
#   ./install-codex.sh              # Install to ~/.codex/skills/
#   ./install-codex.sh --project    # Install to current project's .codex/skills/
#   ./install-codex.sh --list       # List available skills

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALL_MODE="global"
PROJECT_PATH=""

print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --global              Install to ~/.codex/skills/ (default)"
    echo "  --project [PATH]      Install to .codex/skills/ in current or specified directory"
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

        # Copy SKILL.md
        cp "$skill_file" "$target_skill_dir/SKILL.md"

        # Convert examples/ to scripts/ (Codex convention)
        if [ -d "$skill_dir/examples" ]; then
            mkdir -p "$target_skill_dir/scripts"
            cp -r "$skill_dir/examples/"* "$target_skill_dir/scripts/" 2>/dev/null || true
        fi

        # Convert usage-guide.md to references/ (Codex convention)
        if [ -f "$skill_dir/usage-guide.md" ]; then
            mkdir -p "$target_skill_dir/references"
            cp "$skill_dir/usage-guide.md" "$target_skill_dir/references/"
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
    TARGET_DIR="$HOME/.codex/skills"
else
    if [ -n "$PROJECT_PATH" ]; then
        TARGET_DIR="$PROJECT_PATH/.codex/skills"
    else
        TARGET_DIR="$(pwd)/.codex/skills"
    fi
fi

install_skills "$TARGET_DIR"

echo ""
echo "Skills are now available in Codex CLI."
echo "They will be auto-invoked when relevant to your bioinformatics tasks."
