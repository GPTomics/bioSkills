#!/bin/bash
#
# Install bioSkills to Codex CLI
#
# Usage:
#   ./install-codex.sh              # Install to ~/.codex/skills/
#   ./install-codex.sh --project    # Install to current project's .codex/skills/
#   ./install-codex.sh --validate   # Validate all skills before installing
#   ./install-codex.sh --update     # Only update changed skills
#   ./install-codex.sh --uninstall  # Remove all bio-* skills

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALL_MODE="global"
PROJECT_PATH=""
VALIDATE_ONLY=false
UPDATE_MODE=false
UNINSTALL_MODE=false
VERBOSE=false

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --global              Install to ~/.codex/skills/ (default)"
    echo "  --project [PATH]      Install to .codex/skills/ in current or specified directory"
    echo "  --list                List available skills"
    echo "  --validate            Validate all skills before installing"
    echo "  --update              Only update skills that have changed"
    echo "  --uninstall           Remove all bio-* prefixed skills"
    echo "  --verbose             Show detailed output"
    echo "  --help                Show this help message"
}

list_skills() {
    echo "Available bioSkills:"
    echo ""
    find "$SCRIPT_DIR" -name "SKILL.md" -type f | sort | while read -r skill_file; do
        skill_dir=$(dirname "$skill_file")
        skill_name=$(basename "$skill_dir")
        category=$(basename "$(dirname "$skill_dir")")

        if [ "$category" = "bioSkills" ] || [ "$category" = "." ]; then
            continue
        fi

        description=$(grep "^description:" "$skill_file" | head -1 | sed 's/description: //')
        echo "  $category/$skill_name"
        if [ -n "$description" ]; then
            echo "    ${description:0:80}"
        fi
        echo ""
    done

    total=$(find "$SCRIPT_DIR" -name "SKILL.md" -type f | wc -l | tr -d ' ')
    echo "Total skills: $total"
}

validate_all_skills() {
    echo "Validating all skills..."
    echo ""

    local total=0
    local passed=0
    local failed=0
    local failed_skills=()

    while IFS= read -r skill_file; do
        skill_dir=$(dirname "$skill_file")
        skill_name=$(basename "$skill_dir")
        category=$(basename "$(dirname "$skill_dir")")

        if [ "$category" = "bioSkills" ] || [ "$category" = "." ]; then
            continue
        fi

        total=$((total + 1))
        errors=()

        if ! grep -q "^name:" "$skill_file"; then
            errors+=("Missing 'name'")
        fi

        if ! grep -q "^description:" "$skill_file"; then
            errors+=("Missing 'description'")
        else
            desc=$(grep "^description:" "$skill_file" | head -1)
            if ! echo "$desc" | grep -qi "use when"; then
                errors+=("No 'Use when' clause")
            fi
        fi

        if ! grep -q "^tool_type:" "$skill_file"; then
            errors+=("Missing 'tool_type'")
        fi

        if ! grep -q "^primary_tool:" "$skill_file"; then
            errors+=("Missing 'primary_tool'")
        fi

        if [ ! -f "$skill_dir/usage-guide.md" ]; then
            errors+=("No usage-guide.md")
        fi

        if [ ! -d "$skill_dir/examples" ] || [ -z "$(ls -A "$skill_dir/examples" 2>/dev/null)" ]; then
            errors+=("No examples/")
        fi

        if [ ${#errors[@]} -gt 0 ]; then
            failed=$((failed + 1))
            failed_skills+=("$category/$skill_name: ${errors[*]}")
            if [ "$VERBOSE" = true ]; then
                echo -e "  ${RED}FAIL${NC} $category/$skill_name"
            fi
        else
            passed=$((passed + 1))
            if [ "$VERBOSE" = true ]; then
                echo -e "  ${GREEN}PASS${NC} $category/$skill_name"
            fi
        fi
    done < <(find "$SCRIPT_DIR" -name "SKILL.md" -type f | sort)

    echo ""
    echo "Validation complete: $passed/$total passed"

    if [ $failed -gt 0 ]; then
        echo ""
        echo -e "${YELLOW}Failed skills:${NC}"
        for skill in "${failed_skills[@]}"; do
            echo "  - $skill"
        done
        return 1
    fi
    return 0
}

install_skills() {
    local target_dir="$1"

    if [ ! -d "$(dirname "$target_dir")" ]; then
        echo -e "${RED}Error: Parent directory does not exist: $(dirname "$target_dir")${NC}"
        exit 1
    fi

    echo "Installing bioSkills to: $target_dir"
    echo ""

    mkdir -p "$target_dir"

    local installed=0
    local skipped=0
    local errors=0

    while IFS= read -r skill_file; do
        skill_dir=$(dirname "$skill_file")
        skill_name=$(basename "$skill_dir")
        category=$(basename "$(dirname "$skill_dir")")

        if [ "$category" = "bioSkills" ] || [ "$category" = "." ]; then
            continue
        fi

        full_skill_name="bio-${category}-${skill_name}"
        target_skill_dir="$target_dir/$full_skill_name"

        if [ "$UPDATE_MODE" = true ] && [ -d "$target_skill_dir" ]; then
            src_time=$(stat -f %m "$skill_file" 2>/dev/null || stat -c %Y "$skill_file" 2>/dev/null)
            dst_time=$(stat -f %m "$target_skill_dir/SKILL.md" 2>/dev/null || stat -c %Y "$target_skill_dir/SKILL.md" 2>/dev/null || echo 0)

            if [ "$src_time" -le "$dst_time" ]; then
                skipped=$((skipped + 1))
                if [ "$VERBOSE" = true ]; then
                    echo "  Skipped (unchanged): $full_skill_name"
                fi
                continue
            fi
        fi

        if ! mkdir -p "$target_skill_dir" 2>/dev/null; then
            echo -e "  ${RED}Error creating: $full_skill_name${NC}"
            errors=$((errors + 1))
            continue
        fi

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

        installed=$((installed + 1))
        echo "  Installed: $full_skill_name"
    done < <(find "$SCRIPT_DIR" -name "SKILL.md" -type f | sort)

    echo ""
    echo "Installation complete."
    echo "  Installed: $installed"
    if [ "$UPDATE_MODE" = true ]; then
        echo "  Skipped (unchanged): $skipped"
    fi
    if [ $errors -gt 0 ]; then
        echo -e "  ${RED}Errors: $errors${NC}"
    fi
}

uninstall_skills() {
    local target_dir="$1"

    if [ ! -d "$target_dir" ]; then
        echo "No skills directory found at: $target_dir"
        exit 0
    fi

    echo "Removing bioSkills from: $target_dir"
    echo ""

    local removed=0

    for skill_dir in "$target_dir"/bio-*; do
        if [ -d "$skill_dir" ]; then
            skill_name=$(basename "$skill_dir")
            rm -rf "$skill_dir"
            echo "  Removed: $skill_name"
            removed=$((removed + 1))
        fi
    done

    echo ""
    echo "Uninstall complete. Removed $removed skills."
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
        --validate)
            VALIDATE_ONLY=true
            shift
            ;;
        --update)
            UPDATE_MODE=true
            shift
            ;;
        --uninstall)
            UNINSTALL_MODE=true
            shift
            ;;
        --verbose|-v)
            VERBOSE=true
            shift
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

if [ "$VALIDATE_ONLY" = true ]; then
    if validate_all_skills; then
        echo ""
        echo -e "${GREEN}All skills passed validation.${NC}"
        exit 0
    else
        echo ""
        echo -e "${RED}Some skills failed validation.${NC}"
        exit 1
    fi
fi

if [ "$UNINSTALL_MODE" = true ]; then
    uninstall_skills "$TARGET_DIR"
    exit 0
fi

install_skills "$TARGET_DIR"

echo ""
echo "Skills are now available in Codex CLI."
echo "They will be auto-invoked when relevant to your bioinformatics tasks."
