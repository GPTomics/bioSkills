#!/bin/bash
#
# Shared functions for bioSkills install scripts
# Sourced by install-claude.sh, install-codex.sh, install-gemini.sh, install-openclaw.sh
#
# Each installer must define before calling run_installer:
#   TOOL_NAME           - Display name (e.g., "Claude Code")
#   DEFAULT_TARGET_DIR  - Global install path (e.g., "$HOME/.claude/skills")
#   PROJECT_SUBDIR      - Project-level subdirectory (e.g., ".claude/skills")
#   copy_skill_files()  - Platform-specific file copy: copy_skill_files src_dir target_dir
#   print_usage()       - Help text (call print_common_options for shared flags)
#
# Optional overrides:
#   parse_extra_arg()      - Handle platform-specific flags; set EXTRA_SHIFT; return 0 if handled
#   post_install_message() - Custom post-install output
#   check_project_conflicts() - Pre-install conflict detection for project installs

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

INSTALL_MODE="global"
PROJECT_PATH=""
VALIDATE_ONLY=false
UPDATE_MODE=false
UNINSTALL_MODE=false
VERBOSE=false
DRY_RUN=false
CATEGORY_FILTER=""
FORCE=false

print_common_options() {
    echo "  --global              Install to default global location (default)"
    echo "  --project [PATH]      Install to project/workspace directory"
    echo "  --categories CATS     Install only specified categories (comma-separated)"
    echo "  --list                List available skills"
    echo "  --validate            Validate all skills before installing"
    echo "  --update              Only update skills that have changed"
    echo "  --uninstall           Remove all bio-* prefixed skills"
    echo "  --dry-run             Preview what would be installed"
    echo "  --verbose             Show detailed output"
    echo "  --force               Force install even if target has existing content"
    echo "  --help                Show this help message"
}

list_skills() {
    echo "Available bioSkills:"
    echo ""
    find "$SCRIPT_DIR" -name "SKILL.md" -type f | sort | while read -r skill_file; do
        local skill_dir=$(dirname "$skill_file")
        local skill_name=$(basename "$skill_dir")
        local category=$(basename "$(dirname "$skill_dir")")

        [ "$category" = "$(basename "$SCRIPT_DIR")" ] || [ "$category" = "." ] && continue

        local description=$(grep "^description:" "$skill_file" | head -1 | sed 's/description: //')
        echo "  $category/$skill_name"
        [ -n "$description" ] && echo "    ${description:0:80}"
        echo ""
    done

    local total=$(find "$SCRIPT_DIR" -mindepth 3 -name "SKILL.md" -type f | wc -l | tr -d ' ')
    echo "Total skills: $total"
}

validate_all_skills() {
    echo "Validating all skills..."
    echo ""

    local total=0 passed=0 failed=0
    local failed_skills=()

    while IFS= read -r skill_file; do
        local skill_dir=$(dirname "$skill_file")
        local skill_name=$(basename "$skill_dir")
        local category=$(basename "$(dirname "$skill_dir")")

        [ "$category" = "$(basename "$SCRIPT_DIR")" ] || [ "$category" = "." ] && continue

        total=$((total + 1))
        local errors=()

        grep -q "^name:" "$skill_file" || errors+=("Missing 'name'")
        if grep -q "^description:" "$skill_file"; then
            local desc=$(grep "^description:" "$skill_file" | head -1)
            echo "$desc" | grep -qi "use when" || errors+=("No 'Use when' clause")
        else
            errors+=("Missing 'description'")
        fi
        grep -q "^tool_type:" "$skill_file" || errors+=("Missing 'tool_type'")
        grep -q "^primary_tool:" "$skill_file" || errors+=("Missing 'primary_tool'")
        [ -f "$skill_dir/usage-guide.md" ] || errors+=("No usage-guide.md")
        if [ ! -d "$skill_dir/examples" ] || [ -z "$(ls -A "$skill_dir/examples" 2>/dev/null)" ]; then
            errors+=("No examples/")
        fi

        if [ ${#errors[@]} -gt 0 ]; then
            failed=$((failed + 1))
            failed_skills+=("$category/$skill_name: ${errors[*]}")
            if [ "$VERBOSE" = true ]; then
                echo -e "  ${RED}FAIL${NC} $category/$skill_name"
                for err in "${errors[@]}"; do
                    echo "       - $err"
                done
            fi
        else
            passed=$((passed + 1))
            [ "$VERBOSE" = true ] && echo -e "  ${GREEN}PASS${NC} $category/$skill_name"
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

validate_categories() {
    [ -z "$CATEGORY_FILTER" ] && return 0

    local available=()
    for d in "$SCRIPT_DIR"/*/; do
        [ ! -d "$d" ] && continue
        if find "$d" -maxdepth 2 -name "SKILL.md" -print -quit 2>/dev/null | grep -q .; then
            available+=("$(basename "$d")")
        fi
    done

    for cat in ${CATEGORY_FILTER//,/ }; do
        local found=false
        for a in "${available[@]}"; do
            [ "$cat" = "$a" ] && found=true && break
        done
        if [ "$found" = false ]; then
            echo -e "${RED}Error: Unknown category '$cat'${NC}"
            echo "Available categories:"
            printf '  %s\n' "${available[@]}" | sort
            exit 1
        fi
    done
}

category_matches() {
    local category="$1"
    [ -z "$CATEGORY_FILTER" ] && return 0
    echo ",$CATEGORY_FILTER," | grep -q ",$category," && return 0
    return 1
}

install_skills() {
    local target_dir="$1"

    if [ "$DRY_RUN" = false ]; then
        if [ "$INSTALL_MODE" = "project" ]; then
            local base_dir="${PROJECT_PATH:-$(pwd)}"
            if [ ! -d "$base_dir" ]; then
                echo -e "${RED}Error: Project directory does not exist: $base_dir${NC}"
                exit 1
            fi
        elif [ ! -d "$(dirname "$target_dir")" ]; then
            echo -e "${RED}Error: Parent directory does not exist: $(dirname "$target_dir")${NC}"
            exit 1
        fi
        if [ -d "$target_dir" ] && [ ! -w "$target_dir" ]; then
            echo -e "${RED}Error: Cannot write to target directory: $target_dir${NC}"
            exit 1
        fi
        if type check_project_conflicts &>/dev/null; then
            check_project_conflicts "$target_dir"
        fi
    fi

    if [ "$DRY_RUN" = true ]; then
        echo "Dry run - would install to: $target_dir"
    else
        echo "Installing bioSkills to: $target_dir"
    fi
    echo ""

    [ "$DRY_RUN" = false ] && mkdir -p "$target_dir"

    local installed=0 skipped=0 errors=0
    local total_bytes=0 desc_bytes=0

    while IFS= read -r skill_file; do
        local skill_dir=$(dirname "$skill_file")
        local skill_name=$(basename "$skill_dir")
        local category=$(basename "$(dirname "$skill_dir")")

        [ "$category" = "$(basename "$SCRIPT_DIR")" ] || [ "$category" = "." ] && continue
        category_matches "$category" || continue

        local full_skill_name="bio-${category}-${skill_name}"
        local target_skill_dir="$target_dir/$full_skill_name"

        if [ "$UPDATE_MODE" = true ] && [ -d "$target_skill_dir" ]; then
            local src_time=$(stat -f %m "$skill_file" 2>/dev/null || stat -c %Y "$skill_file" 2>/dev/null)
            local dst_time=$(stat -f %m "$target_skill_dir/SKILL.md" 2>/dev/null || stat -c %Y "$target_skill_dir/SKILL.md" 2>/dev/null || echo 0)
            if [ "$src_time" -le "$dst_time" ]; then
                skipped=$((skipped + 1))
                [ "$VERBOSE" = true ] && echo "  Skipped (unchanged): $full_skill_name"
                continue
            fi
        fi

        if [ "$DRY_RUN" = true ]; then
            installed=$((installed + 1))
            total_bytes=$((total_bytes + $(wc -c < "$skill_file")))
            local db=$(grep "^description:" "$skill_file" | head -1 | wc -c)
            desc_bytes=$((desc_bytes + db))
            [ "$VERBOSE" = true ] && echo "  Would install: $full_skill_name"
            continue
        fi

        if ! mkdir -p "$target_skill_dir" 2>/dev/null; then
            echo -e "  ${RED}Error creating: $full_skill_name${NC}"
            errors=$((errors + 1))
            continue
        fi

        if ! copy_skill_files "$skill_dir" "$target_skill_dir"; then
            echo -e "  ${RED}Error copying: $full_skill_name${NC}"
            errors=$((errors + 1))
            continue
        fi

        installed=$((installed + 1))
        echo "  Installed: $full_skill_name"
    done < <(find "$SCRIPT_DIR" -name "SKILL.md" -type f | sort)

    echo ""
    if [ "$DRY_RUN" = true ]; then
        echo "Would install: $installed skills"
        local desc_tokens=$(( (195 + installed * 97 + desc_bytes) / 4 ))
        local avg_body=$(( total_bytes / (installed > 0 ? installed : 1) / 4 ))
        echo "Estimated always-loaded overhead (descriptions): ~${desc_tokens} tokens"
        echo "Average per-skill body (loaded on activation): ~${avg_body} tokens"
    else
        echo "Installation complete."
        echo "  Installed: $installed"
        [ "$UPDATE_MODE" = true ] && echo "  Skipped (unchanged): $skipped"
        [ $errors -gt 0 ] && echo -e "  ${RED}Errors: $errors${NC}"
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
            local skill_name=$(basename "$skill_dir")
            rm -rf "$skill_dir"
            echo "  Removed: $skill_name"
            removed=$((removed + 1))
        fi
    done

    echo ""
    echo "Uninstall complete. Removed $removed skills."

    if [ -d "$target_dir" ] && [ -z "$(ls -A "$target_dir")" ]; then
        rmdir "$target_dir" 2>/dev/null && echo "Removed empty skills directory."
    fi
}

run_installer() {
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
            --categories)
                if [[ -n "$2" && ! "$2" =~ ^-- ]]; then
                    CATEGORY_FILTER=$(echo "$2" | tr -d ' ')
                    shift
                else
                    echo "Error: --categories requires a comma-separated list of categories"
                    exit 1
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
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            --force)
                FORCE=true
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
                EXTRA_SHIFT=0
                if type parse_extra_arg &>/dev/null && parse_extra_arg "$@"; then
                    shift $EXTRA_SHIFT
                else
                    echo "Unknown option: $1"
                    print_usage
                    exit 1
                fi
                ;;
        esac
    done

    if [ "$INSTALL_MODE" = "global" ]; then
        TARGET_DIR="$DEFAULT_TARGET_DIR"
    else
        if [ -n "$PROJECT_PATH" ]; then
            TARGET_DIR="$PROJECT_PATH/$PROJECT_SUBDIR"
        else
            TARGET_DIR="$(pwd)/$PROJECT_SUBDIR"
        fi
    fi

    [ -n "$CATEGORY_FILTER" ] && validate_categories

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
    if type post_install_message &>/dev/null; then
        post_install_message
    else
        echo "Skills are now available in $TOOL_NAME."
        echo "They will be auto-invoked when relevant to your bioinformatics tasks."
    fi
}
