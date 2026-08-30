#!/usr/bin/env bash

# Create a tagged GBS-Typer release.
#
# Usage:
#   ./release.sh 1.0.0
#
# This will:
#   1. Confirm the version does not already exist.
#   2. Confirm the current branch is main.
#   3. Confirm the working tree is clean.
#   4. Update the GBS-Typer container tag in nextflow.config.
#   5. Commit the container-version change.
#   6. Create an annotated Git tag, e.g. v1.0.0.
#   7. Push main and the release tag.

set -euo pipefail


validate_version_format() {
    local version="$1"

    if [[ ! "${version}" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
        echo "Invalid version: ${version}"
        echo "Expected semantic version format, for example: 1.0.0"
        exit 1
    fi
}


validate_new_version() {
    local version="$1"
    local tag="v${version}"

    git fetch --tags

    if git tag --list "${tag}" | grep -Fxq "${tag}"; then
        echo "Version ${tag} already exists."
        exit 1
    fi
}


validate_branch() {
    local branch

    branch="$(git rev-parse --abbrev-ref HEAD)"

    if [[ "${branch}" != "main" ]]; then
        echo "Releases can only be made from the main branch."
        echo "Current branch: ${branch}"
        exit 1
    fi
}


validate_working_tree_clean() {
    if [[ -n "$(git status --porcelain)" ]]; then
        echo "Uncommitted changes were found."
        echo "Commit or stash them before creating a release."
        git status --short
        exit 1
    fi
}


update_container_version() {
    local version="$1"
    local config_file="nextflow.config"
    local temporary_file

    temporary_file="$(mktemp)"

    sed -E \
        "s@(gbs-typer-sanger-nf:)[^'\"[:space:]]+@\1${version}@g" \
        "${config_file}" \
        > "${temporary_file}"

    if cmp -s "${config_file}" "${temporary_file}"; then
        rm -f "${temporary_file}"
        echo "The container version was not changed."
        echo "Check that nextflow.config contains:"
        echo "gbs-typer-sanger-nf:<version>"
        exit 1
    fi

    mv "${temporary_file}" "${config_file}"

    echo "Updated container version:"
    grep -n "gbs-typer-sanger-nf:" "${config_file}"
}


create_release() {
    local version="$1"
    local tag="v${version}"
    local config_file="nextflow.config"

    update_container_version "${version}"

    git add "${config_file}"

    git commit \
        -m "Update container version for release ${tag}" \
        "${config_file}"

    echo "Creating release tag ${tag}..."

    git tag -a \
        "${tag}" \
        -m "GBS-Typer release ${tag}"

    echo "Pushing main and ${tag}..."

    git push origin main
    git push origin "${tag}"

    echo "Release ${tag} completed."
}


main() {
    if [[ $# -ne 1 ]]; then
        echo "Usage: $0 <version without v>"
        echo "Example: $0 1.0.0"
        exit 1
    fi

    local version="$1"

    validate_version_format "${version}"
    validate_new_version "${version}"
    validate_branch
    validate_working_tree_clean
    create_release "${version}"
}


main "$@"