process get_version {
    label 'farm_mid'
    
    output:
    path "version.txt"

    script:
    """
    set -euo pipefail

    echo "version" > version.txt

    if [ -z "${params.version}" ]; then
      v="\$(git -C "${workflow.projectDir}" describe --tags --always --dirty 2>/dev/null || true)"
      [ -z "\$v" ] && v="unknown"
      echo "\$v" >> version.txt
    else
      echo "${params.version}" >> version.txt
    fi
    """
}
