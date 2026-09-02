#!/usr/bin/env bash

set -euo pipefail

script_dir="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(CDPATH= cd -- "${script_dir}/.." && pwd)"

cspice_dir="${CSPICE_DIR:-${project_root}/external/cspice}"
kernel_dir="${SPICE_KERNEL_DIR:-${project_root}/assets/spice}"
naif_root="https://naif.jpl.nasa.gov/pub/naif"
cspice_version="N0067"

spk_name="de442s.bsp"
lsk_name="naif0012.tls"
pck_name="pck00011.tpc"

# NAIF publishes MD5 checksums for planetary SPKs. The text-kernel hashes pin
# the exact files used by this project so a changed download fails explicitly.
spk_md5="cc49327e06088124c0e39d8dde9f0b58"
lsk_md5="25a2fff30b0dedb4d76c06727b1895b1"
pck_md5="3c0bdc01bb4101d6242cf96ece67a586"

download_toolkit=true
download_kernels=true
force=false

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Download the CSPICE N0067 toolkit and the project's pinned generic kernels.

Options:
  --toolkit-only  Download only CSPICE into external/cspice
  --kernels-only  Download only kernels into assets/spice
  --force         Replace existing toolkit and kernel files
  -h, --help      Show this help

Environment overrides:
  CSPICE_DIR              Toolkit installation directory
  SPICE_KERNEL_DIR        Kernel installation directory
  CSPICE_PLATFORM         NAIF platform directory override
  CSPICE_PACKAGE_URL      CSPICE archive URL for an unsupported platform
  CSPICE_ARCHIVE_FORMAT   tar.Z, tar.gz, or zip for an override
EOF
}

fail() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}

for arg in "$@"; do
    case "${arg}" in
        --toolkit-only)
            download_kernels=false
            ;;
        --kernels-only)
            download_toolkit=false
            ;;
        --force)
            force=true
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            fail "unknown option: ${arg}"
            ;;
    esac
done

command -v curl >/dev/null 2>&1 || fail "curl is required"

tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/astrolib-spice.XXXXXX")"
trap 'rm -rf "${tmp_dir}"' EXIT

download_file() {
    local url="$1"
    local output="$2"

    printf 'Downloading %s\n' "${url}"
    curl --fail --location --retry 3 --retry-delay 2 \
        --connect-timeout 20 --output "${output}" "${url}"
    [[ -s "${output}" ]] || fail "downloaded file is empty: ${url}"
}

file_md5() {
    local filepath="$1"

    if command -v md5 >/dev/null 2>&1; then
        md5 -q "${filepath}"
    elif command -v md5sum >/dev/null 2>&1; then
        md5sum "${filepath}" | awk '{print $1}'
    else
        fail "md5 or md5sum is required to verify kernels"
    fi
}

verify_md5() {
    local filepath="$1"
    local expected="$2"
    local actual

    actual="$(file_md5 "${filepath}")"
    [[ "${actual}" == "${expected}" ]] || \
        fail "checksum mismatch for ${filepath}: expected ${expected}, got ${actual}"
}

select_cspice_package() {
    local os
    local arch

    if [[ -n "${CSPICE_PACKAGE_URL:-}" ]]; then
        cspice_platform="custom"
        cspice_url="${CSPICE_PACKAGE_URL}"
        cspice_format="${CSPICE_ARCHIVE_FORMAT:-tar.Z}"
        return
    fi

    os="$(uname -s)"
    arch="$(uname -m)"

    if [[ -n "${CSPICE_PLATFORM:-}" ]]; then
        cspice_platform="${CSPICE_PLATFORM}"
        case "${cspice_platform}" in
            PC_Windows_VisualC_*) cspice_format="zip" ;;
            PC_Cygwin_GCC_*) cspice_format="tar.gz" ;;
            *) cspice_format="tar.Z" ;;
        esac
    else
        case "${os}:${arch}" in
            Darwin:arm64|Darwin:aarch64)
                cspice_platform="MacM1_OSX_clang_64bit"
                cspice_format="tar.Z"
                ;;
            Darwin:x86_64|Darwin:amd64)
                cspice_platform="MacIntel_OSX_AppleC_64bit"
                cspice_format="tar.Z"
                ;;
            Linux:x86_64|Linux:amd64)
                cspice_platform="PC_Linux_GCC_64bit"
                cspice_format="tar.Z"
                ;;
            CYGWIN_NT-*:x86_64|CYGWIN_NT-*:amd64)
                cspice_platform="PC_Cygwin_GCC_64bit"
                cspice_format="tar.gz"
                ;;
            MINGW64_NT-*:x86_64|MINGW64_NT-*:amd64|MSYS_NT-*:x86_64|MSYS_NT-*:amd64)
                cspice_platform="PC_Windows_VisualC_64bit"
                cspice_format="zip"
                ;;
            *)
                fail "unsupported CSPICE platform ${os}/${arch}; set CSPICE_PLATFORM or CSPICE_PACKAGE_URL and CSPICE_ARCHIVE_FORMAT"
                ;;
        esac
    fi

    cspice_url="${naif_root}/toolkit/C/${cspice_platform}/packages/cspice.${cspice_format}"
}

validate_cspice_tree() {
    local directory="$1"

    [[ -f "${directory}/${cspice_version}" ]] || return 1
    [[ -f "${directory}/include/SpiceUsr.h" ]] || return 1
    [[ -f "${directory}/lib/cspice.a" || -f "${directory}/lib/cspice.lib" ]] || return 1
}

install_cspice() {
    local archive
    local extract_dir="${tmp_dir}/cspice-extract"
    local extracted_tree

    if [[ -d "${cspice_dir}" && "${force}" == false ]]; then
        validate_cspice_tree "${cspice_dir}" || \
            fail "existing CSPICE directory is incomplete; rerun with --force"
        printf 'Using existing CSPICE %s: %s\n' "${cspice_version}" "${cspice_dir}"
        return
    fi

    select_cspice_package
    archive="${tmp_dir}/cspice.${cspice_format}"
    mkdir -p "${extract_dir}"
    download_file "${cspice_url}" "${archive}"

    case "${cspice_format}" in
        tar.Z)
            if command -v gzip >/dev/null 2>&1; then
                gzip -dc "${archive}" | tar -xf - -C "${extract_dir}"
            elif command -v uncompress >/dev/null 2>&1; then
                uncompress -c "${archive}" | tar -xf - -C "${extract_dir}"
            else
                fail "gzip or uncompress is required to extract CSPICE"
            fi
            ;;
        tar.gz)
            tar -xzf "${archive}" -C "${extract_dir}"
            ;;
        zip)
            command -v unzip >/dev/null 2>&1 || fail "unzip is required to extract CSPICE"
            unzip -q "${archive}" -d "${extract_dir}"
            ;;
        *)
            fail "unsupported CSPICE archive format: ${cspice_format}"
            ;;
    esac

    extracted_tree="${extract_dir}/cspice"
    validate_cspice_tree "${extracted_tree}" || fail "downloaded CSPICE tree failed validation"

    mkdir -p "$(dirname -- "${cspice_dir}")"
    if [[ -e "${cspice_dir}" ]]; then
        rm -rf "${cspice_dir}"
    fi
    mv "${extracted_tree}" "${cspice_dir}"
    printf 'Installed CSPICE %s for %s: %s\n' \
        "${cspice_version}" "${cspice_platform}" "${cspice_dir}"
}

install_kernel() {
    local name="$1"
    local url="$2"
    local expected_md5="$3"
    local destination="${kernel_dir}/${name}"
    local temporary="${tmp_dir}/${name}"

    if [[ -f "${destination}" && "${force}" == false ]]; then
        verify_md5 "${destination}" "${expected_md5}"
        printf 'Using existing kernel: %s\n' "${destination}"
        return
    fi

    download_file "${url}" "${temporary}"
    verify_md5 "${temporary}" "${expected_md5}"
    mkdir -p "${kernel_dir}"
    mv "${temporary}" "${destination}"
    printf 'Installed kernel: %s\n' "${destination}"
}

if [[ "${download_toolkit}" == true ]]; then
    install_cspice
fi

if [[ "${download_kernels}" == true ]]; then
    install_kernel \
        "${spk_name}" \
        "${naif_root}/generic_kernels/spk/planets/${spk_name}" \
        "${spk_md5}"
    install_kernel \
        "${lsk_name}" \
        "${naif_root}/generic_kernels/lsk/${lsk_name}" \
        "${lsk_md5}"
    install_kernel \
        "${pck_name}" \
        "${naif_root}/generic_kernels/pck/${pck_name}" \
        "${pck_md5}"
fi

printf 'SPICE setup complete.\n'
