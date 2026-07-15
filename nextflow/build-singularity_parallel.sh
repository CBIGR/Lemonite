#!/bin/bash

# Build script for Lemonite Singularity container

set -e

echo "Building Lemonite Analysis Pipeline Singularity image..."
echo "=================================================="

# Check if singularity is installed
if ! command -v singularity &> /dev/null; then
    echo "Error: Singularity is not installed or not in PATH"
    echo "Please install Singularity first: https://sylabs.io/guides/3.0/user-guide/installation.html"
    exit 1
fi

# HPC-UGent's scratch filesystems are mounted with 'nodev', which breaks the proot/ptrace
# fakeroot emulation that Apptainer falls back to for the final squashfs creation step.
# Per https://docs.hpc.ugent.be/Linux/apptainer/ : build using the workernode's local /tmp,
# and move the finished .sif to scratch afterwards.
BUILD_CACHEDIR="/tmp/${USER}/apptainer/cache"
BUILD_TMPDIR="/tmp/${USER}/apptainer/tmpdir"
BUILD_OUTDIR="/tmp/${USER}/apptainer/build"
export SINGULARITY_CACHEDIR="$BUILD_CACHEDIR"
export SINGULARITY_TMPDIR="$BUILD_TMPDIR"
export APPTAINER_CACHEDIR="$BUILD_CACHEDIR"
export APPTAINER_TMPDIR="$BUILD_TMPDIR"
# NOTE: deliberately do NOT export TMPDIR here. Exporting TMPDIR causes Apptainer
# to bind the host build-tmpdir as the container's /tmp, which the --fakeroot user
# namespace cannot write to. apt-key then fails with
#   "Couldn't create temporary file /tmp/apt.conf.XXXXXX"
# and every repo is reported as "not signed", aborting %post with exit 100.
# APPTAINER_TMPDIR (for Apptainer's own staging on local disk) is enough.
mkdir -p "$BUILD_CACHEDIR" "$BUILD_TMPDIR" "$BUILD_OUTDIR"

# This is an unprivileged Apptainer install (no setuid 'starter-suid') and this user has
# no /etc/subuid mapping, so --fakeroot falls back to proot for the final squashfs
# creation. On these RHEL9 workernodes proot's ptrace is blocked by seccomp:
#   proot error: ptrace(TRACEME): Operation not permitted
# Every route through `apptainer build`'s SIF/squashfs step hits this (tested:
# PROOT_NO_SECCOMP=1, --userns, sandbox->sif conversion, system mksquashfs on PATH).
# WORKAROUND: build a --sandbox (a plain rootfs directory, no squashfs => no proot),
# which succeeds, and run the pipeline directly from the sandbox directory. Nextflow/
# Apptainer can `exec` a sandbox dir exactly like a .sif.

echo "Using cache directory: $SINGULARITY_CACHEDIR"
echo "Using temp directory: $SINGULARITY_TMPDIR"

# Check available space on local disk
USER_SPACE=$(df -h "$BUILD_TMPDIR" | awk 'NR==2 {print $4}')
echo "Available space on local disk: $USER_SPACE"

# Build a SANDBOX (rootfs directory) on local disk, then move it to scratch.
# See the note above for why we cannot produce a .sif on these nodes.
SANDBOX_NAME="lemontree-pipeline_v1.0.0_parallel.sandbox"
echo "Building Singularity sandbox: ${SANDBOX_NAME}"
echo "Note: Building without sudo using --fakeroot --sandbox (avoids proot/squashfs)"
# --fix-perms ensures the owner gets rwX on all dirs/files, so the sandbox can be
# relocated and later removed without "Permission denied" on fakeroot-created dirs.
singularity build --force --fakeroot --fix-perms --sandbox "$BUILD_OUTDIR/${SANDBOX_NAME}" Singularity_parallel.def

# Relocate the sandbox from local disk to scratch.
# A --fakeroot sandbox contains dirs created under the fake-root namespace (e.g.
# /var/cache/apt/archives/partial, mode 0700 owned by a fake uid) that our REAL uid
# cannot traverse, and files carry a 'user.rootlesscontainers' xattr that plain
# mv/cp trip over with "Operation not permitted". Both broke a naive `mv`.
#
# Fix:
#   1. Drop the apt partial/cache dirs we don't need (they cause most of the perm errors).
#   2. chmod the tree so our uid can read/traverse everything.
#   3. Use tar with --xattrs to stream the tree to scratch (preserves the
#      rootlesscontainers ownership-faking xattr that the sandbox needs at exec time).
echo "Relocating sandbox to scratch..."
chmod -R u+rwX "$BUILD_OUTDIR/${SANDBOX_NAME}" 2>/dev/null || true
rm -rf "$BUILD_OUTDIR/${SANDBOX_NAME}/var/cache/apt/archives/partial" \
       "$BUILD_OUTDIR/${SANDBOX_NAME}/var/lib/apt/lists/partial" 2>/dev/null || true

rm -rf "${PWD}/${SANDBOX_NAME}"
mkdir -p "${PWD}/${SANDBOX_NAME}"
# tar preserves perms + xattrs across the filesystem boundary where mv/cp fail.
tar --xattrs --xattrs-include='*' -C "$BUILD_OUTDIR/${SANDBOX_NAME}" -cf - . \
  | tar --xattrs --xattrs-include='*' -C "${PWD}/${SANDBOX_NAME}" -xf -
echo "Sandbox relocated to ${PWD}/${SANDBOX_NAME}"

# Cleanup function
cleanup_cache() {
    echo "Cleaning up build cache..."
    if [ -d "$BUILD_CACHEDIR" ]; then
        du -sh "$BUILD_CACHEDIR" 2>/dev/null || echo "Cache directory size: unknown"
        rm -rf "$BUILD_CACHEDIR"
        echo "Cache directory removed"
    fi
    rm -rf "$BUILD_TMPDIR" "$BUILD_OUTDIR"
}

# Verify the build
if [ -d "${SANDBOX_NAME}" ]; then
    echo ""
    echo "✅ Parallel Singularity sandbox built successfully!"
    echo "Sandbox size: $(du -sh "${SANDBOX_NAME}" 2>/dev/null | cut -f1)"
    echo ""
    echo "Smoke-testing the sandbox..."
    singularity exec "${PWD}/${SANDBOX_NAME}" R --version | head -1
    singularity exec "${PWD}/${SANDBOX_NAME}" python3 --version
    echo ""
    cleanup_cache
    echo ""
    echo "Run the pipeline with:"
    echo "nextflow run main.nf -profile singularity --input_dir /path/to/data"
    echo "(conf/singularity.config points 'container' at ${SANDBOX_NAME})"
else
    echo "❌ Failed to build Singularity sandbox"
    exit 1
fi
