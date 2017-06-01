#!/bin/bash

VERSION_FILE="Version"
NAME="nerstrand"

die() {
  echo "ERROR: ${@}" 2>&1
  exit 1
}

VERSION="$(awk '/^#define NERSTRAND_VER_MAJOR/{M=$3}
                /^#define NERSTRAND_VER_MINOR/{m=$3}
                /^#define NERSTRAND_VER_SUBMINOR/{s=$3}
                END {print M "." m "." s}' "include/${NAME}.h")"

echo "Packaging version $VERSION..."

PKGNAME="${NAME}-${VERSION}"
PKGDIR="${PKGNAME}"

mkdir -vp "${PKGDIR}"

echo "${VERSION}" > "${PKGDIR}/${VERSION_FILE}"

tar cfh - --exclude-vcs --exclude-backups --exclude='build' \
      --exclude='*.o' --exclude="${PKGNAME}.tar.gz" --exclude="${0}" \
      --exclude="${PKGDIR}" --exclude="bowstring/domlib" \
      --exclude="bowstring/build" --exclude="bowstring/src/mpi" \
      --exclude="domlib/libdom.a" --exclude="domlib/test" \
      --exclude="domlib/Makefile" --exclude="Makefile" \
      --exclude="api" --exclude="doxygen.conf" \
      ./ | tar xf - -C "${PKGDIR}" \
      || die "Failed to copy package contents to '${PKGDIR}'"

tar czvf "${PKGNAME}.tar.gz" "${PKGDIR}" || die "Failed to create '${PKGNAME}"

rm -rfv "${PKGDIR}"

echo "Finished packaging"

