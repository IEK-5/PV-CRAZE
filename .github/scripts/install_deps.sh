#!/bin/bash
# Install GSL (>= 2.2) and LAPACK for a cibuildwheel build (before-all).
#
# Linux:  runs inside the manylinux_2_28 (AlmaLinux 8) container; GSL 2.5 and
#         OpenBLAS come from the distribution, auditwheel bundles them.
# macOS:  Homebrew's bottles target the runner's macOS version, so GSL is
#         built from source as a static, position-independent library for the
#         wheel's deployment target. LAPACK comes from Accelerate.
set -euo pipefail

GSL_VERSION=2.8
GSL_SHA256=6a99eeed15632c6354895b1dd542ed5a855c0f15d9ad1326c6fe2b2c9e423190

case "$(uname -s)" in
Linux)
	dnf install -y dnf-plugins-core
	dnf config-manager --set-enabled powertools
	dnf install -y gsl-devel openblas-devel
	;;
Darwin)
	prefix="${GSL_PREFIX:-/tmp/gsl}"
	case "$(uname -m)" in
	arm64) export MACOSX_DEPLOYMENT_TARGET=11.0 ;;
	*) export MACOSX_DEPLOYMENT_TARGET=10.13 ;;
	esac
	work=$(mktemp -d)
	curl -sSfL -o "$work/gsl.tar.gz" "https://ftp.gnu.org/gnu/gsl/gsl-$GSL_VERSION.tar.gz"
	echo "$GSL_SHA256  $work/gsl.tar.gz" | shasum -a 256 -c -
	tar -xzf "$work/gsl.tar.gz" -C "$work"
	cd "$work/gsl-$GSL_VERSION"
	./configure --prefix="$prefix" --disable-shared --enable-static --with-pic
	make -j"$(sysctl -n hw.ncpu)"
	make install
	;;
*)
	echo "install_deps.sh: unsupported platform $(uname -s)" >&2
	exit 1
	;;
esac
