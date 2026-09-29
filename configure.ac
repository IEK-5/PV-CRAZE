# pv-craze version
m4_define([pvcraze_version], m4_esyscmd_s([cat VERSION]))
AC_INIT([libpvcraze], [pvcraze_version])

# interface version 
m4_define([pvcraze_major],   [1])
m4_define([pvcraze_revision],[0])
m4_define([pvcraze_patch],   [1])
version_info="pvcraze_major:pvcraze_revision:pvcraze_patch"
AC_SUBST([version_info])


AC_CONFIG_MACRO_DIR([m4])
LT_INIT([])
AC_PROG_SED()
AC_PROG_AWK()
AC_PROG_EGREP()
AC_ENABLE_SHARED
AC_ENABLE_SHARED([])
LT_INIT([win32-dll])
AM_INIT_AUTOMAKE()

# optionally enable linting with clang-tidy
AC_ARG_ENABLE([lint], [AS_HELP_STRING([--enable-lint],[Enable linting])],
  [enable_lint=yes],
  [enable_lint=no])
if test "$enable_lint" = "yes"; then
  AC_PATH_PROG([LINTER], [clang-tidy], [true])
  AS_IF([test "x$LINTER" = xtrue],
    [AC_MSG_WARN([cland-tidy not found! Linting disabled.])])
fi
AM_CONDITIONAL([RUN_LINTER], [test "x$enable_lint" = xyes])

AC_PROG_CC
AC_CONFIG_SRCDIR([src/lib-pv-craze.c])
AC_CONFIG_HEADERS([src/config.h])
AC_CHECK_LIB(m,main)
AX_BLAS()
AC_CHECK_LIB(gslcblas,main)
AC_CHECK_LIB(gsl,main, [], [AC_MSG_ERROR([unable to find the GNU Schientific Library])])
AX_LAPACK()
LIBS="$BLAS_LIBS $FLIBS $LIBS"
     
AC_SUBST(LIBVERSION, [pvcraze_version])
     
AC_CONFIG_FILES([
 Makefile src/Makefile
])
AC_OUTPUT
