m4_define([pvcraze_major_version],    [1])
m4_define([pvcraze_minor_version],    [0])
m4_define([pvcraze_revis_version],    [0])
m4_define([pvcraze_version],[pvcraze_major_version:pvcraze_minor_version:pvcraze_revis_version])
AC_INIT([libpvcraze], [pvcraze_version])

AC_CONFIG_MACRO_DIR([m4])
LT_INIT([])
AC_PROG_SED()
AC_PROG_AWK()
AC_PROG_EGREP()
AC_ENABLE_SHARED
AC_ENABLE_SHARED([])
LT_INIT([win32-dll])
AM_INIT_AUTOMAKE()
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
