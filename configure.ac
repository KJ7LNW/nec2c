#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([nec2c], [1.3], [nkcyham@yahoo.com])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([main.c])
AC_CONFIG_HEADERS([config.h])
AM_MAINTAINER_MODE

# --disable-optimizations sets reproducible-IEEE flags on both C and Fortran
# compilers; omitting the flag applies the standard -Wall -O2 defaults.
AC_ARG_ENABLE([optimizations],
  [AS_HELP_STRING([--disable-optimizations],
    [use reproducible-IEEE flags: -O0 -ffp-contract=off -fexcess-precision=standard])],
  [enable_opts=$enableval],
  [enable_opts=yes])

if test "x$enable_opts" = xno; then
  test "x${CFLAGS+set}" != xset && CFLAGS="-O0 -ffp-contract=off -fexcess-precision=standard"
  test "x${FFLAGS+set}" != xset && FFLAGS="-O0 -ffp-contract=off"
else
  test "x${CFLAGS+set}" != xset && CFLAGS="-Wall -O2"
  # nec2dx must not build at -O2: gfortran emits alignment-required SSE stores
  # in the large MAIN prologue, and stack alignment varies with argv, causing
  # SIGSEGV on helix-geometry decks.
  FFLAGS="-O0 -ffp-contract=off"
fi

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_F77
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Checks for libraries.
AC_CHECK_LIB([m], [sin])
AC_SUBST(PACKAGE_CFLAGS)
AC_SUBST(PACKAGE_LIBS)

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
