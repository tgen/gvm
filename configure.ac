#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([gvm], [0.1], [skulkarni@tgen.org])
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([src])

: ${CFLAGS=""}

# Checks for programs.
AC_PROG_CC

AC_CHECK_PROG(GGO_PROG,gengetopt,yes)
if test x"${GGO_PROG}" != x"yes"
then AC_MSG_ERROR([Please install GNU gengetopt])
fi


# Checks for libraries.
AC_CHECK_LIB([hts], [bam_endpos], [], AC_MSG_ERROR(please install htslib))
AC_CHECK_LIB([hts], [bam_endpos2], [], AC_MSG_WARN(you should consider using the bam_endpos2 patch))
AC_CHECK_LIB([yaml], [yaml_parser_initialize], [], AC_MSG_ERROR(please install libyaml))

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h inttypes.h limits.h stddef.h stdint.h stdlib.h string.h sys/param.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([getpagesize gettimeofday memmove memset strdup strerror strtol])

AC_CONFIG_FILES([Makefile src/Makefile src/gengetopt/Makefile])
AC_OUTPUT
