# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(libkdtree, 0.4.3, astro@joergdietrich.com)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([lib/kdtree_common.c])
AC_CONFIG_HEADER([config.h])

# Optional features
AC_MSG_CHECKING(whether to enable compile time optimization)
AC_ARG_ENABLE([optimization], [AC_HELP_STRING([--enable-optimization],
			      [enable compile time optimization [default=yes] ])],
			       if test $enableval = no; then
			       	       [AC_MSG_RESULT(no)]
			               DO_OPT=no
			       fi
			       ,
			       [AC_MSG_RESULT(yes)]
			       DO_OPT=yes
			       )
AM_CONDITIONAL([OPTIMIZATION], [test "$DO_OPT" = yes])

# Checks for programs and avoid the automatic addition of the -O2 flag
# by autoconf by setting CFLAGS is not set yet.
if test -z "${CFLAGS}"; then
   CFLAGS="-g"
fi
AC_PROG_CC_C99
AC_DISABLE_SHARED
AC_PROG_LIBTOOL


# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdint.h stdlib.h string.h], [], [exit 1])
AC_CHECK_HEADER(sys/time.h)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_HEADER_TIME
AC_TYPE_UINT32_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gettimeofday])

# Integrate Doxygen support
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF) 
DX_MAN_FEATURE(OFF) 
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN($PACKAGE_NAME, Doxyfile)

AC_CONFIG_FILES(Makefile \
	Doxyfile \
	lib/Makefile \
	src/Makefile \
	tests/Makefile)
AC_OUTPUT
