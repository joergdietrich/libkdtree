if OPTIMIZATION
AM_CFLAGS = -O3
AM_LDFLAGS = -O3
else
AM_CFLAGS = -O0
AM_LDFLAGS = -O0
endif
AM_CFLAGS += -W -Wall -Wmissing-prototypes -pedantic

# All of the code, internals included. The tests link against this.
noinst_LTLIBRARIES = libkdtree_core.la
libkdtree_core_la_SOURCES = kdtree_common.c kdtree_cartesian.c \
	kdtree_spherical.c kdtree.h kdtree_internal.h pqueue.c

# The installed library exports only the functions in kdtree.sym,
# which are those declared in kdtree.h.
lib_LTLIBRARIES = libkdtree.la
libkdtree_la_SOURCES =
libkdtree_la_LIBADD = libkdtree_core.la
libkdtree_la_LDFLAGS = $(AM_LDFLAGS) -version-info 1:0:0 \
	-export-symbols $(srcdir)/kdtree.sym
EXTRA_libkdtree_la_DEPENDENCIES = kdtree.sym

include_HEADERS = kdtree.h
EXTRA_DIST = kdtree.sym
