if OPTIMIZATION
AM_CFLAGS = -O3
AM_LDFLAGS = -O3
else
AM_CFLAGS = -O0
AM_LDFLAGS = -O0
endif
AM_CFLAGS += -W -Wall -Wmissing-prototypes -pedantic

# The sources as separate files, with the internals in kdtree_internal.h
# visible. The tests link against this.
noinst_LTLIBRARIES = libkdtree_core.la
libkdtree_core_la_SOURCES = kdtree_common.c kdtree_cartesian.c \
	kdtree_spherical.c kdtree.h kdtree_internal.h pqueue.c

# The installed library, built from kdtree_all.c, which includes the
# same sources as one file with the internals static.
lib_LTLIBRARIES = libkdtree.la
libkdtree_la_SOURCES = kdtree_all.c
libkdtree_la_LIBADD = -lpthread -lm
libkdtree_la_LDFLAGS = $(AM_LDFLAGS) -version-info 1:0:0

include_HEADERS = kdtree.h
