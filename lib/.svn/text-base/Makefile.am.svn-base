if OPTIMIZATION
AM_CFLAGS = -O3
AM_LDFLAGS = -O3
else
AM_CFLAGS = -O0
AM_LDFLAGS = -O0
endif
AM_CFLAGS += -W -Wall -Wmissing-prototypes -pedantic

lib_LTLIBRARIES = libkdtree.la
libkdtree_la_SOURCES = kdtree_common.c kdtree_cartesian.c kdtree_spherical.c \
	kdtree.h pqueue.c pqueue.h qsort.c pmergesort.c

include_HEADERS=kdtree.h

EXTRA_DIST = COPYING.DJ