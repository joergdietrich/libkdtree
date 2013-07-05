AM_CPPFLAGS = -W -Wall -I$(top_srcdir)/lib -I$(top_srcdir)/tests
LDADD = -lkdtree -lpthread -lm
AM_LDFLAGS = -L../lib

EXTRA_DIST = test_cartesian_nn.sh \
	test_cartesian_rectangle.sh \
	test_cartesian_range.sh \
	test_cartesian_qnearest.sh \
	test_spherical_nn.sh \
	test_spherical_rectangle.sh \
	test_spherical_range.sh \
	test_spherical_qnearest.sh

TESTS = test_kd_sqr \
	test_kd_sph_isPointInRect \
	test_kd_sph_isRectInRect \
	test_cartesian_nn.sh \
	test_cartesian_rectangle.sh \
	test_cartesian_range.sh \
	test_cartesian_qnearest.sh \
	test_spherical_nn.sh \
	test_spherical_rectangle.sh \
	test_spherical_range.sh \
	test_spherical_qnearest.sh

check_PROGRAMS = test_kd_sqr \
	test_kd_sph_isPointInRect \
	test_kd_sph_isRectInRect \
	test_cartesian_nn \
	test_cartesian_rectangle \
	test_cartesian_qnearest \
	test_cartesian_range \
	test_spherical_nn \
	test_spherical_rectangle \
	test_spherical_range \
	test_spherical_qnearest

test_kd_sqr_SOURCES =  test_kd_sqr.c
test_kd_sph_isPointInRect_SOURCES = test_kd_sph_isPointInRect.c
test_kd_sph_isRectInRect_SOURCES = test_kd_sph_isRectInRect.c 
test_cartesian_nn_SOURCES = test_cartesian_nn.c create_pointlist.c \
	points_eq.c kdtree_test.h
test_cartesian_rectangle_SOURCES = test_cartesian_rectangle.c \
	create_pointlist.c sorted_queues_eq.c unsorted_queues_eq.c \
	naive_ortRangeSearch.c
test_cartesian_qnearest_SOURCES = test_cartesian_qnearest.c \
	create_pointlist.c sorted_queues_eq.c
test_cartesian_range_SOURCES = test_cartesian_range.c create_pointlist.c \
	sorted_queues_eq.c	
test_spherical_nn_SOURCES = test_spherical_nn.c create_pointlist.c points_eq.c
test_spherical_rectangle_SOURCES = test_spherical_rectangle.c \
	create_pointlist.c sorted_queues_eq.c unsorted_queues_eq.c \
	naive_sph_ortRangeSearch.c
test_spherical_range_SOURCES = test_spherical_range.c create_pointlist.c \
	sorted_queues_eq.c 
test_spherical_qnearest_SOURCES = test_spherical_qnearest.c \
	create_pointlist.c sorted_queues_eq.c
