AM_CPPFLAGS = -W -Wall -I$(top_srcdir)/lib -I$(top_srcdir)/tests
# Most tests call internal functions, so they link against the
# library's code as separate files. The tests that use only kdtree.h
# link against the installed library. The tests call pthread and libm
# functions themselves.
LDADD = ../lib/libkdtree_core.la -lpthread -lm
PUBLIC_LDADD = ../lib/libkdtree.la -lpthread -lm

EXTRA_DIST = test_cartesian_nn.sh \
	test_cartesian_rectangle.sh \
	test_cartesian_range.sh \
	test_cartesian_qnearest.sh \
	test_spherical_nn.sh \
	test_spherical_rectangle.sh \
	test_spherical_range.sh \
	test_spherical_qnearest.sh

TESTS = test_kd_sqr \
	test_duplicates \
	test_build_errors \
	test_edge_cases \
	test_kd_sph_dist \
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
	test_duplicates \
	test_build_errors \
	test_edge_cases \
	test_kd_sph_dist \
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
test_duplicates_LDADD = $(PUBLIC_LDADD)
test_duplicates_SOURCES = test_duplicates.c
test_build_errors_LDADD = $(PUBLIC_LDADD)
test_build_errors_SOURCES = test_build_errors.c
test_edge_cases_SOURCES = test_edge_cases.c
test_kd_sph_dist_LDADD = $(PUBLIC_LDADD)
test_kd_sph_dist_SOURCES = test_kd_sph_dist.c
test_kd_sph_isPointInRect_SOURCES = test_kd_sph_isPointInRect.c
test_kd_sph_isRectInRect_SOURCES = test_kd_sph_isRectInRect.c 
test_cartesian_nn_LDADD = $(PUBLIC_LDADD)
test_cartesian_nn_SOURCES = test_cartesian_nn.c create_pointlist.c \
	points_eq.c kdtree_test.h
test_cartesian_rectangle_SOURCES = test_cartesian_rectangle.c \
	create_pointlist.c sorted_queues_eq.c unsorted_queues_eq.c \
	naive_ortRangeSearch.c naive_item.c
test_cartesian_qnearest_SOURCES = test_cartesian_qnearest.c \
	create_pointlist.c sorted_queues_eq.c naive_item.c
test_cartesian_range_SOURCES = test_cartesian_range.c create_pointlist.c \
	sorted_queues_eq.c naive_item.c
test_spherical_nn_LDADD = $(PUBLIC_LDADD)
test_spherical_nn_SOURCES = test_spherical_nn.c create_pointlist.c points_eq.c
test_spherical_rectangle_SOURCES = test_spherical_rectangle.c \
	create_pointlist.c sorted_queues_eq.c unsorted_queues_eq.c \
	naive_sph_ortRangeSearch.c naive_item.c
test_spherical_range_SOURCES = test_spherical_range.c create_pointlist.c \
	sorted_queues_eq.c naive_item.c
test_spherical_qnearest_SOURCES = test_spherical_qnearest.c \
	create_pointlist.c sorted_queues_eq.c naive_item.c
