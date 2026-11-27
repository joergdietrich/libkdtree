/*! \mainpage libkdtree - A C99 implementation of the kd-tree algorithm
 *
 * \author Jörg Dietrich <astro@joergdietrich.com>
 *
 * Source code and issue tracker:
 * https://github.com/joergdietrich/libkdtree
 *
 * \section intro Introduction
 *
 * libkdtree provides a C99 implementation of kd-trees. Kd-trees are
 * space-partitioning trees that facilitate fast nearest-neighbor and
 * range searches. The present implementation contains two flavors of
 * kd-trees. The first kind of tree is constructed in n-dimensional
 * Euclidean space; the second kind of tree is constructed from data
 * on the 2-d surface of a sphere.
 *
 * Other features are:
 *   - Fast parallelized tree construction using POSIX threads
 *   - Nearest neighbor search
 *   - q-nearest neighbor search
 *   - Range search inside a (hyper-)sphere
 *   - Range search inside a (hyper-)rectangle
 *   - Sorted result lists are stored in an efficient min-max heap
 *   - Storage of arbitrary data associated with each point
 *
 * libkdtree is free software under the GNU Lesser General Public
 * License, version 2.
 *
 * \section install Building and installing
 *
 * Building needs autoconf, automake, libtool and a C99 compiler.
 *
 * \code{.sh}
 * git clone https://github.com/joergdietrich/libkdtree.git
 * cd libkdtree
 * autoreconf -fi
 * ./configure
 * make
 * make check
 * make install
 * \endcode
 *
 * By default only the static library is built; pass
 * <tt>--enable-shared</tt> to configure for a shared one. Programs
 * include kdtree.h and link with <tt>-lkdtree -lpthread -lm</tt>.
 *
 * \section usage Using the library
 *
 * \subsection building Building a tree
 *
 * kd_buildTree() builds a Cartesian tree of any dimension, and
 * kd_sph_buildTree() builds a tree on the sphere. Both take an array
 * of ::kd_point structures, which they reorder in place. The tree copies
 * the coordinates, so the array can be freed afterwards.
 *
 * - \c min and \c max are the corners of a box that holds all points.
 * - \c constr copies a point's data into the tree. Without it, the
 *   tree nodes carry no data. \c destr frees that copy when
 *   kd_destroyTree() destroys the tree.
 * - \c max_threads sets how many threads build the tree. Use a power
 *   of two; 1 builds in the calling thread.
 *
 * The tree needs at least one point. Searches only read the tree, so
 * several threads can search the same tree at once.
 *
 * \subsection radius Search radius
 *
 * Every search takes \c max_dist_sq, a pointer to the \e squared
 * search radius. Nearest-neighbor searches shrink it as they find
 * closer points, so reset it before each search.
 *
 * kd_nearest() and kd_sph_nearest() always return a node. If no
 * point lies within the starting radius, that node is not the
 * nearest neighbor. A radius larger than any possible distance
 * avoids this.
 *
 * \subsection spherical Spherical trees
 *
 * Coordinates are in radians. Element 0 is the longitude (right
 * ascension) in [0, 2&pi;], and element 1 is the latitude
 * (declination) in [-&pi;/2, &pi;/2]. Distances are great-circle
 * angles in radians.
 *
 * A rectangle for kd_sph_ortRangeSearch() that crosses longitude 0
 * starts at a negative longitude: -10&deg; to 10&deg; becomes
 * <tt>min[0] = -10&deg;</tt>, <tt>max[0] = 10&deg;</tt> (in
 * radians). Longitudes above &pi; then count as negative, so
 * <tt>max[0]</tt> must not exceed &pi;.
 *
 * \subsection results Result queues
 *
 * kd_qnearest(), kd_range() and kd_ortRangeSearch(), and their
 * spherical counterparts, return a ::pqueue, or NULL on error.
 *
 * - The results sit at <tt>res->d[1]</tt> to
 *   <tt>res->d[res->size - 1]</tt>, so there are
 *   <tt>res->size - 1</tt> of them. Each ::resItem points to its tree
 *   node and holds the squared distance.
 * - Unordered results (\c KD_UNORDERED and rectangle searches) have
 *   <tt>dist_sq = -1</tt>.
 * - pqremove_min() and pqremove_max() hand out ordered results by
 *   increasing or decreasing distance.
 * - To free a queue, free every ::resItem, then <tt>res->d</tt>, then
 *   \c res.
 *
 * The results point into the tree, so they are only valid until the
 * tree is destroyed.
 *
 * \section examples Examples
 *
 * Both examples are built with the library, from
 * <tt>src/nn_cartesian.c</tt> and <tt>src/nn_spherical.c</tt>.
 *
 * \subsection ex_cartesian Cartesian tree
 *
 * Nearest and 5 nearest neighbors of a point in 3-d space.
 *
 * \include nn_cartesian.c
 *
 * \subsection ex_spherical Spherical tree
 *
 * Nearest neighbor and all neighbors within 10&deg; of a point on the
 * sphere.
 *
 * \include nn_spherical.c
 *
 * The longer programs <tt>src/example_cartesian.c</tt> and
 * <tt>src/example_spherical.c</tt> run every kind of search and time
 * them.
 */
