/*! \mainpage libkdtree - A C99 implementation of the kd-tree algorithm
 *
 * \author Jörg Dietrich <jorgd@umich.edu>
 *
 * \section Introduction
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
 *   - Nearest neighbor search and
 *   - N-nearest neighbor search
 *   - Range search inside a (hyper-)sphere
 *   - Range search inside a (hyper-)rectangle
 *   - Sorted result lists are stored in an efficient min-max heap
 *   - Storage of arbitrary data associated with each point
 *
 */

