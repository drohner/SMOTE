/*
 * vptree.h
 *
 *  Created on: Aug 3, 2014
 *      Author: dorian
 */

#ifndef VPTREE_H_
#define VPTREE_H_

#include "vptree_node.h"


class VPTree{
	public:
		VPTreeNode *root;

		VPTree();
		~VPTree();

		/**
		 * Build a vp-tree on a vector of points (Eigenvector and index)
		 * @param data complete dataset as std::vector<point_t>
		 * @return  root of a vp-tree
		 */
		VPTreeNode * buildVPTree(const std::vector<point_t> &data);


		/**
		 * Search method based on a vp-tree, to find the k nearest neighbors
		 * @param target point, whose neighbors should be found
		 * @param start node, where to start the search. root for certain correct results
		 * @param k number of nearst neighbors to compute
		 * @param indices the matrix containing the indices of the nearest neighbors. will be updated during this procedure
		 * @param distances the matrix containing the distances of the nearest neighbors. will be updated during this procedure
		 * @param postion where the results shall be written in indices and distances
		 */
		void knnSearch(const Eigen::VectorXd &target, const VPTreeNode &start, const unsigned int k, Eigen::MatrixXi &indices, Eigen::MatrixXd &distances, const unsigned int position);

		/**
		 * helper method for parallel search
		 * @param target matrix containing the points for whom the nearst neighbors shall be computed
		 * @param start node, where to start the search. root for certain correct results
		 * @param k number of nearst neighbors to compute
		 * @param indices the matrix containing the indices of the nearest neighbors. will be updated during this procedure
		 * @param distances the matrix containing the distances of the nearest neighbors. will be updated during this procedure
		 * @param startindex where to start in target
		 * @param stopindex where to stop in target
		 */
		void parallelKnn(const Eigen::MatrixXd &target, const VPTreeNode &start, const unsigned int k, Eigen::MatrixXi &indices, Eigen::MatrixXd &distances, const unsigned int startindex, const unsigned int stopindex);

};

#endif /* VPTREE_H_ */
