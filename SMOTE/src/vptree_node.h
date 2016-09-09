/*
 * vptree_node.h
 *
 *  Created on: Aug 6, 2014
 *      Author: dorian
 */

#ifndef VPTREE_NODE_H_
#define VPTREE_NODE_H_

#include "eigen3/Eigen/Dense"

#include <vector>

/**
 * struct for handling the datasets
 * point: d dimesional double vector containing a single point
 * index: index of point in the dataset
 */
struct point_t{
	Eigen::VectorXd point;
	int index;
};


class VPTreeNode{
	public:
		/**
		 * constructor for a vp-tree node
		 * @param left pointer to the left subtree
		 * @param right pointer to the right subtree
		 * @param vp the choosen vantage point
		 * @param median the median of the distances to the other points in this level
		 * @param isLeaf bool, true is this node is a leaf. Otherwise false
		 */
		VPTreeNode(const VPTreeNode *left,const VPTreeNode *right, const point_t vp, const double median, const bool isLeaf);

		const VPTreeNode *left;
		const VPTreeNode *right;
		const point_t vp;
		const double median;
		const bool isLeaf;
};



#endif /* VPTREE_NODE_H_ */
