/*
 * vptree_node.cpp
 *
 *  Created on: Aug 14, 2014
 *      Author: dorian
 */

#include "vptree_node.h"

VPTreeNode::VPTreeNode(const VPTreeNode *left,const VPTreeNode *right, const point_t vp, const double median, const bool isLeaf) :
	left(left),
	right(right),
	vp(vp),
	median(median),
	isLeaf(isLeaf)
	{}
