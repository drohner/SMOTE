/*
 * vptree.cpp
 *
 *  Created on: Aug 13, 2014
 *      Author: dorian
 */

#include "vptree.h"

double calcMedian(std::vector<double> data);
const point_t& selectVP(const std::vector<point_t> &points, int sampleSize);
void deleteNode(VPTreeNode const *node);

//Constructor
VPTree::VPTree(){
	root = nullptr;
}

//Destructor
VPTree::~VPTree(){
	deleteNode(root);
	delete root;
}

//Helpermethod for destructor
void deleteNode(VPTreeNode const  *node){
	if (node->right->isLeaf){
		delete node->right;
	}
	else{
		deleteNode(node->right);
		delete node->right;
	}
	if (node->left->isLeaf){
		delete node->left;
	}
	else{
		deleteNode(node->left);
		delete node->left;
	}

}

VPTreeNode * VPTree::buildVPTree(const std::vector<point_t> &data){
	//get number of points
	const unsigned int dataRows = data.size();

	//nullptr check
	if(dataRows == 0){
		return nullptr;
	}

	//create leaf
	if(dataRows == 1){
		return new VPTreeNode(nullptr, nullptr, data[0], 0.0, true);
	}

	//get number of dimensions
	const unsigned int dataCols = data[0].point.size();

	//get vp
	#ifdef __gnu_linux__
		srand(time(nullptr));
	#else
		srand(dataRows);
	#endif
	const point_t &vp = selectVP(data, (dataRows/10)+25);

	//calculate distances from vp to every other point and mendian
	std::vector<double> distances(dataRows);
	for(int i = 0; i < dataRows; ++i){
		distances[i] = (vp.point-data[i].point).norm();
	}
	//calculate median of distances
	const double median = calcMedian(distances);

	//partition the data
	unsigned int countL = 0;
	unsigned int countR = 0;
	for(unsigned int i = 0; i < dataRows; ++i){
		if(distances[i] < median){
			++countL;
		}
		else {
			++countR;
		}
	}

	std::vector<point_t> left(countL);
	std::vector<point_t> right(countR);

	countL = 0;
	countR = 0;

	for(unsigned int j = 0; j < dataRows; ++j){
		if(distances[j] < median){
			left[countL] = data[j];
			++countL;
		}
		else {
			right[countR] = data[j];
			++countR;
		}
	}

	//create new node with partioned data, vp and median
	VPTreeNode *node = new VPTreeNode(buildVPTree(left), buildVPTree(right), vp, median, false);
	return node;
}

const point_t& selectVP(const std::vector<point_t> &points, int sampleSize){

	const unsigned int size = points.size();
	//check sampleSize
	if(sampleSize > size){
		sampleSize = size;
	}
	if(size == 1){
		return points[0];
	}

	//get sample of points
	std::vector<const point_t*> sampleSetP(sampleSize);
	int random = 0;

	for(int i = 0; i < sampleSize; ++i){
		random = (rand() % (int)(size));
		sampleSetP[i] = &points[random];
	}

	//allocate return values
	double bestSpread = -1;
	point_t const * bestP;

	//for each point p in sampleSetP
	for (const point_t* p : sampleSetP){
		//get sample of points
		std::vector<const point_t*> sampleSetD(sampleSize);
		for(int i = 0; i< sampleSize; ++i){
			random = (rand() % (int)(size));
			sampleSetD[i] = &points[random];
		}

		//calculate the distances for each point in sampleSetD to p
		std::vector<double> distances(sampleSize);
		for(int i = 0; i < sampleSize; ++i){
			distances[i] = (p->point - sampleSetD[i]->point).norm();
		}

		//calculate the median of the distances
		const double median = calcMedian(distances);
		//calculate the mean of the distances
		const double sum = std::accumulate(distances.begin(), distances.end(), 0.0);
		const double mean = sum / distances.size();

		//variance
		double tmp = 0;
		for(int i = 0; i < sampleSize; ++i){
			tmp += std::pow((distances[i]-median)-mean, 2);
		}

		//update spread
		const double spread = tmp/(sampleSize-1);
		if(spread > bestSpread){
			bestSpread = spread;
			bestP = p;
		}
	}
	//return best vantage point
	return *bestP;
}

void VPTree::knnSearch(const Eigen::VectorXd &target, const VPTreeNode &start, const unsigned int k, Eigen::MatrixXi &indices, Eigen::MatrixXd &distances, const unsigned int position) {
	//get search paramters
	const double median = start.median;
	const double dist = (target - start.vp.point).norm();
	const double sigma = distances(position, k - 1);
	//If start  is a leaf, check for knn
	if (start.isLeaf) {
		if (dist < sigma) {
			//insert sorted into the given row
			for (int i = 0; i < k; ++i) {
				if (distances(position, i) >= dist) {
					for (int j = k - 1; j > i; --j) {
						indices(position, j) = indices(position, j - 1);
						distances(position, j) = distances(position, j - 1);
					}
					indices(position, i) = start.vp.index;
					distances(position, i) = dist;
					return;
				}

			}
		}
	}
	else {
		//based on the construction skip specified subtrees
		if (dist < start.median) {
			if (dist < start.median + sigma) {
				knnSearch(target, *start.left, k, indices, distances, position);
			}
			if (dist >= start.median - sigma) {
				knnSearch(target, *start.right, k, indices, distances, position);
			}
		} else {
			if (dist >= start.median - sigma) {
				knnSearch(target, *start.right, k, indices, distances, position);
			}
			if (dist < start.median + sigma) {
				knnSearch(target, *start.left, k, indices, distances, position);
			}
		}

	}
}

void VPTree::parallelKnn(const Eigen::MatrixXd &target, const VPTreeNode &start, const unsigned int k, Eigen::MatrixXi &indices, Eigen::MatrixXd &distances, const unsigned int startindex, const unsigned int stopindex){
	for(int i = startindex; i < stopindex; ++i){
		knnSearch(target.row(i), start,k,indices, distances, i);
	}
}

//calculates the median based on std::nth_element
double calcMedian(std::vector<double> data){
	const unsigned int length = data.size()/2;
	std::nth_element(data.begin(), data.begin()+length, data.end());
	return data[length];
}


