/*
 * smote.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: dorian
 */

#include "smote.h"
#include "vptree.h"

#include "lib/flann/flann.hpp"

#ifdef __gnu_linux__
	#include <thread>
#endif


Rcpp::NumericMatrix Smote::getKnn(const Rcpp::NumericMatrix &data, const Rcpp::NumericMatrix &target,const unsigned int k, const unsigned int method, const unsigned int cores){
	switch (method) {
		case (1):{
			//Flann: KD
			return FlannKnn(data, target, k, method, cores);

		case (2):
			//VP
			return VPKnn(data, target, k, cores);

		case (3):
			//LSH
			return LSHKnn(data, target, k, cores);

		case (4):
			//Flann: Linear Search
			return FlannKnn(data, target, k, method, cores);

		default:
			//Default: KD
			return FlannKnn(data, target, k, method, cores);
		}
	}
}

/*
 * VP Tree
 */
Rcpp::NumericMatrix Smote::VPKnn(const Rcpp::NumericMatrix &data, const Rcpp::NumericMatrix &target, const unsigned int k, const unsigned int cores) {
	//get dimensions
	const unsigned int dataRows = data.nrow();
	const unsigned int dataCols = data.ncol();
	const unsigned int targetRows = target.nrow();
	const unsigned int targetCols = target.ncol();

	//copy data to eigen
	Eigen::MatrixXd eigendat(dataRows, dataCols);
	for (int i = 0; i < dataRows; i++) {
		for (int j = 0; j < dataCols; j++) {
			eigendat(i, j) = data(i, j);
		}
	}

	//create vector of points, paired with index
	std::vector<point_t> points(dataRows);
	for (int i = 0; i < dataRows; i++) {
		point_t tmp;
		tmp.point = eigendat.row(i);
		tmp.index = i;
		points[i] = tmp;
	}

	//copy target to eigen
	Eigen::MatrixXd eigentarget(targetRows, targetCols);
	for (int i = 0; i < targetRows; i++) {
		for (int j = 0; j < targetCols; j++) {
			eigentarget(i, j) = target(i, j);
		}
	}

	//create VP Tree
	VPTree vptree;

	vptree.root = vptree.buildVPTree(points);

	//allocate result matrix
	Eigen::MatrixXi indices (targetRows, k);

	//allocate distance matrix
	Eigen::MatrixXd distances = Eigen::MatrixXd::Constant(targetRows, k, INFINITY);


	//will be alawyas used on windows, because cores is set to 1 in R
	if (cores == 1) {
		//search for each point of target
		for (int i = 0; i < targetRows; i++) {
			vptree.knnSearch(eigentarget.row(i), *vptree.root, k, indices, distances, i);
		}
	}
	//parallel search
	else{
		#ifdef __gnu_linux__
		const unsigned int border = targetRows / cores;
		std::thread threads[cores];
		for (int i = 0; i < cores; ++i) {
			threads[i] = std::thread(&VPTree::parallelKnn, &vptree, std::ref(eigentarget), std::ref(*vptree.root), k, std::ref(indices), std::ref(distances),  i * border, (i + 1) * border);
		}
		for (int j = 0; j < cores; ++j) {
			threads[j].join();
		}
		#endif
	}

	//copy results to return matrix
	Rcpp::NumericMatrix ret(targetRows, k);
	for (int i = 0; i < targetRows; ++i) {
		for (int j = 0; j < k; ++j) {
			ret(i, j) = indices(i, j);
		}
	}
	return ret;
}


/*
 * Locality Sensitive Hashing
 * https://github.com/RSIA-LIESMARS-WHU/LSHBOX
 */
Rcpp::NumericMatrix Smote::LSHKnn(const Rcpp::NumericMatrix &data, const Rcpp::NumericMatrix &target, const unsigned int k, const unsigned int cores){
	//get dimensions
	const unsigned int dataRows = data.nrow();
	const unsigned int dataCols = data.ncol();
	const unsigned int targetRows = target.nrow();
	const unsigned int targetCols = target.ncol();

	//copy data
	std::vector<double> vecdat(dataRows * dataCols);
	for (int i = 0; i < dataRows; ++i) {
		for (int j = 0; j < dataCols; ++j) {
			vecdat[i * dataCols + j] = data(i, j);
		}
	}

	//copy target
	std::vector<double> vectat(targetRows * targetCols);
	for (int i = 0; i < targetRows; ++i) {
		for (int j = 0; j < targetCols; ++j) {
			vectat[i * targetCols + j] = target(i, j);
		}
	}

	//load data and target into lsf matrix format
	lshbox::Matrix<double> mat;
	lshbox::Matrix<double> tat;

	mat.load(vecdat, dataRows, dataCols);
	tat.load(vectat, targetRows, targetCols);

	//initialize lsh method (interative quantization)
	lshbox::itqLsh<double> mylsh;

	//LSH Parameters
	lshbox::itqLsh<double>::Parameter param;

	param.M = 300 * (dataRows / 10000) + k;
	param.L = dataCols / 25 + 1;
	param.D = mat.getDim();
	param.N = 4 * (dataRows / 1000) + k + 1;
	if (dataCols <= param.N) {
		param.N = dataCols;
	}
	param.S = dataRows * 0.01 + 100;
	param.I = 50 + dataRows * 0.00005;

	//train dataset
	mylsh.reset(param);
	mylsh.train(mat);

	//generate k nearest neighbor scanner
	lshbox::Matrix<double>::Accessor accessor(mat);
	lshbox::Metric<double> metric(targetCols, L2_DIST);
	lshbox::Scanner<lshbox::Matrix<double>::Accessor> scanner(accessor, metric,	k, std::numeric_limits<float>::max());

	//allocate return matrix
	Rcpp::NumericMatrix ret(targetRows, k);

	//will be alawyas used on windows, because cores is set to 1 in R
	if (cores == 1) {
		for (int i = 0; i < targetRows; ++i) {
			//Query the nearest neighbors
			scanner.reset(tat[i]);
			mylsh.query(tat[i], scanner);

			//get vector of indices paired with distacnes
			std::vector<std::pair<unsigned, float>> result;
			result = scanner.topk().getTopk();

			//copy to result
			for (int j = 0; j < result.size(); ++j) {
				ret(i, j) = result[j].first;
			}
		}
	}
	//parallel search
	else {
		#ifdef __gnu_linux__
		int border = targetRows / cores;

		//allocate temporal storage for the paired indices and distances
		std::vector<std::vector<std::pair<unsigned, float>>>result(targetRows);

		//generate threads, start and join them
		std::thread threads[cores];
		for (int i = 0; i < cores; ++i) {
			threads[i] = std::thread(&Smote::parallelLSH, this, std::ref(tat), k, std::ref(mylsh), std::ref(result), scanner, i * border, (i + 1) * border);
		}
		for (int j = 0; j < cores; ++j) {
			threads[j].join();
		}

		//copy indices to return matrix
		for (int i = 0; i < targetRows; ++i) {
			for (int j = 0; j < k; ++j) {
				ret(i, j) = result[i][j].first;
			}
		}
		#endif
	}
	return ret;
}


void Smote::parallelLSH(lshbox::Matrix<double> &tat, const unsigned int k, lshbox::itqLsh<double> &mylsh, std::vector<std::vector<std::pair<unsigned, float>>> &result, lshbox::Scanner<lshbox::Matrix<double>::Accessor> scanner, const unsigned int startindex,const unsigned int stopindex){
	//from start to end search based on the lsh structure the k nearest neighbor using a scanner object
	for(int i = startindex; i < stopindex; ++i){
		 scanner.reset(tat[i]);
		 mylsh.query(tat[i], scanner);
		 result[i] = scanner.topk().getTopk();
	}
}

/*
 * FLANN
 * http://www.cs.ubc.ca/research/flann/
 * http://opencv.org/
 */
Rcpp::NumericMatrix Smote::FlannKnn(const Rcpp::NumericMatrix &data, const Rcpp::NumericMatrix &target, const unsigned int k, const unsigned int method, const unsigned int cores){
	//get dimensions
	const unsigned int dataRows = data.nrow();
	const unsigned int dataCols = data.ncol();
	const unsigned int targetRows = target.nrow();
	const unsigned int targetCols = target.ncol();

	//copy data
	flann::Matrix<double> flannMat(new double[dataRows * dataCols], dataRows, dataCols);
	for (int i = 0; i < dataRows; i++) {
		for (int j = 0; j < dataCols; j++) {
			flannMat[i][j] = data(i, j);
		}
	}

	//copy target
	flann::Matrix<double> flannTarget(new double[targetRows * targetCols], targetRows, targetCols);
	for (int i = 0; i < targetRows; i++) {
		for (int j = 0; j < targetCols; j++) {
			flannTarget[i][j] = target(i, j);
		}
	}

	//allocate indices matrix
	flann::Matrix<int> indices(new int[targetRows * k], targetRows, k);
	//allocate distance matrix
	flann::Matrix<double> distances(new double[targetRows * k], targetRows, k);

	if(method == 4){
		//Bruteforce
		//Build index
		flann::Index<flann::L2<double>> index(flannMat, flann::LinearIndexParams());
		index.buildIndex();

		//Specify search parameters
		flann::SearchParams params(flann::FLANN_CHECKS_UNLIMITED);
		params.cores = cores;

		//execute the search
		index.knnSearch(flannTarget, indices, distances, k, params);
	}
	else{
		//kd tree
		//Build index. For exact search only one kd tree is reasonable
		flann::Index<flann::L2<double>> index(flannMat, flann::KDTreeIndexParams(1));
		index.buildIndex();

		//Specify search parameters
		flann::SearchParams params(flann::FLANN_CHECKS_UNLIMITED);
		params.use_heap = flann::FLANN_False;
		params.cores = cores;

		//execute the search
		index.knnSearch(flannTarget, indices, distances, k, params);
	}

	//allocate return matrix
	Rcpp::NumericMatrix ret(targetRows, k);

	//copy data
	for (int i = 0; i < targetRows; i++) {
		for (int j = 0; j < k; j++) {
			ret(i, j) = indices[i][j];
		}
	}

	//clear memory
	delete[] flannMat.ptr();
	delete[] flannTarget.ptr();
	delete[] indices.ptr();
	delete[] distances.ptr();

	return ret;
}
