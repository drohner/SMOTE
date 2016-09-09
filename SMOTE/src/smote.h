/*
 * smote.h
 *
 *  Created on: Sep 16, 2014
 *      Author: dorian
 */

#ifndef SMOTE_H_
#define SMOTE_H_

#include <Rcpp.h>

#include "lib/lsh/lshbox.h"

class Smote{

	public:
		/**
		 * Get the k nearest neighbors of a given target matrix for a given dataset. Different options can be used.
		 * @param data datamatrix, each row is one point
		 * @param target targetmatrix, each row is one point
		 * @param k number of nearest neighbor to compute
		 * @param method which datastructure should be used. 1: kd-tree 2: vp-tree 3: lsh 4: brute force
		 * @param cores the number of processor cores to use for paralleilitation
		 */
		Rcpp::NumericMatrix getKnn(const Rcpp::NumericMatrix &data,const Rcpp::NumericMatrix &target, const unsigned int k, const unsigned int method, const unsigned int cores);

		/**
		 * Get the k nearest neighbors of a given target matrix for a given dataset, using vp-trees
		 * @param data datamatrix, each row is one point
		 * @param target targetmatrix, each row is one point
		 * @param k number of nearest neighbor to compute
		 * @param cores the number of processor cores to use for paralleilitation
		 */
		Rcpp::NumericMatrix VPKnn(const Rcpp::NumericMatrix &data, const Rcpp::NumericMatrix &target, const unsigned int k, const unsigned int cores);

		/**
		 * Get the k nearest neighbors of a given target matrix for a given dataset, using lsh
		 * @param data datamatrix, each row is one point
		 * @param target targetmatrix, each row is one point
		 * @param k number of nearest neighbor to compute
		 * @param cores the number of processor cores to use for paralleilitation
		 */
		Rcpp::NumericMatrix LSHKnn(const Rcpp::NumericMatrix &data, const Rcpp::NumericMatrix &target, const unsigned int k, const unsigned int cores);

		/**
		 * Get the k nearest neighbors of a given target matrix for a given dataset, using either kd-trees oder lineare search
		 * @param data datamatrix, each row is one point
		 * @param target targetmatrix, each row is one point
		 * @param k number of nearest neighbor to compute
		 * @param method: 4: linear search, else kd-tree
		 * @param cores the number of processor cores to use for paralleilitation
		 */
		Rcpp::NumericMatrix FlannKnn(const Rcpp::NumericMatrix &data, const Rcpp::NumericMatrix &target, const unsigned int k, const unsigned int method, const unsigned int cores);

		/**
		 * Helper method for parallel lsh search
		 * @param tat the points in a matrix, for whom the nearst neighbors shall be computed
		 * @param k number of nearest neighbor to compute
		 * @param mylsh the lsh structure based on the complete dataset
		 * @param result the structure to save the indices of the nearest neighbors
		 * @param scanner scanner object for querying the nearst neighbors
		 * @param startindex start value for each thread in tat
		 * @param stopindex	end value for each thread in tat
		 */
		void parallelLSH(lshbox::Matrix<double> &tat, const unsigned int k, lshbox::itqLsh<double> &mylsh, std::vector<std::vector<std::pair<unsigned, float>>> &result, lshbox::Scanner<lshbox::Matrix<double>::Accessor> scanner, const unsigned int startindex,const unsigned int stopindex);

};

#endif /* SMOTE_H_ */
