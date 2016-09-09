/*
 * main.cpp
 *
 *  Created on: Aug 3, 2014
 *      Author: dorian
 */

#include "smote.h"

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector wrapper(SEXP xs, SEXP target, SEXP nn, SEXP m, SEXP cor){

	//get dataset, target points, number of nearest neighbors and method to use
	const Rcpp::NumericMatrix dat =  Rcpp::NumericMatrix(xs);
	const Rcpp::NumericMatrix targ = Rcpp::NumericMatrix(target);
	const unsigned int k = Rcpp::as<int>(nn);
	const unsigned int method = Rcpp::as<int>(m);
	const unsigned int cores = Rcpp::as<int>(cor);

	Smote smote;
	return (Rcpp::wrap(smote.getKnn(dat, targ, k, method, 4)));

}

