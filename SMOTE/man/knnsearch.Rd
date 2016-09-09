\name{knnsearch}
\alias{knnsearch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  k nearest neighbor search
}
\description{
  This function provides several methods to compute the k nearest neighbors for a given dataset.
}
\usage{
knnsearch(data, target = data, k, method = "e", cores = 2)
}

\arguments{
  \item{data}{
    A data.frame or numeric matrix containing the complete dataset. One row describes a point, the columns are the different dimensions/features.
}
  \item{target}{
    A data.frame or numeric matrix containing the points, for whom the k nearest neighbors shall be computed. One row describes a point, the columns are the different dimensions/features. Has to be a subset of data. Default is set to data.
}
  \item{k}{
    How many nearest neighbors shall be computed. Has to be an integer greater 0. Default is set to 5.
}
  \item{method}{
    The method to use for the k nearest neighbor search in the euclidean space. Currently there are four approaches implemented. You can determine by yourself which to use or leave it up to the implementation. The call parameters are
    \itemize{
      \item k for k-d tree: common space partitioning data structure for exact search. Works well for low dimensions.
      \item v for vp tree: vantage point trees for exact nearest neighbor search. Currently not recommended for usage.
      \item b for brute force: for each point the distance to every other point is calculated and afterwards the neighbors are evaluated. Works well for high dimensions, but bad for a high number of points.
      \item l for locality sensitive hashing (lsh): an approximative procedure. The summarized distances are 0.1-0.15 off. Works great in high dimensions.
      \item a for approximative: the program determines which method should be used for the fastest search. lsh can be one of the methods.
      \item e for exact: the program determines which method should be used for the fastest search. Only exact methods will be used.
    }
    Default is set to 'e'.    
}
  \item{cores}{
    How many cores of your processor shall be used for parallelization. Has to be an integer greater 0. Default is set to 2.
}
}
\details{
  If the given input is invalid, an error will occur with further information. \cr
  The vp trees are currently not recommended, because, in most cases, they are slower than the k-d trees or brute force. If there is no special reason to select a search option explicit, leave it up to the implementation. It was carefully evaluated which method is the fastest. \cr
  The k-d tree and the linear seach are provided by the C++ library FLANN. For the lsh method the github project lshbox is used. These libraries were wrapped using the rcpp package. \cr
  Due to datasets with a large number of dimensions and/or samples, two problems can occur: If the dataset is too big, the search can crash, because of limitations of the main memory. In this case, please clean the main memory and terminate applications, which are currently not necessary. Another problem can be the long elapsed time. In this case, please be patient, during the search no output is given. Because of the high computing duration, the proper usage of the parallelization is highly recommended. However, parallelization will only work on unix-like operating systems.
}
\value{
  This function returns a numeric matrix, containing the indices of the nearest neighbors in data. The indices are sorted ascending by the distance. The nearest neighbor is always the index of the point itself.
}
\references{
  FLANN: Marius Muja and David G. Lowe (2013). Fast Library for Approximate Nearest Neighbors. http://www.cs.ubc.ca/research/flann/ \cr
  LSHBOX: Huaxia Wang et. al. (2014). A C++ Toolbox of Locality-Sensitive Hashing for Large Scale Image Retrieval. https://github.com/RSIA-LIESMARS-WHU/LSHBOX and Yunchao Gong and Lazebnik, S. (2011). Iterative Quantization: A Procrustean Approach to Learning Binary Codes\cr
  EIGEN v3: Gael Guennebaud and Benoit Jacob and others (2010), http://eigen.tuxfamily.org/ \cr
  vp tree: Peter N. Yianilos (1993). Data Structures and Algorithms for Nearest Neighbor Search in General Metric Spaces and Ada Wai-chee Fu et. al. (2000). Dynamic vp-tree indexing for n-nearest neighbor search given pair-wise distances
}
\author{
  Dorian Rohner
}

\examples{
#generate data
data = replicate(10, runif(1000))
#get subset of data
target = data[1:100,]

#search nearest neigbors for target in data
knn.target = knnsearch(data = data, target = target, k = 5, method = 'e', cores = 2)

#search nearest neighbors for data
knn.data = knnsearch(data = data, k = 5) #method will be 'e' and number of cores = 2
}

