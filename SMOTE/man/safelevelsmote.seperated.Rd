\name{safelevelsmote.separated}
\alias{safelevelsmote.separated}

\title{
Safe-level-SMOTE
}
\description{
This function provides an implementation of Safe-Level-SMOTE (Synthetic minority Oversampling TEchnique). This procedure was proposed by Bunkhumpornpat et. al. in 2009. Based on the k nearest neighbors and the safe-level, new synthetic instances for the minority class are generated. This call should be used if your data is separated according to the classes.
}
\usage{
safelevelsmote.separated(data, minority, k = 5, method = "e", seed = 321, cores = 2)
}
\arguments{
  \item{data}{
    A data.frame or numeric matrix containing the complete dataset excluding the minority class. One row describes a point, the columns are the different dimensions/features.
  }
  \item{minority}{
    A data.frame or numeric matrix containing the datapoints of the minority class to oversample. One row describes a point, the columns are the different dimensions/features.
}
  \item{k}{
    How many nearest neighbors should be considered in SMOTE. Has to be an integer greater 0. Default is set to 5.
}
  \item{method}{
    The method to use for the k nearest neighbor search in the euclidean space. Currently there are four approaches implemented. You can determine by yourself which to use or leave it up to the implementation. The call parameters are
    \itemize{
      \item k for k-d tree: common space partitioning data structure for exact search. Works well for low dimensions.
      \item v for vp tree: vantage point trees for exact nearest neighbor search. Currently not recommended for usage.
      \item b for brute force: for each point the distance to every other point is calculated and afterwards the. neighbors are evaluated. Works well for high dimensions, but bad for a high number of points.
      \item l for locality sensitive hashing (lsh): an approximative procedure. The summarized distances are 0.1-0.15 off. Works great in high dimensions.
      \item a for approximative: the program determines which method should be used for the fastest search. lsh can be one of the methods.
      \item e for exact: the program determines which method should be used for the fastest search. Only exact methods will be used.
    }
    Default is set to 'e'. For further information see \code{\link{knnsearch}}.
    
}
  \item{seed}{
    Sets the random generator. Has to be an integer greater 0. Default is set to 321.
}
  \item{cores}{
    How many cores of your processor shall be used for parallelization. Has to be an integer greater 0. Default is set to 2.
}
}
\details{
  This implementation will generate at most |minority| (number of points in minority) new samples. \cr
  The procedure does not work for categorial variables. \cr
  If the given input is invalid, an error will occur with further information. \cr
  For further information on the datastructures for the k nearest neighbor search see \code{\link{knnsearch}}. \cr
  Due to datasets with a large number of dimensions and/or samples, two problems can occur: If the dataset is too big, the search can crash, because of limitations of the main memory. In this case, please clean the main memory and terminate applications, which are currently not necessary. Another problem can be the long elapsed time. In this case, please be patient, during the search no output is given.  Because of the high computing duration, the proper usage of the parallelization is highly recommended. However, parallelization will only work on unix-like operating systems.
}
\value{
  The function returns a matrix containing the synthetic instances. One row describes a point, the columns are the different dimensions/features. In some cases you have to set the names or the class of the synthetic instances. 
}
\references{
  Chawla, N. et. al. (2002). Smote: Synthetic minority over-sampling technique \cr
  Bunkhumpornpat, C. et. al. (2009). Safe-Level-SMOTE: Safe-Level-Synthetic Minority Over-Sampling TEchnique for Handling the Class Imbalanced Problem \cr
  See \code{\link{knnsearch}} for the used literature and references for the k nearest neighbors search.
}
\author{
  Dorian Rohner
}

\seealso{
  \code{\link{safelevelsmote.concatinated}},
  \code{\link{smote.separated}},
  \code{\link{lnsmote.separated}},
  \code{\link{knnsearch}}
}
\examples{
#generate data with 2 features
#majority class: 250 samples
data = replicate(2, runif(250))
#first minority class: 50 samples
minority.class1 = replicate(2, runif(50))
#second minority class: 50 samples
minority.class2 = replicate(2, runif(50))

#compute synthetic data for minority.class1
data = rbind(data, minority.class2)
synthetic = safelevelsmote.separated(data = data, minority = minority.class1, k = 5, method = 'e')
}

