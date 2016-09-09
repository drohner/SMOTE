\name{smote.concatinated}
\alias{smote.concatinated}

\title{
SMOTE
}
\description{
This function provides an implementation of SMOTE (Synthetic minority Oversampling TEchnique). This procedure was proposed by Chawla et. al. in 2002. Based on the k nearest neighbors new synthetic instances for the minority class are generated. This call should be used if your data is concatinated in one dataset.
}
\usage{
smote.concatinated(data, col, char, k = 5, o = 1, method = "e", seed = 321, cores = 2)
}
\arguments{
  \item{data}{
    A data.frame or numeric matrix containing all datapoints. One row describes a point, the columns are the different dimensions/features. The identification of the minority class will work with a given column and identificator.
}
  \item{col}{
    The column in which the identificator for the minority class can be found. This will be used numerical, so the name of the column won't work. See details or examples for further information.
}
  \item{char}{
    The character depicting the minority class in data. 
}
  \item{k}{
    How many nearest neighbors should be considered in SMOTE. Has to be an integer greater 0. Default is set to 5.
}
  \item{o}{
    How many times the method shall be used. Has to be either an integer greater 0 or a number between 0 and 1.
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
  This implementation will generate o * |minority| (number of points in minority) new samples. \cr
  The procedure does not work for categorial variables. \cr
  If the given input is invalid, an error will occur with further information. \cr
  To depict col use the numerical identifier for columns. For example, the column containing char is the fifth, col has to be 5. Using the name of the column will result in an error.\cr
  For further information on the datastructures for the k nearest neighbor search see \code{\link{knnsearch}}. \cr
  Due to datasets with a large number of dimensions and/or samples, two problems can occur: If the dataset is too big, the search can crash, because of limitations of the main memory. In this case, please clean the main memory and terminate applications, which are currently not necessary. Another problem can be the long elapsed time. In this case, please be patient, during the search no output is given.  Because of the high computing duration, the proper usage of the parallelization is highly recommended. However, parallelization will only work on unix-line operating systems.
}
\value{
  The function returns a matrix containing the synthetic instances. One row describes a point, the columns are the different dimensions/features. The given structure by 'col' and 'char' will be rebuilt. In some cases you have to set the names or the class of the synthetic instances. 
}
\references{
  Chawla, N. et. al. (2002). Smote: Synthetic minority over-sampling technique \cr
  Bunkhumpornpat, C. et. al. (2009). Safe-Level-SMOTE: Safe-Level-Synthetic Minority Over-Sampling TEchnique for Handling the Class Imbalanced Problem \cr
  Maciejewski, T. et. al. (2011). Local Neighbourhood Extension of SMOTE for Mining Imbalanced Data \cr
  See \code{\link{knnsearch}} for the used literature and references for the k nearest neighbor search.
}
\author{
  Dorian Rohner
}


\seealso{
\code{\link{smote.separated}},
\code{\link{lnsmote.concatinated}},
\code{\link{safelevelsmote.concatinated}},
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

#add labels
data = data.frame(data, '0')
minority.class1 = data.frame(minority.class1, '1')
minority.class2 = data.frame(minority.class2, '2')

#concatinate data
names(minority.class1) = names(data)
names(minority.class2) = names(data)
complete.data = rbind(data, minority.class1,  minority.class2)

#compute synthetic data for minority.class1
synthetic = smote.concatinated(complete.data, col = 3, char = '1', k = 5, o = 4, method = 'e')
}
