\name{SMOTE-package}
\alias{SMOTE-package}
\alias{SMOTE}
\docType{package}
\title{
Implementation of SMOTE algorithms
}
\description{
 In this package three algorithms for oversampling are implemented. These are SMOTE (Synthetic minority over-sampling technique), Safe-Level-SMOTE and Local Neighbourhood SMOTE (LNSMOTE). For these procedures the k nearest neighbors have to be computed. Therefore, a method wrapping several C++ libraries providing this functionality is implemented.
}
\details{
\tabular{ll}{
Package: \tab SMOTE\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2014-11-02\cr
License: \tab GNU GPLv3\cr
}

}
\author{
Dorian Rohner s1dorohn@stmail.uni-bayreuth.de \cr
This package was developed during a bachelor thesis at the University of Bayreuth under supervison of Dr. Christina Bogner (Departement of Ecological Modelling) and Maximilian Sand (Chair for Applied Computer Science III: Robotics and Embedded Systems)
}
\references{
  Chawla, N. et. al. (2002). Smote: Synthetic minority over-sampling technique \cr
  Bunkhumpornpat, C. et. al. (2009). Safe-Level-SMOTE: Safe-Level-Synthetic Minority Over-Sampling TEchnique for Handling the Class Imbalanced Problem \cr
  Maciejewski, T. et. al. (2011). Local Neighbourhood Extension of SMOTE for Mining Imbalanced Data \cr
  See \code{\link{knnsearch}} for the used literature and references for the k nearest neighbor search.
}

\keyword{ package }
\seealso{
\code{\link{smote.separated}}, \cr
\code{\link{smote.concatinated}}, \cr
\code{\link{safelevelsmote.separated}}, \cr
\code{\link{safelevelsmote.concatinated}}, \cr
\code{\link{lnsmote.separated}}, \cr
\code{\link{lnsmote.concatinated}}, \cr
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

data = rbind(data, minority.class2)

#compute synthetic data for minority.class1
synthetic = smote.separated(minority = minority.class1, k = 5, o = 4, method = 'e')

synthetic = safelevelsmote.separated(data = data, minority = minority.class1, k = 5, method = 'e')

synthetic = lnsmote.separated(data = data, minority = minority.class1, k = 5, o = 4, method = 'e')
}