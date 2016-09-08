# install
install.packages("SMOTE_0.2.tar.gz", repos = NULL, type = "source")
library(SMOTE)


?SMOTE
?lnsmote.separated
?knnsearch

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



