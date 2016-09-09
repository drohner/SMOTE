safelevelsmote.separated <-
function(data, minority, k = 5, method = "e", seed = 321, cores = 2){
  #Check Input
  .checkinput_sl(nrow(data), k, seed)
  .checkinput_sl_sep(data, minority)
  
  #save names
  minority.colnames = names(minority)
  
  #Label every minority classmember in the last column with 1, other points with 0  
  minority = cbind(minority, 1)
  data = cbind(data, 0)  
  
  #concatinate data, so the indices of minority are 1-nrow(minority)
  names(minority) = names(data)
  data = rbind(minority, data)
  
  synthetic = .safelevelsmote(data, nrow(minority), k, method, seed, cores)
  
  if(!is.null(synthetic)){
    names(synthetic) = minority.colnames
  }
  return (synthetic)
}

# Function called by the facade methods
# Generates the synthetic instances by calling a lapply function
# @data: the concatinated datatset [matrix]
# @minority.rows: the number of rows of minority [matrix]
# @k: number of nearest neighbors to consider [int]
# @method: the datastructure to use for knn search [int]
# @seed: set the random generator [int]
# @cores: the number of cores to use for parallelization [int]
.safelevelsmote = function(data, minority.rows, k, method, seed, cores){
  #parallelization doesnt work for windows
  if(Sys.info()['sysname'] == 'Windows'){
    cores = 1
  }
  
  #get knn  
  method = .transformMethod(method, ncol(data))
  
  # k +1 because the search will return the point itself as the nearest neighbor. removed after the seach is finished
  knn = knnsearch(as.matrix(data[,-ncol(data)]) , as.matrix(data[,-ncol(data)]), k+1, method, cores)[,-1]
  
  set.seed(seed)  
  data = as.matrix(data)
  
  #call apply function over the indcies. because the minority class samples are the saved in 1:minority.rows, we can use this approach
  smotepoints = do.call(rbind, lapply(c(1:minority.rows), FUN = ".safelevelsmote_apply", data = data, knn = knn))
  
  smotepoints.nrow = 0
  if(!is.null(smotepoints)){
    smotepoints.nrow = nrow(smotepoints)
  }
  
  cat("The minority class got", minority.rows ,"samples.\n")
  cat("Theoretically",minority.rows,"new samples could have been generated.\n")
  cat("Due to noisy samples, only", smotepoints.nrow,"synthetic instances were generated.\n")
  return (smotepoints)
}

# Calculates the synthetic instances based on safe level smote
# This funtion ist thought to be used in an lapply function
# @index: the index of the minority class sample in data [int]
# @data: the complete dataset [matrix]
# @knn: the k nearest neighbor for every point in data [matrix]
.safelevelsmote_apply = function(index, data, knn){  
  #one point of minority
  index_p = index
  p = data[index_p,]
  #remove label
  p = p[-ncol(data)]  
  
  # get n: n one of p's nearest neighbors 
  index_n = sample(knn[index_p,], 1, replace = T)
  n = data[index_n,]
  n = n[-ncol(data)]  
  
  # calculate safelevel(p)
  knn_matrix = data[knn[index_p,],]
  safelevel_p = sum(knn_matrix[,ncol(knn_matrix)])
  # calculate safelevel(n)
  knn_matrix = data[knn[index_n,],]
  safelevel_n = sum(knn_matrix[,ncol(knn_matrix)])
  
  # calculate safelevelratio
  safelevelratio = .calculateSafelevelratio(safelevel_p, safelevel_n)  
  
  # calculate synthetic instance
  if(safelevelratio == Inf && safelevel_p == 0 ){
    # No synthetic point will be generated
    return (NULL)
  }
  else{
    # Genearte new synthetic point
    return (.populate_single_sl(safelevelratio, safelevel_p, safelevel_n, p, n))
  }
}

# Creates based on the given points and safelevels a new point
# @safelevelratio: the safelevelratio of safelevel_p and safelevel_n [double]
# @safelevel_p: safelevel of p [int]
# @safelevel_n: safelevel of n [int]
# @p: point p from the minority class [vector]
# @n: point n, one of ps k nearest neighbors [vector]
# @k: number of nearest neighbors to consider [int]
.populate_single_sl = function(safelevelratio, safelevel_p, safelevel_n, p, n){
  # Calculate difference betwenn p and n  
  dif = n - p
  
  # calculate gap
  # second case
  if(safelevelratio == Inf && safelevel_p != 0 ){
    gap = 0;
  }
  # third case
  else if(safelevelratio == 1){
    # Random betwenn 0 and 1
    gap = runif(1, 0, 1)
  } 
  # fourth case
  else if(safelevelratio > 1){
    # Random betwenn 0 and 1/safelevelratio
    gap = runif(1, 0, 1/safelevelratio)
  }
  # fifth case
  else if(safelevelratio <1){
    # Random betwenn 1-safelevelratio and 1
    gap = runif(1, 1-safelevelratio, 1)
  }
  # Calculate synthetic instance  
  return (p + gap * dif)
}

# Calculates the safelevelratio for two given safelevels
# @safelevel_p: safelevel of p [int]
# @safelevel_n: safelevel of n [int]
.calculateSafelevelratio = function(safelevel_p, safelevel_n){
  # cant divide by 0
  if(safelevel_n == 0){
    return (Inf)
  }
  else{
    return (safelevel_p/safelevel_n)
  }
}

# Checks the given input on validity
.checkinput_sl = function(data.nrow, k, seed){
  if(k <= 0 || length(k)!= 1 || k%%1!=0){
    stop("k has to be a positive integer")
  }  
  if(class(seed) != "numeric" && class(seed) != "integer"){
    stop("seed has to be a number")
  }
  if(k > data.nrow){
    stop("Not enough points to find k neigbhors")
  }
}

# Checks the given input on validity
.checkinput_sl_sep = function(data, minority){
  if(ncol(data) != ncol(minority)){
    stop("data and minority differ in dimension")
  }  
  if(!is.data.frame(data) && !is.matrix(data)){
    stop("data has to be a data.frame or matrix")
  }
  if(!is.data.frame(minority) && !is.matrix(minority)){
    stop("minority has to be a data.frame or matrix")
  }
  if(class(data) != class (minority)){
    stop("data and minority must have the same class")
  }
}