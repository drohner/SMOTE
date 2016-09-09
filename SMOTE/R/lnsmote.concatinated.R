lnsmote.concatinated <-
function(data, col, char, k = 5, o = 1, method = "e", seed = 321, cores = 2){
  #Check Input
  #Check col and char
  .checkinput_ln_con(data, col, char)
  .checkinput_ln(nrow(data),k, o, seed)
  
  #Save the col names for reconstruction
  data.colnames = colnames(data)
  
  #Label every minority classmember in the last column with 1, other points with 0  
  #get minority
  minority = subset(data, data[,col] == char) 
  minority = minority[,-col]
  
  #get data without minority
  data = subset(data, data[,col] != char)
  data = data[,-col]
  
  #add labels
  minority = cbind(minority, 1)
  data = cbind(data, 0)  
  
  #concatinate data, so the indices of minority are 1-nrow(minority)
  names(minority) = names(data)
  data = rbind(minority, data)    
  
  #calculate synthetic points
  smotepoints = .lnsmote(data, nrow(minority), k,  o, method, seed, cores)  

  #no points were created
  if(is.null(smotepoints)){
    return (NULL)
  }

  #Rebuild the given structure
  if (col == 1){
    synthetic = data.frame(char, smotepoints)
  }
  else if(col == (ncol(smotepoints)+1)){
    synthetic = data.frame(smotepoints, char)
  }
  else{
    synthetic = data.frame(smotepoints[,1:col-1], char, smotepoints[,col:ncol(smotepoints)])  
  }
 
  colnames(synthetic) = data.colnames
  return (synthetic)  
}

# Function called by the facade methods
# Generates the synthetic instances by calling a lapply function
# @data: the concatinated datatset [matrix]
# @minority.rows: the number of rows of minority [matrix]
# @k: number of nearest neighbors to consider [int]
# @o: oversamplingratio [int]
# @method: the datastructure to use for knn search [int]
# @seed: set the random generator [int]
# @cores: the number of cores to use for parallelization [int]
.lnsmote = function(data, minority.rows, k, o, method,  seed, cores){  
  #get knn  
  method = .transformMethod(method, ncol(data))
  
  # k +1 because the search will return the point itself as the nearest neighbor. removed after the seach is finished
  # k +1 +1 for modification of gap
  knn = knnsearch(as.matrix(data[,-ncol(data)]) ,as.matrix(data[,-ncol(data)]), k+1+1, method, cores)[,-1]  
  
  set.seed(seed)  
  data = as.matrix(data)
  
  #use oversamplingratio
  indices = rep(c(1:minority.rows), o)
  
  #call apply function over the indcies. because the minority class samples are the saved in 1:minority.rows, we can use this approach
  smotepoints = do.call(rbind, lapply(indices, FUN = ".lnsmote_apply", data = data, knn = knn))
  
  smotepoints.nrow = 0
  if(!is.null(smotepoints)){
    smotepoints.nrow = nrow(smotepoints)
  }
  
  cat("The minority class got", minority.rows ,"samples.\n")
  cat("The oversamplingratio was",o,".\n")
  cat("Theoretically",minority.rows*o,"new samples could have been generated.\n")
  cat("Due to noisy samples, only", smotepoints.nrow,"synthetic instances were generated.\n")
  return (smotepoints)
}

# Calculates the synthetic instances based on safe level smote
# This funtion ist thought to be used in an lapply function
# @index: the index of the minority class sample in data [int]
# @data: the complete dataset [matrix]
# @knn: the k nearest neighbor for every point in data [matrix]
.lnsmote_apply = function(index, data, knn){
  #only consider the k nearest neighbors (the last one may be used for gap modification)
  rk = ncol(knn)-1
  
  #one point of minority without label
  index_p = index
  p = data[index_p,]
  p = p[-length(p)]
  
  # get n: n one of p's nearest neighbors without label
  index_n = sample(knn[index_p,1:rk], 1, replace = T)
  n = data[index_n,] 
  n = n[-length(n)]
  
  # calculate safelevel(p)
  knn_matrix = data[knn[index_p,][1:rk],]
  safelevel_p = sum(knn_matrix[,ncol(knn_matrix)])
  
  
  # calculate safelevel(n)
  knn_index = knn[index_n,][1:rk]
  if(data[index_n, ncol(data)] != 1 && index_p %in% knn_index){
    knn[index_n, match(index_p, knn_index)] = knn[index_n, rk+1]
  }
  safelevel_n = sum(knn_matrix[,ncol(knn_matrix)])
  
  
  if(safelevel_p == 0 && safelevel_n == 0 ){
    # No synthetic instance will be generated    
    return (NULL)
  }
  else{
    # Genearte new synthetic point
    return (.populate_single_ln(safelevel_p, safelevel_n, p, n, ncol(knn)-1))
  }
}

# Creates based on the given points and safelevels a new point
# @safelevel_p: safelevel of p [int]
# @safelevel_n: safelevel of n [int]
# @p: point p from the minority class [vector]
# @n: point n, one of ps k nearest neighbors [vector]
# @k: number of nearest neighbors to consider [int]
.populate_single_ln = function(safelevel_p, safelevel_n, p, n, k){
  #calculate gap
  if(safelevel_n == 0 && safelevel_p > 0){
    gap = 0
  }
  else{
    safelevelratio = safelevel_p/safelevel_n
    if(safelevelratio == 1){
      gap = runif(1, 0, 1)
    }
    else if(safelevelratio > 1){
      gap = runif(1, 0, 1/safelevelratio)
    }
    else{
      gap = 1-runif(1, 0, safelevelratio)
    }
  }
  if(n[length(n)] == 0){
    gap = gap * (safelevel_n/k)
  }
  
  dif = n-p
  # Calculate synthetic instance    
  return (p + gap * dif)
}

# Checks the given input on validity
.checkinput_ln = function(nrow.data, k, o, seed){
  if(k <= 0 || length(k)!= 1 || k%%1!=0){
    stop("k has to be a positive integer")
  }  
  if(o <= 0 || length(o)!= 1 || o%%1!=0){
    stop("o has to be a positive integer")
  }  
  if(class(seed) != "numeric" && class(seed) != "integer"){
    stop("seed has to be a number")
  }
  if(k > nrow.data){
    stop("Not enough points to find k neigbhors")
  }
}

# Checks the given input on validity
.checkinput_ln_con = function(data, col, char){
  if(!(char %in% data[,col])){
    stop("character cannot be found in the given column")
  }  
  if(col > ncol(data)){
    stop("col has to be less or equal to the number of columns of data")
  }
  if(!is.data.frame(data) && !is.matrix(data)){
    stop("data has to be a data.frame or matrix")
  }
  if(!is.character(char)){
    stop("char has to be a character")
  }
}