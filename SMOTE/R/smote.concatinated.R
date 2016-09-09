smote.concatinated <-
function(data, col, char, k = 5, o = 1, method = "e", seed = 321, cores = 2){
  #Check Input
  #Check col and char
  .checkinput_con(data, col, char)
  .checkinput(nrow(data),k, o, seed)
  
  #Save the col names for reconstruction
  data.colnames = colnames(data)
  
  #get minority
  minority = subset(data, data[,col] == char) 
  minority = minority[,-col]
    
  #calculate synthetic points
  smotepoints = .smote(as.matrix(minority), k, o, method, seed, cores)  
    
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
# @minority: the instances of the minority class [matrix]
# @k: number of nearest neighbors to consider [int]
# @o: oversamplingratio [int]
# @method: the datastructure to use for knn search [int]
# @seed: set the random generator [int]
# @cores: the number of cores to use for parallelization [int]
.smote = function(minority, k, o, method, seed , cores){    
  set.seed(seed)
  
  #if o < 1, randomize the data for oversampling
  if(o < 1){
    randomindices = sample(c(1:nrow(minority)), nrow(minority)*o)
    minority = minority[randomindices,]
    indices = 1:nrow(minority)
  }
  else{
    indices = rep(c(1:nrow(minority)), o)
  }  
  
  #parallelization doesnt work for windows
  if(Sys.info()['sysname'] == 'Windows'){
    cores = 1
  }
  
  #get knn  
  method = .transformMethod(method, ncol(minority))  
  # k +1 because the search will return the point itself as the nearest neighbor. removed after the seach is finished
  knn = knnsearch(as.matrix(minority[,-ncol(minority)]) ,as.matrix(minority[,-ncol(minority)]), k+1, method, cores)[,-1]  
  
  synthetic= do.call(rbind, lapply(indices, FUN = ".smote_apply", minority = minority, knn = knn))
  
  cat("The minority class got", nrow(minority) ,"samples.\n")
  cat("The oversamplingratio was",o,".\n")
  cat(nrow(synthetic),"synthetic instances were generated.\n")
  
  return (synthetic)
}

# Calculates the synthetic instances based on safe level smote
# This funtion ist thought to be used in an lapply function
# @index: the index of the minority class sample in data [int]
# @data: the complete dataset [matrix]
# @knn: the k nearest neighbor for every point in data [matrix]
.smote_apply = function(index, minority, knn){
  #one point of minority
  index_p = index
  p = minority[index_p,]
  
  #get n: n one of p's nearest neighbors 
  index_n = sample(knn[index_p,], 1)
  n = minority[index_n,]
  
  #genrate random gap between 0 and 1
  gap = runif(1, 0, 1)
  dif = n - p
  return (p + gap * dif)
}

# Checks the given input on validity
.checkinput = function(data.nrow, k, o, seed){
  if(k <= 0 || length(k)!= 1 || k%%1!=0){
    stop("k has to be a positive integer")
  }  
  if(o <= 0 || length(o)!= 1 || o%%1!=0){
    stop("o has to be a positive integer")
  }  
  if(class(seed) != "numeric" && class(seed) != "integer"){
    stop("seed has to be a number")
  }
  if(k > data.nrow){
    stop("Not enough points to find k neigbhors")
  }
}

# Checks the given input on validity
.checkinput_con = function(data, col, char){
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
