knnsearch <-
function(data, target = data, k, method = 'e', cores = 2){  
  method = .transformMethod(method, ncol(data))
  if(Sys.info()['sysname'] == 'Windows'){
    cores = 1
  }
  return(.Call("SMOTE_wrapper", as.matrix(data), as.matrix(target), k, method, cores, PACKAGE = "SMOTE")+1)
}

#transforms the char depicting the method to use for knn into an int
.transformMethod = function(char, data.cols){
  method = 0
  if(char == 'k'){
    return (1)
  }
  if(char == 'v'){
    return(2)
  }
  if(char == 'l'){
    return(3)
  }
  if(char == 'b'){
    return (4)
  }
  if(char == 'e'){
    if(data.cols <= 10){
      return (1)
    }
    else{
      return (4)
    }
  }
  if(char == 'a'){
    if(data.cols <= 10){
      return (1)
    }
    if(data.cols <= 500){
      return (3)
    }
    else{
      return (4)
    }
  }
  if(data.cols <= 10){
    return (1)
  }
  else{
    return (4)
  }
}