#judgeOrder() will be used to judge the result according to
#the given weight and coordinate characters,
#returning the final results of target vector of characters.
judgeOrder <- function( target = c("B","A"),
                        W = c("A","B","O"),wt =c(1,1,0)){
  Maxw <- 0
  reT <- vector()
  L <- length(target)
  for (i in 1:L){       #to avoid the invalid targets and W.
    if (length(which(W == target[i])) != 1) {
      return("Error!")
    }
    else{
      fw <- wt[which(W == target[i])]
      if(fw > Maxw){
        reT <- target[i]
        Maxw = fw
      }
      else if (fw == Maxw){
        # there is some paralleling elements like "A" & "B" in blood group
        if(length(which(reT == target[i])) == 0)reT <- c(reT,target[i])
        #add the new elements only
      }
      else {
        next
      }
    }
  }
  if (length(reT) > 1){
    #???ڲ???Ԫ??ʱ??????W?г??ֵ??Ⱥ?˳?????򣬼?AB??BA??ͬ
    i <- 1
    while(i < length(reT)){
      j <- length(reT)
      while (j > i){
        if(which(W == reT[j]) <  which(W == reT[j-1])){
          t <- reT[j]
          reT[j] <- reT[j-1]
          reT[j-1] <- t
        }
        j <- j-1
      }
      i <- i+1
    }
  }
  reT <- paste(reT,collapse = "")
  return(reT)
}
generate_samples_simulation <- function(size = 100,total = 2,
                               W=c("A","B","O"),wt=c(1,1,0),wp=c(0.8,0.1,0.1)){
  re.lastcol <- vector("character",size)
  re.newrow.f <- vector("character",total)
  for(i in 1:size){
    re.newrow.f <- sample(x=W,prob = wp,size = total,replace = TRUE)
    re.lastcol[i] <- judgeOrder(target = re.newrow.f,W = W,wt = wt)  }
  return(table(re.lastcol))
}


