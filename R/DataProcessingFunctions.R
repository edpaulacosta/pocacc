

#' @export
generateAccumulationValuesFromListBMDValues = function(bmd){
  
  bmd = sort(as.numeric(as.character(bmd)))
  
  tm = table(bmd)
  
  if(length(tm)>1){
    
    length(bmd)
    accumulated = rep(0,length(tm))
    accumulated[1]=tm[1]
    for(j in 2:length(tm)){
      accumulated[j]=accumulated[j-1]+tm[j]
    }	
  }else{
    accumulated = rep(0,length(tm))
    accumulated[1]=tm[1]
    
  }	
  
  bmd_Values_col = as.numeric(names(tm))
  
  #plot(bmdValues,accumulated,log="x")
  
  return(cbind(bmd_Values_col,accumulated))
  
  
}
