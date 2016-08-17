        FUN = function(InputMat, log=TRUE){
          if(log==TRUE) return(c( exp(InputMat[,1]-1.96*InputMat[,2]), rev(exp(InputMat[,1]+1.96*InputMat[,2]))))
          if(log==FALSE) return(c( InputMat[,1]-1.96*InputMat[,2], rev(InputMat[,1]+1.96*InputMat[,2])))
        } 