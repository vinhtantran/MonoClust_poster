cv.plot <- function(cptable, ...) {
    plot(1:9, cp.table[,1], type="l", lty=2, 
         ylim=c(0, max(cp.table[,1] + cp.table[,2])), 
         main = "MSE for CV of monothetic clustering",
         ylab = "MSE +/- 1 SE",
         xlab = "Number of clusters",
         ...)
    
    error.bar(1:9, cp.table[,1], cp.table[,2], col="red")
    
}

# From Monkey's Uncle blog at http://monkeysuncle.stanford.edu/?p=485
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}