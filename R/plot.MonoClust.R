plot.MonoClust<-function(x,margin,which,abbrev=4,text=TRUE,...){
    ## This function sets some defaults and changes things a bit, but is mostly a 
    ## wrapper for our slightly modified version of rpart's plot function (see plots.R).
    
    if(missing(margin)){margin<-c(.12,.02,0,.05)}
    if(missing(which)){which <- 4}
    
    plot.rpart(x,margin=margin,...)
    
    ## REMOVE: Tan, 3/1/15, Remove Inertia line
#     lines(x=c(.88,.88),y=c(0,1))
    
#     for(i in seq(0,1,.1)){
#         lines(x=c(.86,.88),y=c(i,i))
#         text(.73,i,i)
#     }
    
    if(text){
        text.MonoClust(x,which=which,abbrev=abbrev)
    }
    
    
}

Nclustplot<-function(x, main, type, ylab, xlab,...){
    
    if(missing(main)){main<-"Marginal Cluster Analysis"}
    if(missing(type)){type<-"b"}
    if(missing(ylab)){ylab<-"Proportion of Deviance Explained"}
    if(missing(xlab)){xlab<-"Number of Clusters"}
    
    inds <- seq(from=2,to=nrow(x$frame), by =2)
    plot(inds,round((1-as.numeric(x$frame$yval[inds])/1),digits=2),type=type, xaxt="n", ylab=ylab, xlab=xlab, main=main)
    axis(1, at=inds, labels= as.character(2 + 0:(length(inds)-1)))
    
}



