perm.plot <- function(object, mds.data, node = 1, colors){
    
    if (!inherits(object, "MonoClust")) stop("Not a legitimate \"MonoClust\" object")
    
    if (length(colors) < 3) colors <- c(1, 2, 5)
    
    frame <- object$frame
    nodes <- as.numeric(row.names(frame))
    # TODO: Check if node is <leaf> --> error
    if (node > max(as.numeric(row.names(frame))) |
            !(node %in% nodes)) stop("The chosen node is too large for the current tree")
    
    children <- find.children(nodes, node)
    members <- object$Membership
    newmember <- members
    newmember[members %in% children$left] <- colors[1]
    newmember[members %in% children$right] <- colors[2]
    newmember[!(members %in% c(children$left, children$right))] <- colors[3]
    
    plot(mds.data, col=newmember)
}

find.children <- function(node.list, node) {
    left.node <- node*2
    right.node <- node*2 + 1
    max.power <- floor(log(max(node.list), base=2))
    start.power <- floor(log(node, base=2)) + 2
    
    if (start.power > max.power) {
        l <- left.node
        r <- right.node
    } else {
        l <- left.node
        r <- right.node
        for (i in 1:(max.power - start.power + 1)) {
            l <- c(l, seq(left.node * 2^i , left.node * 2^i + 2^i - 1))
            r <- c(r, seq(right.node * 2^i , right.node * 2^i + 2^i - 1))
        }
    }
    return(list(left = l, right=r))
}