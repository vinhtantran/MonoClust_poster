
PCAmix <- function (X.quanti = NULL, X.quali = NULL, ndim = 5, graph = TRUE) 
{
    cl <- match.call()
    rec <- recod(X.quanti, X.quali)
    n <- rec$n
    p <- rec$p
    p1 <- rec$p1
    X <- rec$X
    Z <- rec$Z
    G <- rec$G
    indexj <- rec$indexj
    Ztilde <- Z/sqrt(n)
    e <- svd(Ztilde)
    Tot <- sum(e$d^2)
    if (ndim <= 1) 
        stop("'ndim' must be an interger greater or equal to 2")
    if (ndim > length(e$d)) 
        stop("'ndim' must be an integer smaller than min(nrow, p1+p2-m-1)")
    U <- e$u[, 1:ndim] * sqrt(n)
    rownames(U) <- rownames(Z)
    d <- e$d[1:ndim]
    F <- U %*% diag(d)
    A <- e$v[, 1:ndim] %*% diag(d)
    rownames(A) <- colnames(Z)
    A1 <- NULL
    A2 <- NULL
    A2coord <- NULL
    C <- NULL
    ps <- NULL
    if (p1 > 0) {
        A1 <- A[1:p1, ]
        colnames(A1) <- paste("dim", 1:ndim, sep = "")
    }
    if (p1 != p) {
        if (p1 > 0) 
            A2 <- A[-(1:p1), ]
        else A2 <- A
        ns <- apply(G, 2, sum)
        ps <- ns/n
        A2coord <- sweep(A2, MARGIN = 1, STATS = sqrt(ps), FUN = "/")
        colnames(A2coord) <- paste("dim", 1:ndim, sep = "")
        C <- matrix(NA, (p - p1), length(d))
        rownames(C) <- colnames(X.quali)
        colnames(C) <- paste("dim", 1:ndim, sep = "")
        for (j in 1:(p - p1)) {
            C[j, ] <- apply(A2[which(indexj == (j + p1)) - p1, 
                ], 2, FUN = function(x) {
                sum(x^2)
            })
        }
    }
    sload <- rbind(A1^2, C)
    names(d) <- colnames(U) <- colnames(F) <- colnames(A) <- colnames(sload) <- paste("dim", 
        1:ndim, sep = "")
    res <- list(call = cl, eig = d^2, Tot = Tot, scores.stand = U, 
        scores = F, sload = sload, A = A, categ.coord = A2coord, 
        quanti.cor = A1, quali.eta2 = C, rec = rec, ndim = ndim)
    class(res) <- "PCAmix"
    if (graph) {
        plot.PCAmix(res, main = "Scores")
        if (p1 != p) {
            dev.new()
            plot.PCAmix(res, choice = "categ", main = "Categories")
        }
        if (p1 != 0) {
            dev.new()
            plot.PCAmix(res, choice = "cor", main = "Correlation circle")
        }
        dev.new()
        plot.PCAmix(res, choice = "var", main = "Squared loadings")
    }
    return(res)
}


PCArot <- function (obj, dim, itermax = 100, graph = TRUE) 
{
    cl <- match.call()
    if (!inherits(obj, "PCAmix")) 
        stop("use only with \"PCAmix\" objects")
    if (dim <= 1) 
        stop("'dim' must be an interger greater or equal to 2")
    if (dim > length(obj$eig)) 
        stop("'dim' must be an integer smaller or equal to 'ndim' specified in 'PCAmix'")
    A <- obj$A[, 1:dim]
    scores.stand <- obj$scores.stand[, 1:dim]
    indexj <- obj$rec$indexj
    p1 <- obj$rec$p1
    p <- obj$rec$p
    Z <- obj$rec$Z
    G <- obj$rec$G
    X.quali <- obj$rec$X.quali
    n <- obj$rec$n
    Tot <- obj$Tot
    if (dim == 2) {
        res <- sol.2dim(A, indexj, p, p1)
        theta <- res$theta
        iter <- 1
        T <- res$T
        A.rot <- A %*% T
        scores.stand.rot <- scores.stand %*% T
    }
    else {
        matcombn <- combn(1:dim, 2)
        nbpaires <- ncol(matcombn)
        A.rot <- A
        scores.stand.rot <- scores.stand
        cptthetazero <- 0
        iter <- 0
        while (cptthetazero < nbpaires) {
            cptthetazero <- 0
            iter <- iter + 1
            for (j in 1:nbpaires) {
                indicescol <- matcombn[, j]
                res <- sol.2dim(A.rot[, indicescol], indexj, 
                  p, p1)
                theta <- res$theta
                T <- res$T
                A.rot[, indicescol] <- A.rot[, indicescol] %*% 
                  T
                scores.stand.rot[, indicescol] <- scores.stand.rot[, 
                  indicescol] %*% T
                if (round(theta, digit = 3) == 0) 
                  cptthetazero <- cptthetazero + 1
            }
            if (iter > itermax) 
                stop("Stop: maximum number of iterations reached.")
        }
        theta <- NULL
        T <- NULL
    }
    A1 <- NULL
    A2 <- NULL
    A2coord <- NULL
    C <- NULL
    ps <- NULL
    if (p1 > 0) {
        A1 <- A.rot[1:p1, ]
        colnames(A1) <- paste("dim", 1:dim, sep = "", ".rot")
    }
    if (p1 != p) {
        if (p1 > 0) 
            A2 <- A.rot[-(1:p1), ]
        else A2 <- A.rot
        ns <- apply(G, 2, sum)
        ps <- ns/n
        A2coord <- sweep(A2, MARGIN = 1, STATS = sqrt(ps), FUN = "/")
        colnames(A2coord) <- paste("dim", 1:dim, sep = "", ".rot")
        C <- matrix(NA, (p - p1), dim)
        rownames(C) <- colnames(X.quali)
        colnames(C) <- paste("dim", 1:dim, sep = "", ".rot")
        for (j in 1:(p - p1)) {
            C[j, ] <- apply(A2[which(indexj == (j + p1)) - p1, 
                ], 2, FUN = function(x) {
                sum(x^2)
            })
        }
    }
    sload.rot <- rbind(A1^2, C)
    var.rot <- apply(sload.rot, 2, sum)
    scores.rot <- sweep(scores.stand.rot, 2, STATS = sqrt(var.rot), 
        FUN = "*")
    colnames(sload.rot) <- names(var.rot) <- colnames(scores.rot) <- colnames(scores.stand.rot) <- paste("dim", 
        1:dim, sep = "", ".rot")
    res <- list(call = cl, theta = theta, T = T, eig = var.rot, 
        Tot = Tot, sload = sload.rot, scores.stand = scores.stand.rot, 
        scores = scores.rot, categ.coord = A2coord, quanti.cor = A1, 
        quali.eta2 = C, ndim = dim, rec = obj$rec, iter = iter)
    class(res) <- "PCAmix"
    if (graph) {
        plot.PCAmix(res, main = "Scores after rotation")
        if (p1 != p) {
            dev.new()
            plot.PCAmix(res, choice = "categ", main = "Categories after rotation")
        }
        if (p1 != 0) {
            dev.new()
            plot.PCAmix(res, choice = "cor", main = "Correlation circle after rotation")
        }
        dev.new()
        plot.PCAmix(res, choice = "var", main = "Squared loadings after rotation")
    }
    return(res)
}

summary.PCAmix <- function (object, ...) 
{
    x <- object
    if (!inherits(x, "PCAmix")) 
        stop("use only with \"PCAmix\" objects")
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    n <- x$rec$n
    p1 <- x$rec$p1
    p <- x$rec$p
    p2 <- p - p1
    if (colnames(x$sload)[1] == "dim1.rot") 
        cat("Method = rotation after ")
    else cat("Method = ")
    if (p1 == p) 
        cat("Principal Component Analysis (PCA)")
    if (p1 == 0) 
        cat("Multiple Correspondence Analysis (MCA)")
    if ((p1 != 0) && (p1 != p)) 
        cat("Factor Analysis of mixed data (FAmix)")
    cat("\n")
    cat("\n")
    cat("Data:", "\n")
    cat(paste("   number of observations: ", n), sep = " ")
    cat("\n")
    if ((p1 != 0) && (p2 == 0)) {
        cat(paste("   number of variables: ", p1), sep = " ")
        cat("\n")
    }
    if ((p1 == 0) && (p2 != 0)) {
        cat(paste("   number of variables: ", p2), sep = " ")
        cat("\n")
    }
    if ((p1 != 0) && (p2 != 0)) {
        cat(paste("   number of  variables: ", p), sep = " ")
        cat("\n")
        cat(paste("        number of numerical variables: ", 
            p1), sep = " ")
        cat("\n")
        cat(paste("        number of categorical variables: ", 
            p2), sep = " ")
        cat("\n")
    }
    cat("\n")
    if (colnames(x$sload)[1] == "dim1.rot") 
        cat("Squared loadings after rotation:")
    else cat("Squared loadings :")
    cat("\n")
    print(round(x$sload, digit = 2))
    cat("\n")
    cat("\n")
}


plot.PCAmix <- function (x, axes = c(1, 2), choice = "ind", stand = FALSE, label = TRUE, 
    quali = NULL, posleg = "topleft", xlim = NULL, ylim = NULL, 
    ...) 
{
    if (!inherits(x, "PCAmix")) 
        stop("use only with \"PCAmix\" objects")
    if (max(axes) > x$ndim) 
        stop(paste("axes must be between 1 and ", x$ndim, sep = ""))
    if (!(choice %in% c("ind", "var", "categ", "cor"))) 
        stop("the argument 'choice' must be 'ind','var','cor' or 'categ'")
    dim1 <- axes[1]
    dim2 <- axes[2]
    if (stand) {
        lab.x <- paste("Dimension ", dim1, sep = "")
        lab.y <- paste("Dimension ", dim2, sep = "")
        scores <- x$scores.stand
    }
    else {
        px <- round(x$eig[dim1]/x$Tot * 100, digit = 2)
        py <- round(x$eig[dim2]/x$Tot * 100, digit = 2)
        lab.x <- paste("Dim ", dim1, " (", px, "%)", sep = "")
        lab.y <- paste("Dim ", dim2, " (", py, "%)", sep = "")
        scores <- x$scores
    }
    p1 <- x$rec$p1
    p <- x$rec$p
    quanti.coord <- x$quanti.cor
    if (choice == "ind") {
        if (is.null(xlim)) {
            xmin <- min(scores[, dim1])
            xmax <- max(scores[, dim1])
            xlim <- c(xmin, xmax) * 1.2
        }
        if (is.null(ylim)) {
            ymin <- min(scores[, dim2])
            ymax <- max(scores[, dim2])
            ylim <- c(ymin, ymax) * 1.2
        }
        if (is.null(quali)) {
            plot(scores[, axes], xlim = xlim, ylim = ylim, xlab = lab.x, 
                ylab = lab.y, pch = 20, ...)
            abline(h = 0, lty = 2)
            abline(v = 0, lty = 2)
            if (label) 
                text(scores[, axes], labels = rownames(scores), 
                  pos = 3, ...)
        }
        else {
            if (is.numeric(quali)) 
                stop("quali must be categorical")
            if (length(quali) != nrow(x$scores)) 
                stop("the length of quali is inapproriate")
            quali <- as.factor(quali)
            plot(scores[, axes], xlim = xlim, ylim = ylim, xlab = lab.x, 
                ylab = lab.y, pch = 20, col = as.numeric(quali), 
                ...)
            abline(h = 0, lty = 2)
            abline(v = 0, lty = 2)
            legend(posleg, legend = levels(quali), text.col = c(1:length(levels(quali))), 
                cex = 0.8)
            if (label) 
                text(scores[, axes], labels = rownames(x$scores), 
                  pos = 3, col = as.numeric(quali), ...)
        }
    }
    if (choice == "var") {
        if (is.null(xlim)) {
            xmax <- max(x$sload[, dim1])
            xlim <- c(-0.1, xmax * 1.2)
        }
        if (is.null(ylim)) {
            ymax <- max(x$sload[, dim2])
            ylim <- c(-0.1, ymax * 1.2)
        }
        plot(0, 0, type = "n", xlab = lab.x, ylab = lab.y, xlim = xlim, 
            ylim = ylim, ...)
        abline(v = 0, lty = 2)
        abline(h = 0, lty = 2)
        for (j in 1:nrow(x$sload)) {
            arrows(0, 0, x$sload[j, dim1], x$sload[j, dim2], 
                length = 0.1, angle = 15, code = 2)
            if (label) {
                if (x$sload[j, dim1] > x$sload[j, dim2]) {
                  pos <- 4
                }
                else pos <- 3
                text(x$sload[j, dim1], x$sload[j, dim2], labels = rownames(x$sload)[j], 
                  pos = pos)
            }
        }
    }
    if (choice == "categ") {
        if (p1 == p) 
            stop("Argument 'categ' is inappropriate for only quantitative data")
        if (is.null(xlim)) {
            xmin <- min(x$categ.coord[, dim1])
            xmax <- max(x$categ.coord[, dim1])
            xlim <- c(xmin, xmax) * 1.2
        }
        if (is.null(ylim)) {
            ymin <- min(x$categ.coord[, dim2])
            ymax <- max(x$categ.coord[, dim2])
            ylim <- c(ymin, ymax) * 1.2
        }
        plot(x$categ.coord[, axes], xlim = xlim, ylim = ylim, 
            xlab = lab.x, ylab = lab.y, pch = 17, ...)
        abline(h = 0, lty = 2)
        abline(v = 0, lty = 2)
        if (label) 
            text(x$categ.coord[, axes], labels = rownames(x$categ.coord), 
                pos = 3, ...)
    }
    if (choice == "cor") {
        if (p1 == 0) 
            stop("Argument 'cor' is inappropriate for only qualitative data")
        if (is.null(xlim)) 
            xlim <- c(-1, 1) * 1.3
        if (is.null(ylim)) 
            ylim <- c(-1, 1) * 1.3
        plot(0, 0, type = "n", xlab = lab.x, ylab = lab.y, xlim = xlim, 
            ylim = ylim, ...)
        x.cercle <- seq(-1, 1, by = 0.01)
        y.cercle <- sqrt(1 - x.cercle^2)
        lines(x.cercle, y = y.cercle)
        lines(x.cercle, y = -y.cercle)
        abline(v = 0, lty = 2)
        abline(h = 0, lty = 2)
        for (j in 1:nrow(quanti.coord)) {
            arrows(0, 0, quanti.coord[j, dim1], quanti.coord[j, 
                dim2], length = 0.1, angle = 15, code = 2)
            if (label) {
                if (abs(quanti.coord[j, dim1]) > abs(quanti.coord[j, 
                  dim2])) {
                  if (quanti.coord[j, dim1] >= 0) 
                    pos <- 4
                  else pos <- 2
                }
                else {
                  if (quanti.coord[j, dim2] >= 0) 
                    pos <- 3
                  else pos <- 1
                }
                text(quanti.coord[j, dim1], quanti.coord[j, dim2], 
                  labels = rownames(quanti.coord)[j], pos = pos, 
                  ...)
            }
        }
    }
}


recod <- function (X.quanti, X.quali) 
{
    G <- NULL
    Gcod <- NULL
    if (!is.null(X.quanti)) {
        if (is.factor(X.quanti)) 
            stop("All variables in X.quanti must be numerical")
        if (is.numeric(X.quanti)) 
            X.quanti <- data.frame(X.quanti)
        for (v in 1:ncol(X.quanti)) {
            if (!is.numeric(X.quanti[, v])) 
                stop("All variables in X.quanti must be numeric")
        }
        n1 <- nrow(X.quanti)
        p1 <- ncol(X.quanti)
        Z1 <- recodquant(X.quanti)
        Xn1 <- unique(as.matrix(X.quanti), MARGIN = 2)
    }
    if (!is.null(X.quali)) {
        if (is.numeric(X.quali)) 
            stop("All variables in X.quali must be categorical")
        if (is.factor(X.quali)) 
            X.quali <- data.frame(X.quali)
        for (v in 1:ncol(X.quali)) {
            if (is.numeric(X.quali[, v])) 
                stop("All variables in X.quali must be categorical")
        }
        n2 <- nrow(X.quali)
        p2 <- ncol(X.quali)
        G <- recodqual(X.quali)
        ns <- apply(G, 2, sum)
        ps <- ns/nrow(G)
        Gcod <- sweep(G, MARGIN = 2, STATS = sqrt(ps), FUN = "/")
        moy <- apply(Gcod, 2, mean)
        Z2 <- sweep(Gcod, MARGIN = 2, STATS = moy, FUN = "-")
        Xn2 <- unique(as.matrix(X.quali), MARGIN = 2)
        nb.moda <- function(moda) {
            moda <- as.factor(moda)
            length(levels(moda))
        }
        nbmoda <- apply(X.quali, 2, nb.moda)
        indexj2 <- NULL
        for (j in 1:ncol(X.quali)) {
            indexj2 <- c(indexj2, rep(j, nbmoda[j]))
        }
    }
    if (!is.null(X.quanti) && !is.null(X.quali)) 
        if (n1 == n2) {
            n <- n1
            p <- p1 + p2
            Z <- cbind(Z1, Z2)
            X <- cbind.data.frame(X.quanti, X.quali)
            Xn <- cbind.data.frame(Xn1, Xn2)
            indexj <- c(1:p1, indexj2 + p1)
        }
        else stop("number of objects in X.quanti and X.quali must be the same")
    if (!is.null(X.quanti) && is.null(X.quali)) {
        n <- n1
        p <- p1
        Z <- Z1
        X <- X.quanti
        Xn <- Xn1
        indexj <- 1:p1
    }
    if (is.null(X.quanti) && !is.null(X.quali)) {
        n <- n2
        p <- p2
        p1 <- 0
        Z <- Z2
        X <- X.quali
        Xn <- Xn2
        indexj <- indexj2
    }
    if (is.null(X.quanti) && is.null(X.quali)) 
        stop("A data matrix must be given")
    if (is.null(colnames(X))) 
        colnames(X) <- paste("V", 1:ncol(X), sep = "")
    for (j in 1:ncol(X)) if (colnames(X)[j] == "") 
        colnames(X)[j] <- paste("V", j, sep = "")
    return(list(X = X, Xn = Xn, Z = Z, n = n, p = p, p1 = p1, 
        indexj = indexj, G = G, Gcod = Gcod, X.quanti = X.quanti, 
        X.quali = X.quali))
}


recodquant <- function (X) 
{
    X <- as.matrix(X)
    Xcod <- apply(X, 2, missing.mean)
    red <- sqrt((nrow(X) - 1)/nrow(X))
    Z <- scale(Xcod, scale = (apply(Xcod,2,sd) * red))
    apply(Z, 1, function(x) sum(is.na(x)))
    if (sum(is.na(Z)) != 0) 
        stop("There are columns in X.quanti where all the values are identical")
    return(Z)
}


recodqual <- function (X) 
{
    X <- as.matrix(X)
    GNA <- tab.disjonctif.NA(X)
    G <- replace(GNA, is.na(GNA), 0)
    ns <- apply(G, 2, sum)
    nmiss <- apply((apply(GNA, 2, is.na)), 2, sum)
    n <- nrow(X)
    if (sum((n - nmiss) == ns) != 0) 
        stop("There are columns in X.quali where all the categories are identical")
    return(G)
}

missing.mean <- function (C1) 
{
    moy <- mean(C1, na.rm = T)
    ind <- which(is.na(C1) == T)
    if (length(ind) >= 1) {
        C1[ind] <- moy
    }
    return(C1)
}

tab.disjonctif.NA <- function (tab) 
{
    tab <- as.data.frame(tab)
    modalite.disjonctif <- function(i) {
        moda <- tab[, i]
        nom <- names(tab)[i]
        n <- length(moda)
        moda <- as.factor(moda)
        x <- matrix(0, n, length(levels(moda)))
        ind <- (1:n) + n * (unclass(moda) - 1)
        indNA <- which(is.na(ind))
        x[(1:n) + n * (unclass(moda) - 1)] <- 1
        x[indNA, ] <- NA
        if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda), 
            "n", "N", "y", "Y"))) 
            dimnames(x) <- list(row.names(tab), paste(nom, levels(moda), 
                sep = "."))
        else dimnames(x) <- list(row.names(tab), levels(moda))
        return(x)
    }
    if (ncol(tab) == 1) 
        res <- modalite.disjonctif(1)
    else {
        res <- lapply(1:ncol(tab), modalite.disjonctif)
        res <- as.matrix(data.frame(res, check.names = FALSE))
    }
    return(res)
}

