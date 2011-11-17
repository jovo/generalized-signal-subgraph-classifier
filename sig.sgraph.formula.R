sig.sgraph.formula <- function(formula, data = NULL, ..., subset, na.action = na.omit, scale = TRUE) 
{
	call <- match.call()
    # if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
        # "term.labels")) != 1L)) 
        # stop("'formula' missing or incorrect")
    # print(call)
    if (!inherits(formula, "formula")) 
        stop("method is only for formula objects")
    m <- match.call(expand.dots = FALSE)
    print(names(m))
    # names(m)[2] <- "formula"
    print(call)
    if (identical(class(eval.parent(m$data)), "matrix")) 
        m$data <- as.data.frame(eval.parent(m$data))
    # print(m$...)
    m$... <- NULL
    # m$scale <- NULL
    m[[1]] <- as.name("model.frame")
    # m$na.action <- na.action
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
#    attr(Terms, "intercept") <- 0
    ### making x into a data set
    x <- model.matrix(Terms, m)
	## getting y
    y <- model.extract(m, "response")
    attr(x, "na.action") <- attr(y, "na.action") <- attr(m, "na.action")
    
    #taken from svm
    if (length(scale) == 1) 
        scale <- rep(scale, ncol(x))
    if (any(scale)) {
        remove <- unique(c(which(labels(Terms) %in% names(attr(x, 
            "contrasts"))), which(!scale)))
        scale <- !attr(x, "assign") %in% remove
    }
    
    # print(x)
    print("In Formula")
    ret <- sig.sgraph.default(x=x, y=y, scale = scale, ..., na.action = na.action, formula=formula, data=data)
    ret$call <- call
    ret$call[[1]] <- as.name("sig.sgraph")
    ret$terms <- Terms
    if (!is.null(attr(m, "na.action"))) 
        ret$na.action <- attr(m, "na.action")
    class(ret) <- c("sig.sgraph.formula", class(ret))
    return(ret)
}



    # m <- model.frame(terms(reformulate(attributes(Terms)$term.labels)), 
        # data.frame(m))
    # for (i in seq(along = ncol(m))) {
        # if (is.ordered(m[[i]])) 
            # m[[i]] <- as.numeric(m[[i]])
    # }