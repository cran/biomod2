.somers2 = function (x, y, weights = NULL, normwt = FALSE, na.rm = TRUE)
{
    if (length(y) != length(x))
        stop("y must have same length as x")
    y <- as.integer(y)
    wtpres <- length(weights)
    if (wtpres && (wtpres != length(x)))
        stop("weights must have same length as x")
    if (na.rm) {
        miss <- if (wtpres)
            is.na(x + y + weights)
        else is.na(x + y)
        nmiss <- sum(miss)
        if (nmiss > 0) {
            miss <- !miss
            x <- x[miss]
            y <- y[miss]
            if (wtpres)
                weights <- weights[miss]
        }
    }
    else nmiss <- 0
    u <- sort(unique(y))
    if (any(!(y %in% 0:1) ))
        stop("y must be binary")
    if (wtpres) {
        if (normwt)
            weights <- length(x) * weights/sum(weights)
        n <- sum(weights)
    }
    else n <- length(x)
    if (n < 2)
        stop("must have >=2 non-missing observations")
    n1 <- if (wtpres)
        sum(weights[y == 1])
    else sum(y == 1)
    if (n1 == 0 || n1 == n)
        return(c(C = NA, Dxy = NA, n = n, Missing = nmiss))
    mean.rank <- if (wtpres){
      require(Hmisc,quietly=TRUE)
      wtd.mean(wtd.rank(x, weights, na.rm = FALSE), weights *
            y)
    } else mean(rank(x)[y == 1])
    c.index <- (mean.rank - (n1 + 1)/2)/(n - n1)
    dxy <- 2 * (c.index - 0.5)
    r <- c(c.index, dxy, n, nmiss)
    names(r) <- c("C", "Dxy", "n", "Missing")
    r
}

# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# # Hmisc internal functions to avoid Hmisc dependencies
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# .wtd.mean <- function (x, weights = NULL, normwt = "ignored", na.rm = TRUE) 
# {
#     if (!length(weights)) 
#         return(mean(x, na.rm = na.rm))
#     if (na.rm) {
#         s <- !is.na(x + weights)
#         x <- x[s]
#         weights <- weights[s]
#     }
#     sum(weights * x)/sum(weights)
# }
# 
# .wtd.rank <- function (x, weights = NULL, normwt = FALSE, na.rm = TRUE) 
# {
#     if (!length(weights)) 
#         return(rank(x), na.last = if (na.rm) NA else TRUE)
#     tab <- .wtd.table(x, weights, normwt = normwt, na.rm = na.rm)
#     freqs <- tab$sum.of.weights
#     r <- cumsum(freqs) - 0.5 * (freqs - 1)
#     approx(tab$x, r, xout = x)$y
# }
# 
# .wtd.table <- function (x, weights = NULL, type = c("list", "table"), normwt = FALSE, 
#     na.rm = TRUE) 
# {
#     type <- match.arg(type)
#     if (!length(weights)) 
#         weights <- rep(1, length(x))
#     isdate <- .testDateTime(x)
#     ax <- attributes(x)
#     ax$names <- NULL
#     x <- if (is.character(x)) 
#         .as.category(x)
#     else oldUnclass(x)
#     lev <- levels(x)
#     if (na.rm) {
#         s <- !is.na(x + weights)
#         x <- x[s, drop = FALSE]
#         weights <- weights[s]
#     }
#     n <- length(x)
#     if (normwt) 
#         weights <- weights * length(x)/sum(weights)
#     i <- order(x)
#     x <- x[i]
#     weights <- weights[i]
#     if (any(diff(x) == 0)) {
#         weights <- tapply(weights, x, sum)
#         if (length(lev)) {
#             levused <- lev[sort(unique(x))]
#             if ((length(weights) > length(levused)) && any(is.na(weights))) 
#                 weights <- weights[!is.na(weights)]
#             if (length(weights) != length(levused)) 
#                 stop("program logic error")
#             names(weights) <- levused
#         }
#         if (!length(names(weights))) 
#             stop("program logic error")
#         if (type == "table") 
#             return(weights)
#         x <- .all.is.numeric(names(weights), "vector")
#         if (isdate) 
#             attributes(x) <- c(attributes(x), ax)
#         names(weights) <- NULL
#         return(list(x = x, sum.of.weights = weights))
#     }
#     xx <- x
#     if (isdate) 
#         attributes(xx) <- c(attributes(xx), ax)
#     if (type == "list") 
#         list(x = if (length(lev)) lev[x] else xx, sum.of.weights = weights)
#     else {
#         names(weights) <- if (length(lev)) 
#             lev[x]
#         else xx
#         weights
#     }
# }
# 
# .as.category <- function (x) 
# {
#     x <- as.factor(x)
#     class(x) <- NULL
#     x
# }
# 
# .testDateTime <- function (x, what = c("either", "both", "timeVaries")) 
# {
#     what <- match.arg(what)
#     cl <- class(x)
#     if (!length(cl)) 
#         return(FALSE)
#     dc <- if (.R.) 
#         c("Date", "POSIXt", "POSIXct", "dates", "times", "chron")
#     else c("timeDate", "date", "dates", "times", "chron")
#     dtc <- if (.R.) 
#         c("POSIXt", "POSIXct", "chron")
#     else c("timeDate", "chron")
#     switch(what, either = any(cl %in% dc), both = any(cl %in% 
#         dtc), timeVaries = {
#         if ("chron" %in% cl || "Date" %in% cl || !.R.) {
#             y <- as.numeric(x)
#             length(unique(round(y - floor(y), 13))) > 1
#         } else if (.R.) length(unique(format(x, "%H%M%S"))) > 
#             1 else FALSE
#     })
# }
#     
# .all.is.numeric <- function (x, what = c("test", "vector"), extras = c(".", "NA")) 
# {
#     what <- match.arg(what)
#     old <- options(warn = -1)
#     on.exit(options(old))
#     x <- sub("[[:space:]]+$", "", x)
#     x <- sub("^[[:space:]]+", "", x)
#     xs <- x[!(x %in% c("", extras))]
#     isnum <- !any(is.na(as.numeric(xs)))
#     if (what == "test") 
#         isnum
#     else if (isnum) 
#         as.numeric(x)
#     else x
# }