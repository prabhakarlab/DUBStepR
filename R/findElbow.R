#' @title
#' Find the Elbow in a Curve
#'
#' @description
#' This utility function finds the elbow in a curve which is concave
#' relative to a line drawn between the first and last points.
#' The elbow is defined as the point with the greatest
#' orthogonal distance from that line.
#'
#' @param y Numeric vector of y values for the curve.
#'@param ylab Y-axis label.
#' @param plot Logical. Should a plot be made?
#'
#' @param returnIndex Logical. Should the return value
#' be the index of the elbow point?
#'
#' @return If \code{returnIndex = TRUE}, the index of
#' the elbow point.  If \code{returnIndex = FALSE},
#' a data frame containing an index values (x),
#' the y values passed to the function, and the
#' the orthogonal distances of the y values from
#' the line connecting the first and last points.
#' \code{which.max(data_frame_name$dist)} will give the index of
#' the elbow point.
#'
#' @references The concept of this function is based on the
#' clever idea in the
#' Stackoverflow post at stackoverflow.com/a/2022348/633251
#' and relies on code posted at
#' paulbourke.net/geometry/pointlineplane/pointline.r
#' (referenced in the SO post).  Minor modifications
#' to the code were made to that code in order to vectorize it.
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @importFrom stats lm coef
#'
#' @section Warning:
#' This function makes some simple checks that the data is concave as defined above.  Even so, it may give  answers in some cases that are not valid.  Please check on typical data that you encounter to verify that it works in your cases.
#'
#'

findElbow <- function(y, ylab = "y values", plot = FALSE, returnIndex = TRUE) {

    # The following helper functions were found at
    # paulbourke.net/geometry/pointlineplane/pointline.r
    # via the SO reference below.  The short segment check
    # was modified to permit vectorization.

    ##========================================================
    ##
    ##  Credits:
    ##  Theory by Paul Bourke http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
    ##  Based in part on C code by Damian Coventry Tuesday, 16 July 2002
    ##  Based on VBA code by Brandon Crosby 9-6-05 (2 dimensions)
    ##  With grateful thanks for answering our needs!
    ##  This is an R (http://www.r-project.org) implementation by Gregoire Thomas 7/11/08
    ##
    ##========================================================

    distancePointLine <- function(x, y, slope, intercept) {
        ## x, y is the point to test.
        ## slope, intercept is the line to check distance.
        ##
        ## Returns distance from the line.
        ##
        ## Returns 9999 on 0 denominator conditions.
        x1 <- x-10
        x2 <- x+10
        y1 <- x1*slope+intercept
        y2 <- x2*slope+intercept
        distancePointSegment(x,y, x1,y1, x2,y2)
    }

    distancePointSegment <- function(px, py, x1, y1, x2, y2) {
        ## px,py is the point to test.
        ## x1,y1,x2,y2 is the line to check distance.
        ##
        ## Returns distance from the line, or if the intersecting point on the line nearest
        ## the point tested is outside the endpoints of the line, the distance to the
        ## nearest endpoint.
        ##
        ## Returns 9999 on 0 denominator conditions.
        lineMagnitude <- function(x1, y1, x2, y2) sqrt((x2-x1)^2+(y2-y1)^2)
        ans <- NULL
        ix <- iy <- 0   # intersecting point
        lineMag <- lineMagnitude(x1, y1, x2, y2)
        if(any(lineMag < 0.00000001)) { # modified for vectorization by BAH
            #warning("short segment")
            #return(9999)
            warning("At least one line segment given by x1, y1, x2, y2 is very short.")
        }
        u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
        u <- u / (lineMag * lineMag)
        if(any(u < 0.00001) || any(u > 1)) { # BAH added any to vectorize
            ## closest point does not fall within the line segment, take the shorter distance
            ## to an endpoint
            ix <- lineMagnitude(px, py, x1, y1)
            iy <- lineMagnitude(px, py, x2, y2)
            if(ix > iy)  ans <- iy
            else ans <- ix
        } else {
            ## Intersecting point is on the line, use the formula
            ix <- x1 + u * (x2 - x1)
            iy <- y1 + u * (y2 - y1)
            ans <- lineMagnitude(px, py, ix, iy)
        }
        ans
    }

    # End of helper functions by PB

    ### Now for the actual findElbow function!

    # Find the elbow using the method described in
    # stackoverflow.com/a/2022348/633251
    # but translated to R (see above).

    # Add an index to argument values for easy plotting
    DF <- data.frame(x = 1:length(y), y = y)
    fit <- lm(y ~ x, DF[c(1,nrow(DF)),]) # 2 point 'fit'
    m <- coef(fit)[2]
    b <- coef(fit)[1]

    # Check to make sure the data is concave as described
    # in the documentation, as arbitrary trends could give
    # misleading answers.  The following approach simply
    # checks to make sure all values are either above or
    # below the reference line.  This allows the values
    # to vary quite a bit and still return an answer.

    concave <- FALSE
    use <- 2:(nrow(DF)-1)
    refpts <- m*DF$x[use] + b
    if (all(refpts > DF$y[use]) | all(refpts < DF$y[use]))
        concave <- TRUE
    else {
        warning("Your curve doesn't appear to be concave")
    }
    # if (!concave) stop("Your curve doesn't appear to be concave")

    # Calculate the orthogonal distances
    use <- 2:(nrow(DF)-1)
    elbowd <- distancePointLine(DF$x[use], DF$y[use], coef(fit)[2], coef(fit)[1])
    DF$dist <- c(NA, elbowd, NA) # first & last points don't have a distance

    if (plot) {
        edm <- which.max(DF$dist)
        plot(DF[,1:2], type = "b", xlab = "No. of steps", ylab = ylab,
             main = "Looking for the Elbow")
        graphics::segments(DF$x[1], DF$y[1],
                 DF$x[nrow(DF)], DF$y[nrow(DF)], col = "red")
        graphics::points(DF$x[edm], DF$y[edm], cex = 1.5, col = "red")
        graphics::points(DF$x[edm], DF$y[edm], pch = 20)
    }

    if (returnIndex) return(which.max(DF$dist))
    if (!returnIndex) return(DF)

}

