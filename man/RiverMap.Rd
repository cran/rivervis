\name{RiverMap}
\alias{RiverMap}

\title{
  River Layout Coordinates Calculation and Plotting
}
\description{
  This function calculates plotting coordinates for rivers and draws the river chart accordingly.
}
\usage{
RiverMap(river, length, parent, position, distance, 
         row = NA, direction = 1, margin = 0.5,
         bd.col = "black", ln.col = "grey40", 
         ln.lty = 3, ln.lwd = 1, bg.col = "grey80",
         pt.shw = TRUE, pt.col = "black", 
         pt.pch = 20, pt.bg = "black", pt.cex = 1, pt.lwd = 1,
         mar.t = 0.05, mar.b = 0.05, mar.l = 0.2, mar.r = 0.1)
}

\arguments{
 \item{river}{
    a vector of river names.
  }
  \item{length}{
    a vector of river lengths.
  }
  \item{parent}{
    a vector of river parents. The parent of a river is the river into which it flows. 
    The parent of the main stream is NA.
  }
  \item{position}{
    a vector of river positions. The river position indicates its position relative to its parent - whether it is a left bank river, right bank river or main stream. The left bank river is on the left when looking downstream of its parent. The right bank river is on the right when looking downstream of its parent.
    The value of position can be "R", "L" or "M". For the main stream, the value is "M".
  }
  \item{distance}{
    a vector of distances denotes the distance between the mouths of each river and the mouths of each river's parent. 
  }
  \item{row}{
    a vector of row numbers. The main stream is on row 0. In the river chart, rivers with negative row numbers are plotted below the main stream while rivers with positive row numbers are plotted above the main stream. If a value for row is provided, the rivers will be plotted according to the provided row numbers. If a value for row is not provided, a vector of best fit row numbers will be calculated before plotting.
  }
  
  \item{direction}{
    a value. In the river chart, rivers flow from right to left (direction = 1), or from left to right (direction = -1). By default, direction = 1.
  }
  \item{margin}{
    a value. The margin height between rivers in the topological plot. By default, margin = 0.5 and margin height is 0.5 times the river height in the river chart.
    }

  \item{bd.col}{
    river border colour.
  }
  \item{bg.col}{
    background colour.
  }
  \item{ln.col}{
    lead line colour.
  }
  \item{ln.lty}{
    lead line style.
  }
  \item{ln.lwd}{
    lead line width.
  }
  \item{pt.shw}{
    show anchor point (\code{TRUE}) or not (\code{FALSE}). Anchor points represent the locations of the river mouths.
  }
  \item{pt.col}{
    anchor point colour.
  }
  \item{pt.pch}{
    anchor point character.
  }
  \item{pt.bg}{
    anchor point background(fill) colour when \code{pch=21:25}.
  }
  \item{pt.cex}{
    anchor point size.
  }
  \item{pt.lwd}{
    anchor point border width.
  }
  \item{mar.t}{
    top margin size. This ranges in [0, 1] where 1 is the total height of the diagram region.
  }

  \item{mar.b}{
    bottom margin size.  This ranges in [0, 1] where 1 is the total height of the diagram region.
  }

  \item{mar.l}{
    left margin size.  This ranges in [0, 1] where 1 is the total width of the diagram region.
  }

  \item{mar.r}{
    right margin size. This ranges in [0, 1] where 1 is the total width of the diagram region.
  }

}


\value{
  The \code{RiverMap} returns a list containing data for river map drawing, and plots the river map accordingly. The output list can be used for further plotting. The output list includes,
  \item{riverdata}{a data frame. This contains input vectors \code{river}, \code{length}, \code{parent}, \code{position} and \code{distance}. It also includes calculated x-coordinates of river mouths (\code{rmouth}) and sources (\code{rsource}), defined in the same units as the inputs \code{length} and \code{distance}. The last included vector is the row number  for each river (\code{row}), in which the main stream has a fixed row number of 0.}
  \item{H.MAX}{the number of rows.}
  \item{H.SIZE}{the height of each row in the topological plot.}
  \item{W.MAX}{the width of river layout, in the same units as \code{length} and \code{distance}.}
  \item{W.SIZE}{the reciprocal of \code{W.MAX}.}
  \item{X1}{normalised x-coordinate of river mouths.}  
  \item{X2}{normalised x-coordinate of river sources.}
  \item{Y}{normalised y-coordinate of rivers.}
  \item{direction}{flow direction. Flow from right to left (\code{direction = 1}), or from left to right (\code{direction = -1}).}

}

\author{
  Feng Mao
}

\seealso{
  \code{\link{RiverLayout}},
  \code{\link{RiverDraw}},
  \code{\link{par}}.
}

\examples{

data(Ballinderry)

riverlayout <- RiverLayout(B.river$River,B.river$Length,
                           B.river$Parent,B.river$Position,
                           B.river$Distance, direction = -1)
str(riverlayout)

RiverMap(B.river$River,B.river$Length,B.river$Parent,
         B.river$Position, B.river$Distance)[[1]]

RiverMap(B.river$River,B.river$Length,B.river$Parent,
         B.river$Position, B.river$Distance, 
         row = c(5,-1,6,3,-4,2,-6,7), direction = -1)
}

\keyword{hplot}
