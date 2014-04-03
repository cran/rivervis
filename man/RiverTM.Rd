\name{RiverTM}
\alias{RiverTM}

\title{
  Tick Marks on River Charts
}
\description{
  This adds tick marks to the river chart.
}
\usage{
RiverTM(tickmark, value, riverlayout,
        range = NA, side = "L", pos = 1, 
        tm.l = 1, tm.col = "black",
        lbl.shw = TRUE, lbl.col = "black", 
        lbl.cex = 0.7, lbl.row = TRUE, label = NA)
}

\arguments{
  \item{tickmark}{
    a vector of tick mark values.
  }
  \item{value}{
    the variables which the tick marks are for.
  }
  \item{riverlayout}{
    the output list of \code{RiverLayout} or \code{RiverMap}.
  }
  \item{range}{
    bar-chart value range. A vector of two values indicating lower limit and upper limit.
  }
  \item{lbl.shw}{
    show labels of tick marks (\code{TRUE}) or not (\code{FALSE}).
  }
  
  \item{lbl.col}{
    label colour.
  }
  
  \item{lbl.cex}{
    label size.
  }
  
  \item{lbl.row}{
    show one label per row (\code{TRUE}) or not (\code{FALSE}).
  }
    
  \item{label}{
    a vector of tick mark labels.
  }
  
  \item{side}{
    position of tick marks. \code{"l"} for left and \code{"r"} for right.
  }
  \item{pos}{
    position of tick marks. \code{-1} for in and \code{1} for out.
  }
  \item{tm.l}{
    tick mark length.
  }
  \item{tm.col}{
    tick mark colour.
  }
}


\author{
  Feng Mao
}

\seealso{
  \code{\link{RiverLayout}},
  \code{\link{RiverDraw}},
  \code{\link{RiverMap}}.
  \code{\link{par}}.
}

\examples{

data(Ballinderry)

riverlayout <- RiverLayout(B.river$River,B.river$Length,
                           B.river$Parent,B.river$Position,
                           B.river$Distance, direction = -1)
RiverDraw(riverlayout)

RiverTM(c(0,2,4,6,8,10), B.bio[4:5], riverlayout, 
        pos=-1, side = "L", range = c(0,8))

RiverTM(c(0,3,8), B.bio[4:5], riverlayout, 
        pos=-1, side = "R", range = c(0,8))
}

\keyword{hplot}

