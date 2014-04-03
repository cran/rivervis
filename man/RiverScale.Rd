\name{RiverScale}
\alias{RiverScale}

\title{
  River Chart Scale
}
\description{
  This plots the scale of river charts.
}
\usage{
RiverScale(length, label, riverlayout,
           loc = NA, 
           scl.col = "black", scl.lwd = 1,
           lbl.cex = 0.5, lbl.pos = 4, lbl.ofs = 0.5)
}

\arguments{
  \item{length}{
    the length of the scale. The length is defined in the same units as the river length. The function will convert this real length into a segment with the same scale as the rivers, and plot it on the river chart.
  }
  \item{label}{
    a string defining a scale label indicating the real length the scale segment represents.
  }
  \item{riverlayout}{
    the output list of \code{RiverLayout} or \code{RiverMap}.
  }
  \item{loc}{
    location of scale. One or two values in the range [0, 1] to define left and bottom margin sizes. If \code{loc = NA}, use mouse to locate the arrow. ESC to confirm.
  }

  \item{scl.col}{
    scale colour.
  }
  \item{scl.lwd}{
    scale line width.
  }
  \item{lbl.cex}{
    scale label size.
  }
  \item{lbl.pos}{
    scale label position. \code{1} for below, \code{2} for left, \code{3} for above, and \code{4} for right. See \code{par} for details.
  }
  \item{lbl.ofs}{
    scale label position offset.
  }
}

\author{
  Feng Mao
}

\seealso{
  \code{\link{RiverLayout}},
  \code{\link{RiverDraw}},
  \code{\link{RiverMap}},
  \code{\link{par}},
  \code{\link{locator}}.
}
\examples{

data(Ballinderry)

riverlayout <- RiverLayout(B.river$River,B.river$Length,
                           B.river$Parent,B.river$Position,
                           B.river$Distance, direction = -1)
RiverDraw(riverlayout)

RiverScale(2, "2 km", riverlayout, loc = c(0.8, 0.3), lbl.cex = 0.8)

# Or use mouse to allocate the river scale. Try the command below.
# RiverScale(2, "2 km", riverlayout, lbl.cex = 0.8)
}

\keyword{hplot}