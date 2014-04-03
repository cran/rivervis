\name{RiverDirection}
\alias{RiverDirection}

\title{
  River Direction Arrow
}
\description{
  This plots river flow direction arrow on river charts.
}
\usage{
RiverDirection(riverlayout, 
               loc = NA, 
               arw.length = 0.05, arw.lty = 1, 
               arw.lwd = 1, arw.angle = 30, arw.col = "black",
               label = "Flow direction", lbl.cex = 0.5, lbl.pos = 4, lbl.ofs = 0.5)
}

\arguments{
  \item{riverlayout}{
    the output list of \code{RiverLayout} or \code{RiverMap}.
  }
  \item{loc}{
    location of arrow. One or two values in the range [0, 1] for left and bottom margin sizes. If \code{loc = NA}, use mouse to locate the arrow. ESC to confirm.
  }
  \item{arw.length}{
    arrow length.
  }  
  \item{arw.lty}{
    arrow line style.
  }
  \item{arw.lwd}{
    arrow line width.
  }    
  \item{arw.angle}{
    arrow head angle.
  }
  \item{arw.col}{
    arrow colour.
  }
  \item{label}{
    label of the arrow.
  }  
  \item{lbl.cex}{
    label size.
  }
  \item{lbl.pos}{
    label position.
  }
  \item{lbl.ofs}{
    label position offset.
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
  \code{\link{locator}},
  \code{\link{arrows}}.
}

\examples{

data(Ballinderry)

riverlayout <- RiverLayout(B.river$River,B.river$Length,
                           B.river$Parent,B.river$Position,
                           B.river$Distance, direction = -1)
RiverDraw(riverlayout)

RiverDirection(riverlayout, loc = c(0.8, 0.3), arw.length = 0.03, lbl.cex = 0.8)

# Or use mouse to allocate the flow direction sign. Try the command below.
# RiverDirection(riverlayout, arw.length = 0.03, lbl.cex = 0.8)

}

\keyword{hplot}

