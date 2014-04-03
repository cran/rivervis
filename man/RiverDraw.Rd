\name{RiverDraw}
\alias{RiverDraw}

\title{
  River Chart Plotting
}
\description{
This plots the river chart according to the output list provided by \code{RiverLayout}.
}

\usage{
RiverDraw(riverlayout, 
          bd.col = "black", 
          ln.col = "grey40", ln.lty = 3, ln.lwd = 1, bg.col = "grey80",
          pt.shw = TRUE, pt.col = "black", pt.pch = 20,
          pt.bg = "black", pt.cex = 1, pt.lwd = 1,
          mar.t = 0.05, mar.b = 0.05, mar.l = 0.2, mar.r = 0.1)
}

\arguments{
  \item{riverlayout}{
    the output list of \code{RiverLayout}.
  }
  
  \item{bd.col}{
    border colour.
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
    anchor point style.
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
    top margin size.
  }

  \item{mar.b}{
    bottom margin size.
  }

  \item{mar.l}{
    left margin size.
  }

  \item{mar.r}{
    right margin size.
  }

}



\author{
  Feng Mao
}



\seealso{
  \code{\link{RiverLayout}},
  \code{\link{RiverMap}},
  \code{\link{par}}.
}

\examples{

data(Ballinderry)

# River flows right
riverlayout <- RiverLayout(B.river$River,B.river$Length,
                           B.river$Parent,B.river$Position,
                           B.river$Distance, direction = -1)
RiverDraw(riverlayout)

# River flows left
riverlayout.left <- RiverLayout(B.river$River,B.river$Length,
                                B.river$Parent,B.river$Position,
                                B.river$Distance)

RiverDraw(riverlayout.left)
}

\keyword{hplot}

