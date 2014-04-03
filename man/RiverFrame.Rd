\name{RiverFrame}
\alias{RiverFrame}

\title{
River Frame Plotting
}
\description{
This plots river frames, lead lines and archor points.
}

\usage{
RiverFrame(riverlayout,
           ln.shw = T, ln.col = "grey40", ln.lty = 3, ln.lwd = 1,
           pt.shw = T, pt.col = "black", pt.pch = 20, pt.bg = "black",
           pt.cex = 1, pt.lwd = 1,
           bd.shw = T, bd.col = "black")
}

\arguments{
  \item{riverlayout}{
    the output list of \code{RiverLayout}.
  }
  \item{bd.col}{
    border colour.
  }
  \item{bd.shw}{
    show boarders (\code{TRUE}) or not (\code{FALSE}).
  }
  \item{ln.col}{
    lead line colour.
  }
  \item{ln.shw}{
    show lead lines (\code{TRUE}) or not (\code{FALSE}).
  }
  \item{ln.lty}{
    lead line style.
  }
  \item{ln.lwd}{
    lead line width.
  }
  \item{pt.shw}{
    show anchor points (\code{TRUE}) or not (\code{FALSE}). Anchor points represent the locations of the river mouths.
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

