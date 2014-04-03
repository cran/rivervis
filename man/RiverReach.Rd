\name{RiverReach}
\alias{RiverReach}

\title{
  River Reach Plotting
}

\description{
  This highlights river reaches on the river chart.
}

\usage{
RiverReach(reach, river, from, to, group, style, riverlayout,
           rea.pos = NA, rea.col = "lightblue", rea.lty = 1,
           rea.lwd = 1, rea.den = NULL, bd.col = "black",
           ln.shw = T, ln.col = "grey40", ln.lty = 3, ln.lwd = 1,
           pt.shw = T, pt.col = "black", pt.pch = 20,
           pt.bg = "black", pt.cex = 1, pt.lwd = 1)
}

\arguments{
  \item{reach}{
    a vector of reach names.
  }
  \item{river}{
    a vector of rivers to which the reaches belong.
  }
  \item{from}{
    a numeric vector of starting points.
  }
  \item{to}{
    a numeric vector of ending points.
  }
  \item{group}{
    a vector of reach group names. This indicates to which group the reaches belong.
  }
  \item{style}{
    a vector of reach styles. The value of "style" denotes the location of reach lines. Especially, \code{0} denotes "on axis" and \code{99} means "the reach is presented as a band rather than a line". 
  }
  \item{riverlayout}{
    the output list of \code{RiverLayout} or \code{RiverMap}.
  }

  \item{rea.pos}{
    a vector of absolute positions of lines. The values range in [0,1].
  }
  \item{rea.col}{
    line colour.
  }
  \item{rea.lty}{
    line style.
  }
  \item{rea.lwd}{
    line width.
  }
  \item{rea.den}{
    the density of shading lines, in lines per inch. See \code{rect}.
  }
  \item{bd.col}{
    colour of river chart frames.
  }
  \item{ln.shw}{
    show lead lines (\code{TRUE}) or not (\code{FALSE}).
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
  
}

\author{
  Feng Mao
}


\seealso{
  \code{\link{RiverLayout}},
  \code{\link{RiverDraw}},
  \code{\link{RiverMap}},
  \code{\link{par}}.
}

\examples{

data(Ballinderry)

riverlayout <- RiverLayout(B.river$River,B.river$Length,
                           B.river$Parent,B.river$Position,
                           B.river$Distance, direction = -1)
RiverDraw(riverlayout)

RiverReach(B.reach$Reach, B.reach$River, B.reach$From, 
           B.reach$To, B.reach$Group, B.reach$Style, riverlayout,
           rea.lwd = 5)

RiverReach(B.reach$Reach, B.reach$River, B.reach$From, 
           B.reach$To, B.reach$Group, c(99,1,2,3,4,5), riverlayout,
           rea.lwd = 5)

RiverReach(B.reach$Reach, B.reach$River, B.reach$From, 
           B.reach$To, B.reach$Group, B.reach$Style, riverlayout,
           rea.lwd = 4, rea.lty = c(3,1,1,1),
           rea.col = c("#FFC944FF", "#FFC944FF", "#E87539FF", "light blue"))
}

\keyword{hplot}

