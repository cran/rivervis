\name{RiverPoint}
\alias{RiverPoint}

\title{
  Points on River Charts
}
\description{
  This function plots scatter points or broken lines on the river chart.
}
\usage{
RiverPoint(site, river, distance, value, riverlayout,
           type = "l", pt.col = "grey40", pt.bg = "black", pt.pch = 20, pt.cex = 1, 
           lbl.cex = 0.7, lbl.adj = c(0.5,2), lbl.ofs = 0.5, lbl.col = "black",
           lbl.srt = 0, lbl.pos = NULL, lbl.shw = FALSE, ln.lwd = 1)
}

\arguments{
  \item{site}{
    a vector of site IDs.
  }
  \item{river}{
    a vector of river names.
  }
  \item{distance}{
    a vector of distances from sites to the river mouth.
  }
  \item{value}{
    a vector of values.
  }
  \item{riverlayout}{
    the output list of \code{RiverLayout} or \code{RiverMap}.
  }  
  \item{type}{
    type of plot. See \code{plot} for details. The default value is "l", which means "lines".
  }
  \item{ln.lwd}{
    line width.
  }
  \item{pt.col}{
    point or point border colour.
  }
  \item{pt.bg}{
    point point background(fill) colour when \code{pt.pch=21:25}.
  }
  \item{pt.pch}{
    point style.
  }
  \item{pt.cex}{
    point size.
  }
  \item{lbl.cex}{
    label size.
  }
  \item{lbl.adj}{
    label adjustment. One or two values in [0,1] for x and y (optional) adjustment.
  }
  \item{lbl.ofs}{
    label position offset.
  }
  \item{lbl.col}{
    label colour.
  }
  \item{lbl.srt}{
    label angle.
  }
  \item{lbl.pos}{
    label position. \code{1} for below, \code{2} for left, \code{3} for above, and \code{4} for right. See \code{par} for details.
  }
  \item{lbl.shw}{
    show labels (\code{TRUE}) or not (\code{FALSE}).
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

RiverPoint(B.elevation$Site, B.elevation$River, B.elevation$Distance, 
           B.elevation$Elevation, riverlayout)

RiverPoint(B.elevation$Site, B.elevation$River, B.elevation$Distance, 
           B.elevation$Elevation, riverlayout, lbl.shw=TRUE)

RiverPoint(B.elevation$Site, B.elevation$River, 
           B.elevation$Distance, NA, riverlayout, lbl.shw=TRUE)

}


\keyword{hplot}

