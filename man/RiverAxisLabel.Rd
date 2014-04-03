\name{RiverAxisLabel}
\alias{RiverAxisLabel}

\title{
  River Axis Labels
}
\description{
  This adds left or right axis labels to the river chart.
}
\usage{
RiverAxisLabel(label, riverlayout,
               cex = 0.7, adj = c(0.5, -2), srt = 90,
               col = "black", pos = NULL, offset = 0.5,
               side = "L", mainonly = TRUE)
}

\arguments{
  \item{label}{
    the axis label to be shown on the river chart.
  }
  \item{riverlayout}{
    the output list of \code{RiverLayout}.
  }
  \item{cex}{
    text size.
  }
  \item{adj}{
    text adjustment. One or two values in the range [0,1] for x and y (optional) adjustment. 
  }
  \item{srt}{
    text angle.
  }
  \item{col}{
    text colour.
  }
  \item{pos}{
    text position. \code{1} for below, \code{2} for left, \code{3} for above, and \code{4} for right. See \code{par} for details.
  }
  \item{offset}{
    text position offset.
  }
  \item{side}{
    left (\code{"L"}) or right (\code{"R"}) axis.
  }
  \item{mainonly}{
    the axis title is only shown for the main stream only (\code{"TRUE"}) or not (\code{"FALSE"}).
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

RiverAxisLabel("Score Value", riverlayout, adj = c(0.5, 1))

}

\keyword{hplot}
