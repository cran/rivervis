\name{RiverLabel}
\alias{RiverLabel}

\title{
  River Labels on River Charts
}
\description{
  This adds the name labels to the plotted rivers.
}

\usage{
RiverLabel(riverlayout,
          cex = 0.7, adj = c(0, -1), srt = 90, col = "black", 
          pos = NULL, offset = 0.5, corner = "lb")
}

\arguments{
  \item{riverlayout}{
    the output list of \code{RiverLayout} or \code{RiverMap}.
  }
    
  \item{corner}{
    river label position, which can be at any of the four river chart corners. \code{"lt"} for left-top, \code{"lb"} for left-bottom, \code{"rt"} for right-top, \code{"rb"} for right-bottom.
  }
  \item{cex}{
    text size.
  }  
  \item{adj}{
    text adjustment. One or two values in [0,1] for x and y (optional) adjustment. 
  }  
  \item{srt}{
    text angle.
  }  
  \item{col}{
    text colour. 
  }
  \item{pos}{
    text position. \code{1} for below, \code{2} for left, \code{3} for above, and \code{4} for right. See \code{par} for details. It overrides \code{adj} if it is not \code{NULL}.
  }  
  \item{offset}{
    text position offset. 
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

RiverLabel(riverlayout, corner = "lt", srt = 0)

RiverLabel(riverlayout, corner = "lb")

RiverLabel(riverlayout, corner = "rt", srt = -90)
}

\keyword{hplot}

