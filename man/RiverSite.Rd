\name{RiverSite}
\alias{RiverSite}

\title{
  Site of Interest
}
\description{
  This plots sites of interest on the river chart.
}
\usage{
RiverSite(site, river, distance, group, riverlayout,
              pt.pch = 21, pt.col = NA, pt.bg = "red", pt.cex = 1,
              lbl.cex = 0.5, lbl.srt = 0, lbl.adj = c(0.5,2), lbl.col = "black",
              lbl.pos = 1, lbl.ofs = 0.5, lbl.shw = TRUE)
}

\arguments{
  \item{site}{
    a character vector of site names.
  }
  \item{river}{
    a vector of rivers on which the sites are located.
  }  
  \item{distance}{
    a vector. The along-the-river distance between the site and the mouth of the river.
  }  
  \item{group}{
    a vector. Group names of river locations.
  }  
  \item{riverlayout}{
    the output list of \code{RiverLayout} or \code{RiverMap}.
  }
  \item{pt.pch}{
    point style.
  }  
  \item{pt.col}{
    point border colour.
  }  
  \item{pt.bg}{
    point background(fill) colour when \code{pt.pch=21:25}.
  }    
  \item{pt.cex}{
    point size.
  }
  \item{lbl.shw}{
    show labels (\code{TRUE}) or not (\code{FALSE}).
  }  
  \item{lbl.cex}{
    label size.
  }  
  \item{lbl.srt}{
    label angle.
  }  
  \item{lbl.adj}{
    label adjustment. One or two values in [0,1] for x and y (optional) adjustment. 
  }  
  \item{lbl.col}{
    label colour.
  }  
  \item{lbl.pos}{
    label position. \code{1} for below, \code{2} for left, \code{3} for above, and \code{4} for right. See \code{par} for details.
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
  \code{\link{par}}.
}


\examples{

data(Ballinderry)

riverlayout <- RiverLayout(B.river$River,B.river$Length,
                           B.river$Parent,B.river$Position,
                           B.river$Distance, direction = -1)
RiverDraw(riverlayout)

RiverSite(B.soi$SOI, B.soi$River, B.soi$Distance, B.soi$Group, riverlayout, 
          pt.bg = c("red","green","yellow", "blue"), lbl.shw = FALSE)

RiverSite(B.soi$SOI, B.soi$River, B.soi$Distance, 1, riverlayout, pt.bg = "red")
}

\keyword{hplot}