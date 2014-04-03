\name{RiverBlock}
\alias{RiverBlock}

\title{
  River Block-Chart
}
\description{
  This plots blocks to display qualitative data on the river chart.
}
\usage{
RiverBlock(site, river, distance, value, riverlayout, arrangement,
           pt.shw = FALSE, hw.rat = 1.5, h.gap = 0.05, w.gap = 0.025, 
           block.col = NA, block.lwd = 1, bd.col = "grey20",                    
           par.shw = TRUE, par.pos = 2,
           par.ofs = 1, par.cex =0.6, par.adj = c(1,0.5),
           lbl.shw = TRUE, lbl.cex = 0.7, lbl.adj = c(0.5,2), lbl.ofs = 0.5,
           lbl.col = "black", lbl.srt = 0, lbl.pos = NULL)
}

\arguments{
  \item{site}{
    a vector of site names.
  }
  \item{river}{
    a vector of rivers on which the sites are located.
  }  
  \item{distance}{
    a vector. The along-the-river distances between the sites and the mouth of the river.
  }  
  \item{value}{
    a data frame containing the qualitative variables to be shown on the block-chart.
  }  
  \item{riverlayout}{
    the output list of RiverLayout or RiverMap.
  }
  \item{arrangement}{
    a vector indicating the block number for each line.
  } 
  \item{pt.shw}{
    show location point (\code{TRUE}) or not (\code{FALSE}).
  }  
  \item{hw.rat}{
    the ratio of block height and width in the plotted diagram.
  }
  \item{h.gap}{
    vertical gap size between blocks. By default, the vertical gap is river height * 0.05 in each river chart.
  }  
  \item{w.gap}{
    horizontal gap size between blocks when there is more than one block in each line. By default, the horizontal gap is largest block width * 0.025.
  }  
  \item{block.col}{
    block colour.
  }
  \item{block.lwd}{
    block line width.
  }  
  \item{bd.col}{
    block border col.
  }  
  \item{par.shw}{
    show parameter names (\code{TRUE}) or not (\code{FALSE}).
  }  
  
  \item{par.pos}{
    parameter label position. \code{1} for below, \code{2} for left, \code{3} for above, and \code{4} for right. See \code{par} for details.
  }  
  
  \item{par.ofs}{
    parameter label position offset.
  }   
    
  \item{par.cex}{
    parameter label size. 
  }  
  
  \item{par.adj}{
    parameter label adjustment. One or two values in [0,1] for x and y (optional) adjustment. 
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

RiverBlock(B.hm$Site, B.hm$River, B.hm$Distance, B.hm[4:8], riverlayout, 
           c(1,1,1,2), block.col = c("#5381FFFF", "#7BE859FF", 
           "#FFC944FF", "#E87539FF", "#FF3931FF"))

}

\keyword{hplot}
