\name{RiverBlockChart}
\alias{RiverBlockChart}

\title{
  River Block Chart Plotting
}
\description{
  This function plots a river block chart to display qualitative data without the topological structure of the river network. The function does not require the output list from \code{RiverLayout} or \code{RiverMap}.
}

\usage{
RiverBlockChart(site, river, distance, value, arrangement,
                h.gap = 0.1, w.gap = 0.25, w.gap.s = 0.1, r.gap = 0.25, 
                block.col = NA, block.lwd = 1, border.col = "grey20",
                bg.col = "lightgrey", mar = 0.1, hw.rat = 1.5,
                site.shw = TRUE, site.pos = 1, site.ofs = 1.5, 
                site.cex = 0.5, site.col = "black",site.order = "A",
                rvr.shw = TRUE, rvr.ofs = 1.5, rvr.cex = 0.7, 
                rvr.col = "black", rvr.t.b = "b", rvr.order = NA,
                par.shw = TRUE, par.pos = 2, par.ofs = 1, 
                par.cex = 0.6, par.adj = c(1, 0.5), par.col = "black")
}

\arguments{
  \item{site}{
    a character vector of site names.
  }
  \item{river}{
    a vector of rivers to which the sites belong.
  }  
  \item{distance}{
    a vector. The along-the-river distance between the site and the mouth of the river.
  }
  \item{value}{
    a data frame containing the qualitative variables to be shown on the river block-chart.
  }  
  \item{arrangement}{
    a vector indicating the block number for each line.
  }  
  \item{h.gap}{
    vertical gap size between blocks in the plot. By default, the vertical gap is block height * 0.1
  }  
  \item{w.gap}{
    horizontal gap size between sites. By default, the horizontal gap is largest block width * 0.25.
  } 
  \item{w.gap.s}{
    horizontal gap size between small blocks when there are more than one block in each line. By default, the horizontal gap is largest block width * 0.1.
  } 
  \item{r.gap}{
    horizontal gap size between rivers. By default, gap between rivers is block width * 0.25.
  }
  \item{block.col}{
    a vector of block colours. The length of vector should be as the same as the the number of levels of \code{value}.
  }  
  \item{block.lwd}{
    a value of block line width.
  }
  \item{border.col}{
    a value of block border colour.
  }  
  \item{bg.col}{
    a value of river background colour.
  }  
  \item{mar}{
    a value of smallest margin size.
  }  
  \item{hw.rat}{
    the ratio of block height and block width. 
  } 
  \item{site.shw}{
    show site names (\code{TRUE}) or not (\code{FALSE}).
  } 
  \item{site.pos}{
    site position. \code{1} for below, \code{2} for left, \code{3} for above, and \code{4} for right. See \code{par} for details.
  } 
  \item{site.ofs}{
    site position offset.
  } 
  \item{site.cex}{
    site name size.
  }   
  \item{site.col}{
    site colour.
  } 
  \item{site.order}{
    order of sites within each river. Alphabetical order (\code{"A"}), river flow left (\code{"L"}), river flow right (\code{"R"}).
  } 
  \item{rvr.shw}{
    show river labels (\code{TRUE}) or not (\code{FALSE}).
  } 
  \item{rvr.t.b}{
    location of river label. \code{"t"} for top (above) and \code{"b"} for below (bottom).
  } 
  \item{rvr.ofs}{
    river label position offset.
  } 
  \item{rvr.cex}{
    river label size.
  } 
  \item{rvr.col}{
    river label colour.
  } 
  \item{rvr.order}{
    order of rivers. Alphabetical order (\code{NA}) or a vector of custom order.
  } 
  \item{par.shw}{
    show parameter labels (\code{TRUE}) or not (\code{FALSE}).
  } 
  \item{par.pos}{
    parameter label position. \code{1} for below, \code{2} for left, \code{3} for above, and \code{4} for right. See \code{par} for details. It overrides \code{par.adj} if par.pos is not \code{NULL}.
  }     
  \item{par.ofs}{
    parameter label position offset.
  } 
  \item{par.cex}{
    parameter label size.
  } 
  \item{par.col}{
    parameter label colour.
  } 
  \item{par.adj}{
    parameter label adjustment. One or two values in [0,1] for x and y (optional) adjustment. 
  }  

}

\author{
  Feng Mao
}


\seealso{
  \code{\link{par}}.
}

\examples{

data(Ballinderry)

RiverBlockChart(B.hm$Site, B.hm$River, B.hm$Distance, B.hm[4:8],  
                c(1,1,1,2), mar = 0.15, 
                block.col = c("#5381FFFF", "#7BE859FF", "#FFC944FF", "#E87539FF", "#FF3931FF"))

RiverBlockChart(B.hm$Site, B.hm$River, B.hm$Distance, B.hm[4:8],  
                rvr.order = c("Rock", "KillymoonClaggan", 
                "Ballinderry", "Ballymully", "Kildress", 
                "Kingsmill", "Lissan", "Tulnacross"),
                site.order = "R", mar = 0.15, c(1,1,1,2), 
                block.col = c("#5381FFFF", "#7BE859FF", 
                "#FFC944FF", "#E87539FF", "#FF3931FF"))

RiverBlockChart(B.hm$Site, B.hm$River, B.hm$Distance, B.hm[4:8],  
                rvr.order = c("Rock", "KillymoonClaggan", 
                "Ballinderry", "Ballymully", "Kildress", 
                "Kingsmill", "Lissan", "Tulnacross"),
                site.order = "L", mar = 0.15, c(1,1,1,2), 
                block.col = c("#5381FFFF", "#7BE859FF", 
                "#FFC944FF", "#E87539FF", "#FF3931FF"))
}

\keyword{hplot}

