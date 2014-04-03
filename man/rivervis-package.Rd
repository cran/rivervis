\name{rivervis-package}
\alias{rivervis}
\docType{package}

\title{
  rivervis Package: River Visualisation Tool
}
\description{
  The \bold{rivervis} package is designed to visualise river ecosystem data.
}

\details{
  In general, the \bold{rivervis} package draws two types of diagrams - river charts and river block charts. River charts can present points, lines, bars and blocks in relation to the topological structure of the river network. River block charts show qualitative data without the river network structure. It is recommended to run the examples below and in each function manual. The \bold{rivervis} package contains 15 functions in total.
  \describe{
    \item{\code{RiverLayout}}{
      This calculates best fit plotting coordinates for rivers to be shown on river charts. The output is a list, which is used when plotting the river chart and the information on the river chart. It provides an opportunity to change the coordinates and other plotting parameters before actually plotting.
    }
    \item{\code{RiverDraw}}{
      This plots the river charts according to the output list of \code{RiverLayout}.
    }
    \item{\code{RiverMap}}{
      This can be understood as a combination of \code{RiverLayout} and \code{RiverDraw}. It not only calculates best fit plotting coordinates for rivers to be shown on river charts, but also plots the river charts according to the calculated coordinates. This implies that the coordinates cannot be changed before river chart plotting.
    }
    \item{\code{RiverFrame}}{
      This plots river frames, lead lines and archor points. 
    }    
    \item{\code{RiverPoint}}{
      This plots points or broken lines on the river chart.
    }
    \item{\code{RiverBar}}{
      This plots bars for quantitative data on the river chart.
    }
    \item{\code{RiverBlock}}{
      This plots blocks for qualitative data on the river chart.
    }
    \item{\code{RiverSite}}{
      This plots sites of interest on the river chart.
    }
    \item{\code{RiverLabel}}{
      This adds the name labels to the plotted rivers.
    }
    \item{\code{RiverTM}}{
      This adds tick marks to the river chart.
    }
    \item{\code{RiverAxisLabel}}{
      This adds left or right axis labels to the river chart.
    }
    \item{\code{RiverReach}}{
      This highlights river reaches on the river chart.
    }
    \item{\code{RiverDirection}}{
      This adds a flow direction arrow on the river chart.
    }
    \item{\code{RiverScale}}{
      This adds a plotting scale on the river chart.
    }
    \item{\code{RiverBlockChart}}{
      This function plots a river block chart for qualitative data without the topological structure of the river network. The function does not require the output list from \code{RiverLayout} or \code{RiverMap}.
    }
  }
}

\author{
  Feng Mao,
  Yichuan Shi, and
  Keith Richards
}


\examples{

data(Ballinderry)

# Example 1

riverlayout <- RiverMap(B.river$River,B.river$Length,
                           B.river$Parent,B.river$Position,
                           B.river$Distance, direction = -1)

RiverReach(B.reach$Reach[1], B.reach$River[1], B.reach$From[1], 
            B.reach$To[1], B.reach$Group[1], B.reach$Style[1],
            riverlayout, rea.col = c("#51B0A8FF"))

RiverPoint(B.elevation$Site, B.elevation$River, 
           B.elevation$Distance, B.elevation$Elevation, riverlayout)

RiverBar(B.bio$Site, B.bio$River, B.bio$Distance, B.bio[4:5], riverlayout, 
         range = c(0,8), bar.col = c("#FFC944FF", "#FF3931FF"))

RiverLabel(riverlayout, offset = -1, corner = "lt", srt = 0)

RiverTM(c(0,2,4,6,8,10), B.bio[4:5], riverlayout, 
        pos=-1, side = "L", range = c(0,8))

RiverTM(c(0, 50, 100, 150, 200, 250), B.elevation[3], riverlayout, 
        pos=-1, side = "R", range = c(0,250))

RiverAxisLabel("ASPT score", riverlayout, adj = c(0.5, -3))

RiverAxisLabel("Elevation (m)", riverlayout, side = "R", srt = 270, adj = c(0.5, -4))

title(main = "Ballinderry River", line = -2)

legend(x = 0.8, y = 0.45, inset=0.05, 
       title = "Legend", c("Elevation", "ASPT Spring", "ASPT Autumn", "Upper Ballinderry"), 
       lty = c(1,-1,-1, -1), 
       pch = c(-1, 22,22, 22),
       col= c("grey40", "black", "black", "black"),
       border = c(NA, rep("black", 4)),
       pt.bg = c(NA, "#FFC944FF", "#FF3931FF", "#51B0A8FF"),
       pt.cex = c(NA, 2,2,2.5),
       cex = 0.8)

RiverScale(2, "2 km", riverlayout, loc = c(0.8, 0.2), lbl.cex = 0.8)

RiverDirection(riverlayout, loc = c(0.8, 0.15), arw.length = 0.03, lbl.cex = 0.8)

# Example 2

riverlayout <- RiverLayout(B.river$River,B.river$Length,
                           B.river$Parent,B.river$Position,
                           B.river$Distance, direction = -1)
RiverDraw(riverlayout)

RiverBlock(B.hm$Site, B.hm$River, B.hm$Distance, B.hm[4:8], riverlayout, 
           c(1,1,1,2), lbl.shw = TRUE,
           block.col = c("#5381FFFF", "#7BE859FF", "#FFC944FF", "#E87539FF", "#FF3931FF"))

RiverSite(B.soi$SOI[21:27], B.soi$River[21:27], B.soi$Distance[21:27], B.soi$Group[21:27], 
          riverlayout, pt.pch = c(25, 24), lbl.shw = FALSE, pt.bg = "blue")

RiverLabel(riverlayout, offset = -1, corner = "lt", srt = 0)

title(main = "Ballinderry River Hydromorphological Conditions", line = -2)

legend(x = 0.8, y = 0.43, inset=0.05, 
       title = "Legend", 
       c("High", "Good", "Moderate", "Poor", "Bad", "Other left tribs", "Other right tribs"), 
       border = c(rep("black", 5), NA,NA),
       pch = c(rep(NA, 5), 25, 24),
       fill = c("#5381FFFF", "#7BE859FF", "#FFC944FF", "#E87539FF", "#FF3931FF", NA, NA),
       pt.bg = c(rep(NA, 5), "blue", "blue"),
       cex = 0.8)

RiverScale(2, "2 km", riverlayout, loc = c(0.8, 0.1), lbl.cex = 0.8)

RiverDirection(riverlayout, arw.length = 0.03, loc = c(0.8, 0.05), lbl.cex = 0.8)


# Example 3

RiverBlockChart(B.hm$Site, B.hm$River, B.hm$Distance, B.hm[4:8],  
                c(1,1,1,2), mar = 0.15, 
                block.col = c("#5381FFFF", "#7BE859FF", "#FFC944FF", "#E87539FF", "#FF3931FF"))

legend(x = 0.88, y = 0.6, inset=0.05, 
       title = "Quality Levels", c("High", "Good", "Moderate", "Poor", "Bad"), 
       border = rep("black", 5), 
       fill = c("#5381FFFF", "#7BE859FF", "#FFC944FF", "#E87539FF", "#FF3931FF"),
       cex = 0.8)

title(main = "Ballinderry River Hydromorphological Summary", line = -4)

}

\keyword{hplot}