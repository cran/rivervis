\name{Ballinderry}
\alias{B.bio}
\alias{B.elevation}
\alias{B.hm}
\alias{B.reach}
\alias{B.river}
\alias{B.soi}
\docType{data}
\title{
Ballinderry River Basin Dataset
}
\description{
A ballinderry river basin dataset for demonstrating purposes.
}
\usage{data(Ballinderry)}
\details{
  This contains the following 6 datasets.   
    \describe{
        \item{\code{B.bio}}{  
          ASPT Scores for Ballinderry River Basin in Spring and Autumn, 2009. It is a data frame with 90 observations on the following 4 variables.
          \describe{
            \item{\code{Site}}{Elevation site ID}
            \item{\code{River}}{Rivers on which the elevation sites are located}
            \item{\code{Distance}}{The along-the-river distance between the elevation sites and the mouth of the river.}
            \item{\code{Elevation}}{A numeric vector of elevation values}
          }
        }
        \item{\code{B.elevation}}{
          This dataset can be used to plot elevation profile of Ballinderry Rivers. It is a data frame with 90 observations on the following 4 variables.
          \describe{
            \item{\code{Site}}{Elevation site ID}
            \item{\code{River}}{Rivers on which the elevation sites are located}
            \item{\code{Distance}}{The along-the-river distance between the elevation sites and the mouth of the river.}
            \item{\code{Elevation}}{A numeric vector of elevation values}
          }
        }
        \item{\code{B.hm}}{
          Selected hydromorphological results from RHAT. The hydromorphological variables are ordinary factors which have five grades: High, Good, Moderate, Poor and Bad. It is a data frame with 17 observations on the following 8 variables.
          \describe{
            \item{\code{Site}}{Monitorings sites}
            \item{\code{River}}{Rivers on which the monitoring sites are located}
            \item{\code{Distance}}{The along-the-river distance between the sites and the mouth of the river.}    
            \item{\code{ChanMorph}}{Channel morphology condition}
            \item{\code{Substrate}}{Substrate condition}
            \item{\code{ChanFlow}}{Channel flow condition}
            \item{\code{BankVegLeft}}{Left bank vegetation condition}
            \item{\code{BankVegRight}}{Right bank vegetation condition}
          }
        }    
        \item{\code{B.reach}}{
          Selected Ballinderry River reaches. It is a data frame with 14 observations on the following 5 variables.
          \describe{
            \item{\code{Reach}}{Reach names.}
            \item{\code{River}}{Rivers on which the monitoring sites are located.}
            \item{\code{From}}{A numeric vector of starting points of reaches.}
            \item{\code{To}}{A numeric vector of ending points of reaches.}
            \item{\code{Group}}{A vector of reach group names. This indicates to which group the reaches belong.}
            \item{\code{Style}}{A vector of reach styles and the location of reach lines.}    
          }
        }    
        \item{\code{B.river}}{
          Main rivers and tributaries in Ballinderry Basin for \code{RiverMap} and \code{RiverLayout}. It is a data frame with 8 observations on the following 5 variables.
          \describe{
            \item{\code{River}}{River names.}
            \item{\code{Length}}{Length of rivers.}
            \item{\code{Position}}{Relative relations between rivers and their parent rivers.}
            \item{\code{Parent}}{Parent rivers.}
            \item{\code{Distance}}{Distance between the mouths of each river and the mouths of each river's parent.}
          }
        }    
        \item{\code{B.soi}}{
          This dataset provides 27 sites in Ballinderry River Basin. The sites have four types: monitoring sites, towns, conjunctions of left and right tributaries. It is a data frame with 27 observations on the following 4 variables.
          \describe{
            \item{\code{SOI}}{Sites of interest.}
            \item{\code{River}}{Rivers on which the sites are located.}
            \item{\code{Distance}}{The along-the-river distance between the sites and the mouth of the river.}
            \item{\code{Group}}{Groups of the sites.}
          }
        }       
    }
    


  
}

\source{
North Ireland Environment Agency
}

\keyword{datasets}
