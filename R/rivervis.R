# Version -------------------
# rivervis_v0.38.3
# R 3.0.1

# Path building -----------------------------------------------------------

PathBuild <- function(river, parent, OBN){
  p <- cbind(river, parent)
  for (i in 1:OBN) {
    p <- cbind(p,p[,2][match(p[,1+i],p[,1])])
    if ((sum(!is.na(p[,i+2])) == 0)) break
  }
  p <- p[,1:(ncol(p)-1)]
  p
}

# Calculate the sum without NA -------------------

SumNotNA <- function(x){
  sum(x,na.rm=TRUE)
}

# Calculation the location of river Mouth and river Source ----------------------------

MouthSource <- function(path, riverlayout, OBN){
  p <- path
  for (i in 1:OBN){
    p[p == riverlayout$river[i]] <- riverlayout$distance[i]
  }
  class(p) <- "numeric"
  m <- apply(p,1,SumNotNA)
  s <- m + riverlayout$length
  mouthsource <- cbind(rmouth = m, rsource = s)
  row.names(mouthsource) <- riverlayout$river
  mouthsource
}

# Relative position matrix -------------

RelPos <- function(path, riverlayout, OBN, DIGITMAX){
  p <- path
  for (i in 1:OBN){
    p[p == riverlayout$river[i]] <- riverlayout$position[i]
    p[i,1:sum(!is.na(p[i,]))] <- rev(p[i,][!is.na(p[i,])])
  }
  row.names(p) <- riverlayout$River
  colnames(p) <- paste("digit",c(0:DIGITMAX),sep="")
  p
}

# Digit weight-----------------

DigitWeight <- function(DIGITMAX){
  digitweight <- matrix(NA,(DIGITMAX-1),1)
  for (i in 1:(DIGITMAX-1)){
    digitweight[i,1] <- 10^(9-i) # Actually 4 is a much better choice than 10.
  }
  digitweight
}

# Relative position numeric matrix ---------------------  

RelPosMatrix <- function(pos, DIGITMAX){
  p <- pos[,2:(DIGITMAX+1)]
  #  for (i in 2:(DIGITMAX-1)){
  p[,1][is.na(p[,1])] <- 0
  p[,1][p[,1] == "L"] <- -1
  p[,1][p[,1] == "R"] <- 1
  p[,2:DIGITMAX][p[,2:DIGITMAX] == "L"] <- 1
  p[,2:DIGITMAX][is.na(p[,2:DIGITMAX])] <- 2
  p[,2:DIGITMAX][p[,2:DIGITMAX] == "R"] <- 3
  #  }
  class(p) <- "numeric"
  p
}

# Sorting ----------------

RowCal <- function(posmatrix, digitweight, riverlayout, path, OBN, DIGITMAX){
  s <- posmatrix[,2:DIGITMAX] %*% digitweight * posmatrix[,1]
  colnames(s) <- "s"
  OBNLEFT <- length(s[s<0])
  OBNRIGHT <- length(s[s>0])
  row <- c(c(1:length(s[s<0])),0,c(-1:-length(s[s>0])))
  k <- data.frame(MouthSource(path, riverlayout, OBN),s)
  k <- data.frame(k[order(k$s,k$rmouth),],row)
  
  for (i in 2:OBNLEFT){
    j = i
    while (all(k$rmouth[i] > k$rsource[k$row == (j-1)]) | 
             all(k$rsource[i] < k$rmouth[k$row == (j-1)])){
      k$row[i] <- j - 1
      j = j - 1
    }
  }
  
  for (i in (OBNLEFT+3):OBN){
    j = i
    while (all(k$rmouth[i] > k$rsource[k$row == (OBNLEFT+2-j)]) | 
             all(k$rsource[i] < k$rmouth[k$row == (OBNLEFT+2-j)])){
      k$row[i] <- OBNLEFT+2-j
      j = j - 1
    }
  }                    
  
  #  row <- matrix(k$row, dimnames = list(rownames(k),"Row"))
  row <- data.frame(river = rownames(k), row = -k$row)
  
  row
}



# RiverLayout Function ------------

RiverLayout <- function(river, length, parent, position, distance,
                        direction = 1, 
                        margin = 0.5){
  
  riverlayout <- data.frame(river = river, length = length, parent = parent, position = position, distance = distance, stringsAsFactors = FALSE)
  
  OBN <- nrow(riverlayout) # Observation number, or number or rows
  
  riverlayout <- cbind(rivercode = paste("river",c(1:OBN),sep=""), riverlayout, stringsAsFactors = FALSE) # Allocate rivercode for rivers
  
  riverlayout$parent[riverlayout$position == "M"] <- NA # Make sure the Parent of mainstream is NA
  
  path <- PathBuild(riverlayout$river, riverlayout$parent, OBN)
  
  DIGITMAX <- ncol(path)-1
  
  pos <- RelPos(path, riverlayout, OBN, DIGITMAX)
  
  digitweight <- DigitWeight(DIGITMAX)
  
  posmatrix <- RelPosMatrix(pos, DIGITMAX)
  
  # Calculate Mouth and Source
  
  riverlayout <- cbind(riverlayout, MouthSource(path, riverlayout, OBN))
  
  # Calculate Row
  
  riverlayout <- merge(riverlayout, RowCal(posmatrix, digitweight, riverlayout, path, OBN, DIGITMAX), by = "river", sort = FALSE)
  
  row <- riverlayout$row
  rsource <- riverlayout$rsource
  rmouth <- riverlayout$rmouth
  
  # Judge flow direction
  if (direction == -1){
    MAX.SOURCE <- max(rsource)
    row <- -row
    rsource <- MAX.SOURCE - rsource
    rmouth <- MAX.SOURCE - rmouth
    riverlayout$rsource <- rsource
    riverlayout$rmouth <- rmouth
    riverlayout$row <- row
  }
  
  # Calculate unit height
  H.MAX <- max(row)-min(row)+1 # total number of rows
  
  H.SIZE <- 1/(H.MAX + H.MAX*margin + 1) # define the unit height of each river. Assume the margin between rows is HSIZE/2
  
  Y.ZERO <- abs(min(row[row<=0])) * (margin*H.SIZE + H.SIZE) + margin*H.SIZE # get the y coordinate for row 0 as a reference line
  
  Y <- row * (margin*H.SIZE + H.SIZE) + Y.ZERO # Y of left bottom points of river rectangles
  
  # Calculate unit width
  W.MAX <- max(rsource, rmouth) - min(rsource, rmouth) # maximum width in original units (km)
  
  W.SIZE <- 1/W.MAX # leave some space for the right
  
  X1 <- rmouth * W.SIZE # X of Mouth location
  
  X2 <- rsource * W.SIZE # X of Source location
  
  list(riverdata = riverlayout, H.MAX = H.MAX, H.SIZE = H.SIZE, W.MAX = W.MAX, W.SIZE = W.SIZE, X1 = X1, X2 = X2, Y = Y, direction = direction)
  
}



# RiverDraw - draw rivers based on riverlayout =========================

RiverDraw <- function(riverlayout,
                      bd.col = "black",
                      ln.col = "grey40",
                      ln.lty = 3,
                      ln.lwd = 1,
                      bg.col = "grey80",
                      pt.shw = TRUE,
                      pt.col = "black",
                      pt.pch = 20,
                      pt.bg = "black",
                      pt.cex = 1,
                      pt.lwd = 1,
                      mar.t = 0.05,  
                      mar.b = 0.05,
                      mar.l = 0.2,
                      mar.r = 0.1
){
  
  # Data transfer
  RIVER.DATA <- riverlayout[[1]]
  H.MAX <- riverlayout[[2]]
  H.SIZE <- riverlayout[[3]]
  W.MAX <- riverlayout[[4]]
  W.SIZE <- riverlayout[[5]]
  X1 <- riverlayout[[6]]
  X2 <- riverlayout[[7]]
  Y <- riverlayout[[8]]
  DIRECTION <- riverlayout[[9]]
  
  
  par(mar=c(0,0,0,0))
  
  # Plot new sheet
  plot.new()
  
  # plotting margin
  par(usr = c(-mar.l, 1+mar.r, -mar.b, 1+mar.t))
  
  # Plot river rectangles
  rect(X1, Y, X2, Y+H.SIZE, col = bg.col, border = bd.col) # draw river rectangles
  
  # Plot lead line
  Y.PARENT <- Y[match(RIVER.DATA$parent, RIVER.DATA$river)] # Y of Parent of each river, using dictionary technique
  
  segments(X1, Y, X1, Y.PARENT, col = ln.col, lty = ln.lty, lwd = ln.lwd)
  
  # Plot river rectangle-frames
  rect(X1, Y, X2, Y+H.SIZE, col = NA, border = bd.col) # draw the frame of river rectangles, in case they have been covered by leadlines
  
  # Plot anchor points
  if (pt.shw){
    points(X1, Y.PARENT, pch=pt.pch, col = pt.col, bg = pt.bg, cex = pt.cex, lwd = pt.lwd) # plot anchor points
  }
  
}





# RiverMap - Draw river rectangles ===========================

RiverMap <- function(river, length, parent, position, distance, 
                     row = NA, 
                     direction = 1, # when Direction=1, river mouth on left. when Direction=-1, river mouth on right.
                     margin = 0.5,
                     bd.col = "black",
                     ln.col = "grey40",
                     ln.lty = 3,
                     ln.lwd = 1,
                     bg.col = "grey80",
                     pt.shw = TRUE,
                     pt.col = "black",
                     pt.pch = 20, 
                     pt.bg = "black",
                     pt.cex = 1,
                     pt.lwd = 1,
                     mar.t = 0.05,  
                     mar.b = 0.05,
                     mar.l = 0.2,
                     mar.r = 0.1
){ 
  
  par(mar=c(0,0,0,0))
  
  riverlayout <- data.frame(river = river, length = length, parent = parent, position = position, distance = distance, stringsAsFactors = FALSE)
  
  OBN <- nrow(riverlayout) # Observation number, or number or rows
  
  riverlayout <- cbind(rivercode = paste("river",c(1:OBN),sep=""), riverlayout, stringsAsFactors = FALSE) # Allocate rivercode for rivers
  
  riverlayout$parent[riverlayout$position == "M"] <- NA # Make sure the Parent of mainstream is NA
  
  path <- PathBuild(riverlayout$river, riverlayout$parent, OBN)
  
  DIGITMAX <- ncol(path)-1
  
  pos <- RelPos(path, riverlayout, OBN, DIGITMAX)
  
  digitweight <- DigitWeight(DIGITMAX)
  
  posmatrix <- RelPosMatrix(pos, DIGITMAX)
  
  # Calculate Mouth and Source
  
  riverlayout <- cbind(riverlayout, MouthSource(path, riverlayout, OBN))
  
  # Calculate Row
  
  
  if(all(is.na(row))){
    riverlayout <- merge(riverlayout, RowCal(posmatrix, digitweight, riverlayout, path, OBN, DIGITMAX), by = "river", sort = FALSE)
    
    row <- riverlayout$row
    
  } else{  
    riverlayout <- data.frame(riverlayout, row = row)
  }
  
  row <- riverlayout$row
  rsource <- riverlayout$rsource
  rmouth <- riverlayout$rmouth
  
  # Judge flow direction
  if (direction == -1){
    MAX.SOURCE <- max(rsource)
    row <- -row
    rsource <- MAX.SOURCE - rsource
    rmouth <- MAX.SOURCE - rmouth
    riverlayout$rsource <- rsource
    riverlayout$rmouth <- rmouth
    riverlayout$row <- row
  }
  
  # Calculate unit height
  H.MAX <- max(row)-min(row)+1 # total number of rows
  
  H.SIZE <- 1/(H.MAX + H.MAX*margin + 1) # define the unit height of each river. Assume the margin between rows is HSIZE/2
  
  Y.ZERO <- abs(min(row[row<=0])) * (margin*H.SIZE + H.SIZE) + margin*H.SIZE # get the y coordinate for row 0 as a reference line
  
  Y <- row * (margin*H.SIZE + H.SIZE) + Y.ZERO # Y of left bottom points of river rectangles
  
  # Calculate unit width
  W.MAX <- max(rsource, rmouth) - min(rsource, rmouth) # maximum width in original units (km)
  
  W.SIZE <- 1/W.MAX # leave some space for the right
  
  X1 <- rmouth * W.SIZE # X of Mouth location
  
  X2 <- rsource * W.SIZE # X of Source location
  
  # Plot new sheet
  plot.new()
  
  # plotting margin
  par(usr = c(-mar.l, 1+mar.r, -mar.b, 1+mar.t))
  
  # Plot river rectangles
  rect(X1, Y, X2, Y+H.SIZE, col = bg.col, border = bd.col) # draw river rectangles
  
  # Plot lead line
  Y.PARENT <- Y[match(riverlayout$parent, riverlayout$river)] # Y of Parent of each river, using dictionary technique
  
  segments(X1, Y, X1, Y.PARENT, col = ln.col, lty = ln.lty, lwd = ln.lwd)
  
  # Plot river rectangle-frames
  rect(X1, Y, X2, Y+H.SIZE, col = NA, border = bd.col) # draw the frame of river rectangles, in case they have been covered by leadlines
  
  # Plot anchor points
  if (pt.shw){
    points(X1, Y.PARENT, pch=pt.pch, bg = pt.bg, cex = pt.cex, lwd = pt.lwd, col = pt.col) # plot anchor points
  }
  
  list(riverdata = riverlayout, H.MAX = H.MAX, H.SIZE = H.SIZE, W.MAX = W.MAX, W.SIZE = W.SIZE, X1 = X1, X2 = X2, Y = Y, direction = direction)
  
}


# RiverFrame - draw boarders, leadlines and anchor points =======================

RiverFrame <- function(riverlayout,
                       ln.shw = T,
                       ln.col = "grey40",
                       ln.lty = 3,
                       ln.lwd = 1,
                       pt.shw = T,
                       pt.col = "black",
                       pt.pch = 20,
                       pt.bg = "black",
                       pt.cex = 1,
                       pt.lwd = 1,
                       bd.shw = T,
                       bd.col = "black"){
  # Data transfer
  RIVER.DATA <- riverlayout[[1]]
  H.MAX <- riverlayout[[2]]
  H.SIZE <- riverlayout[[3]]
  W.MAX <- riverlayout[[4]]
  W.SIZE <- riverlayout[[5]]
  X1 <- riverlayout[[6]]
  X2 <- riverlayout[[7]]
  Y <- riverlayout[[8]]
  DIRECTION <- riverlayout[[9]]
  
  # Plot frame
  if (bd.shw){
    rect(X1, Y, X2, Y+H.SIZE, border = bd.col) # draw river rectangles    
  }
  
  # Plot lead line
  Y.PARENT <- Y[match(RIVER.DATA$parent, RIVER.DATA$river)] # Y of Parent of each river, using dictionary technique
  
  if (ln.shw){
    segments(X1, Y, X1, Y.PARENT, col = ln.col, lty = ln.lty, lwd = ln.lwd)    
  }
  
  # Plot anchor points
  if (pt.shw){
    points(X1, Y.PARENT, pch=pt.pch, col = pt.col, bg = pt.bg, cex = pt.cex, lwd = pt.lwd) # plot anchor points
  }
}


# RiverPoints - draw lines and elevation profiles =======================

RiverPoint <- function(site, river, distance, value, riverlayout,
                       type = "l",
                       pt.col = "grey40",
                       pt.bg = "black",
                       pt.pch = 20,
                       pt.cex = 1,
                       lbl.cex = 0.7,
                       lbl.adj = c(0.5,2),
                       lbl.ofs = 0.5,
                       lbl.col = "black",
                       lbl.srt = 0,
                       lbl.pos = NULL,
                       lbl.shw = FALSE,
                       ln.lwd = 1){ 
  
  # Data transfer
  RIVER.DATA <- riverlayout[[1]]
  H.MAX <- riverlayout[[2]]
  H.SIZE <- riverlayout[[3]]
  W.MAX <- riverlayout[[4]]
  W.SIZE <- riverlayout[[5]]
  X1 <- riverlayout[[6]]
  X2 <- riverlayout[[7]]
  Y <- riverlayout[[8]]
  DIRECTION <- riverlayout[[9]]
  
  # Point plotting 
  VALUE.MAX <- max(value) # the largest value
  VALUE.MIN <- min(value) # the smallest value
  VALUE.SIZE <- H.SIZE * 0.9/(VALUE.MAX - VALUE.MIN) # a ratio to turn real elevation to plotting scale, assuming the largest value is 0.9*HSIZE
  
  # Direction converting
  if (DIRECTION == -1){
    length  <- RIVER.DATA$length[match(river, RIVER.DATA$river)]
    distance <- length - distance
    X.VALUE <- X2[match(river, RIVER.DATA$river)] + distance * W.SIZE # use dictionary technique to calculate the X of Elev points on plot
  }else{  
    # Calculate X and Y
    X.VALUE <- X1[match(river, RIVER.DATA$river)] + distance * W.SIZE # use dictionary technique to calculate the X of Elev points on plot
  }
  
  Y.VALUE <- Y[match(river, RIVER.DATA$river)] + value * VALUE.SIZE + H.SIZE * 0.05 # use dictionary technique to calculate the Y of Elev points on plot
  
  
  V <- data.frame(river=factor(river), X.VALUE, Y.VALUE)
  
  for (i in RIVER.DATA$river){
    points(V[which(river==i),]$X.VALUE, 
           V[which(river==i),]$Y.VALUE, type=type, col = pt.col, bg = pt.bg, pch = pt.pch, lwd = ln.lwd, cex = pt.cex)
  }
  
  if (lbl.shw){
    X.SITE <- X.VALUE
    Y.SITE <- Y[match(river, RIVER.DATA$river)] # Y of sampling sites
    text(X.SITE, Y.SITE, labels = site, cex = lbl.cex, adj = lbl.adj, srt = lbl.srt, offset = lbl.ofs, col = lbl.col, pos = lbl.pos)
  }
  
}



# Bar-charts - draw bar-charts ============================

RiverBar <- function(site, river, distance, value, riverlayout,
                     range = NA,
                     bar.w = 1,
                     bar.col = NA,
                     bd.col = "black",
                     lbl.cex = 0.7,
                     lbl.adj = c(0.5,2),
                     lbl.ofs = 0.5,
                     lbl.col = "black",
                     lbl.srt = 0,
                     lbl.pos = NULL,
                     lbl.shw = TRUE,
                     pt.shw = FALSE){ 
  
  # Data transfer
  RIVER.DATA <- riverlayout[[1]]
  H.MAX <- riverlayout[[2]]
  H.SIZE <- riverlayout[[3]]
  W.MAX <- riverlayout[[4]]
  W.SIZE <- riverlayout[[5]]
  X1 <- riverlayout[[6]]
  X2 <- riverlayout[[7]]
  Y <- riverlayout[[8]]
  DIRECTION <- riverlayout[[9]]
  
  
  # Direction converting
  if (DIRECTION == -1){
    length  <- RIVER.DATA$length[match(river, RIVER.DATA$river)]
    distance <- length - distance
    X.SITE <- X2[match(river, RIVER.DATA$river)] + distance * W.SIZE # use dictionary technique to calculate the X of sampling sites
  }else{  
    # Site location
    X.SITE <- X1[match(river, RIVER.DATA$river)] + distance * W.SIZE # use dictionary technique to calculate the X of sampling sites
  }
  
  Y.SITE <- Y[match(river, RIVER.DATA$river)] # Y of sampling sites
  
  
  # Site plot
  if (pt.shw){
    points(X.SITE, Y.SITE) # plot sampling sites
  }
  
  # Bar plotting
  if (all(is.na(range))){
    VALUE.MAX <- max(value)
    VALUE.MIN <- min(value)
  } else{
    VALUE.MAX <- max(range)
    VALUE.MIN <- min(range)
  }
  
  VALUE.SIZE <- H.SIZE * 0.9/(VALUE.MAX-VALUE.MIN)
  
  # Bar-charts
  W.BAR <- 0.01 * bar.w # bar width, default bar size is 0.01
  N.VALUE <- length(value) # Number of quantitative variables
  
  for (i in 1:N.VALUE){
    rect((X.SITE - W.BAR * N.VALUE/2 + (i-1)*W.BAR), 
         Y.SITE, 
         (X.SITE - W.BAR * N.VALUE/2 + i * W.BAR), 
         (Y.SITE + (value[,i] - VALUE.MIN) * VALUE.SIZE), 
         col = bar.col[i],
         border = bd.col)
  }
  
  if (lbl.shw){
    text(X.SITE, Y.SITE, labels = site, cex = lbl.cex, adj = lbl.adj, srt = lbl.srt, offset = lbl.ofs, col = lbl.col, pos = lbl.pos)
  }
} 



# Block-charts - draw block charts ====================================


RiverBlock <- function(site, river, distance, value, riverlayout, arrangement,
                       pt.shw = FALSE, # show=1, hide=0
                       hw.rat = 1.5,
                       h.gap = 0.05, # by default, H.GAP is H.SIZE * 0.05
                       w.gap = 0.025, # by default, W.GAP is W.SIZE * 0.025
                       block.col = NA,
                       block.lwd = 1,
                       bd.col = "grey20",                    
                       par.shw = TRUE,
                       par.pos = 2,
                       par.ofs = 1,
                       par.cex =0.6,
                       par.adj = c(1,0.5),
                       lbl.shw = TRUE,
                       lbl.cex = 0.7,
                       lbl.adj = c(0.5,2),
                       lbl.ofs = 0.5,
                       lbl.col = "black",
                       lbl.srt = 0,
                       lbl.pos = NULL){
  
  # Data transfer
  RIVER.DATA <- riverlayout[[1]]
  H.MAX <- riverlayout[[2]]
  H.SIZE <- riverlayout[[3]]
  W.MAX <- riverlayout[[4]]
  W.SIZE <- riverlayout[[5]]
  X1 <- riverlayout[[6]]
  X2 <- riverlayout[[7]]
  Y <- riverlayout[[8]]
  DIRECTION <- riverlayout[[9]]
  
  # Direction converting
  if (DIRECTION == -1){
    length  <- RIVER.DATA$length[match(river, RIVER.DATA$river)]
    distance <- length - distance
    X.SITE <- X2[match(river, RIVER.DATA$river)] + distance * W.SIZE # use dictionary technique to calculate the X of sampling sites
  }else{  
    # Site location
    X.SITE <- X1[match(river, RIVER.DATA$river)] + distance * W.SIZE # use dictionary technique to calculate the X of sampling sites
  }
  
  Y.SITE <- Y[match(river, RIVER.DATA$river)] # Y of sampling sites
  
  # Site plot
  if (pt.shw){
    points(X.SITE, Y.SITE) # plot sampling sites
  }
  
  # Plot block-charts
  N.SITE <- nrow(value) # site/observation number
  N.LINE <- length(arrangement)
  
  H.GAP <- H.SIZE * h.gap # "GAP" is the distance between blocks
  W.GAP <- H.SIZE * w.gap
  
  H.BLOCK <- (H.SIZE-H.GAP)/N.LINE - H.GAP
  # "b" is the height of the rectangle
  W.BLOCK <- H.BLOCK/hw.rat # "a" is the width of the rectangle
  
  Y.BLOCK <- rep(N.LINE:1, arrangement) # repeat the sequence VStr:1, repeating times is provided by VStr
  Y.BLOCK <- H.GAP + (Y.BLOCK-1)*(H.GAP+H.BLOCK)
  X.BLOCK <- sequence(arrangement)-1 # VStr provides the "to"s of the sequence
  
  N.PERLINE <- rep(arrangement, arrangement) # number of blocks per line
  
  for (i in 1:N.SITE){   # draw small blocks site by site
    
    rect((X.SITE[i]-W.BLOCK/2+X.BLOCK*(W.BLOCK+W.GAP)/N.PERLINE),
         (Y.SITE[i]+Y.BLOCK),
         (X.SITE[i]-W.BLOCK/2+X.BLOCK*(W.BLOCK+W.GAP)/N.PERLINE+(W.BLOCK-(N.PERLINE-1)*W.GAP)/N.PERLINE),
         (Y.SITE[i]+Y.BLOCK + H.BLOCK), 
         col = block.col[as.numeric(value[i,])], border = bd.col, lwd = block.lwd)
  }
  
  # Parameter names
  PAR.NAME.LIST <- split(colnames(value), rep(1:length(arrangement), arrangement))
  
  PAR.NAMES <- NA
  
  PAR.NAMES <- TitlePaste(PAR.NAME.LIST)
  
  X.PAR <- rep(min(X1[RIVER.DATA$row==0],X2[RIVER.DATA$row==0]), N.LINE)
  
  Y.PAR <- Y[RIVER.DATA$row==0] + sort(unique(Y.BLOCK), decreasing = TRUE) + H.BLOCK/2
  
  if (par.shw){
    text(X.PAR, Y.PAR, PAR.NAMES, pos = par.pos, offset = par.ofs, cex = par.cex, adj = par.adj)
  }
  
  if (lbl.shw){
    text(X.SITE, Y.SITE, labels = site, cex = lbl.cex, adj = lbl.adj, srt = lbl.srt, offset = lbl.ofs, col = lbl.col, pos = lbl.pos)
  }
  
  
}



# River names ==========================


RiverLabel <- function(riverlayout,
                       cex = 0.7,
                       adj = c(0, -1),
                       srt = 90,
                       col = "black",
                       pos = NULL,
                       offset = 0.5,
                       corner = "lb"     # left-top = lt, left-bottom = lb, right-top = rt, right-bottom = rb
){ 
  
  # Data transfer
  RIVER.DATA <- riverlayout[[1]]
  H.MAX <- riverlayout[[2]]
  H.SIZE <- riverlayout[[3]]
  W.MAX <- riverlayout[[4]]
  W.SIZE <- riverlayout[[5]]
  X1 <- riverlayout[[6]]
  X2 <- riverlayout[[7]]
  Y <- riverlayout[[8]]
  DIRECTION <- riverlayout[[9]]
  
  if (DIRECTION == -1){
    if (corner == "lb"){
      X.RIVER <- X2
      Y.RIVER <- Y
    }
    if (corner == "lt"){
      X.RIVER <- X2
      Y.RIVER <- Y + H.SIZE
    }
    if (corner == "rb"){
      X.RIVER <- X1
      Y.RIVER <- Y
    }
    if (corner == "rt"){
      X.RIVER <- X1
      Y.RIVER <- Y + H.SIZE
    }
    
  }else{
    if (corner == "lb"){
      X.RIVER <- X1
      Y.RIVER <- Y
    }
    if (corner == "lt"){
      X.RIVER <- X1
      Y.RIVER <- Y + H.SIZE
    }
    if (corner == "rb"){
      X.RIVER <- X2
      Y.RIVER <- Y
    }
    if (corner == "rt"){
      X.RIVER <- X2
      Y.RIVER <- Y + H.SIZE
    }
    
  }
  
  
  
  text(X.RIVER, Y.RIVER, labels = RIVER.DATA$river, cex = cex, adj = adj, srt = srt, col = col, pos = pos, offset = offset)
  
}




# Tick marks - drawing tick marks ========================


RiverTM <- function(tickmark, # a vector of tick mark values
                    value, # original data
                    riverlayout,
                    range = NA,
                    side = "L", # tick mark on left or right
                    pos = 1, # in (-1) or out (1)
                    tm.l = 1,
                    tm.col = "black",
                    lbl.shw = TRUE,
                    lbl.col = "black",
                    lbl.cex = 0.7,
                    lbl.row = TRUE,
                    label = NA){ # relative length of the tick mark
  
  # Data transfer
  RIVER.DATA <- riverlayout[[1]]
  H.MAX <- riverlayout[[2]]
  H.SIZE <- riverlayout[[3]]
  W.MAX <- riverlayout[[4]]
  W.SIZE <- riverlayout[[5]]
  X1 <- riverlayout[[6]]
  X2 <- riverlayout[[7]]
  Y <- riverlayout[[8]]
  DIRECTION <- riverlayout[[9]]
  
  #TickMark
  if (all(is.na(range))){
    tickmark <- tickmark[which(tickmark>(min(value)-(max(value)-min(value))*0.05/0.9) & tickmark<(max(value)+(max(value)-min(value))*0.05/0.9))]
  } else{
    tickmark <- tickmark[which(tickmark>(min(range)-(max(range)-min(range))*0.05/0.9) & tickmark<(max(range)+(max(range)-min(range))*0.05/0.9))]
  }
  
  
  L.TM <- 0.005 * tm.l # the default length of tick mark is 0.005
  
  # left or right?
  
  if (DIRECTION == -1){
    X <- X1
    X1 <- X2
    X2 <- X 
  }
  
  
  if(all(is.na(label))){
    
    label <- tickmark
  }
  
  
  if (all(is.na(range))){
    VALUE.MAX <- max(value)
    VALUE.MIN <- min(value)
  } else{
    VALUE.MAX <- max(range)
    VALUE.MIN <- min(range)
  }
  
  VALUE.SIZE <- H.SIZE * 0.9/(VALUE.MAX-VALUE.MIN)
  
  if (side == "L"){
    
    LOC.ROW <- data.frame(X.ROW = X1, Y.ROW = Y)
    
    if (lbl.row){
      for (i in (1:length(Y))[duplicated(Y)]){
        LOC.ROW <- LOC.ROW[LOC.ROW$X.ROW != max(LOC.ROW$X.ROW[LOC.ROW$Y.ROW==Y[i]]) | LOC.ROW$Y.ROW!=Y[i],]
      }
    }
    
    
    
    X.TM.LBL <- rep(LOC.ROW$X.ROW, each=length(tickmark))
    
    Y.TM.LBL <- rep(LOC.ROW$Y.ROW, each=length(tickmark)) + rep((tickmark - VALUE.MIN) * VALUE.SIZE, nrow(LOC.ROW)) + 0.05 * H.SIZE
    
    X.TM <- rep(X1, each=length(tickmark))
    
    Y.TM <- rep(Y, each=length(tickmark)) + rep((tickmark - VALUE.MIN) * VALUE.SIZE, nrow(RIVER.DATA)) + 0.05 * H.SIZE
    
    
    segments(X.TM, 
             Y.TM, 
             X.TM - pos * L.TM, 
             Y.TM, col = tm.col)
    
    if (lbl.shw){
      text(X.TM.LBL - L.TM,
           Y.TM.LBL,
           labels = label,
           cex = lbl.cex,
           adj = c(1,0.5),
           col = lbl.col)
    }
    
  }
  
  if (side == "R"){
    
    LOC.ROW <- data.frame(X.ROW = X2, Y.ROW = Y)
    
    if (lbl.row){
      for (i in (1:length(Y))[duplicated(Y)]){
        LOC.ROW <- LOC.ROW[LOC.ROW$X.ROW != min(LOC.ROW$X.ROW[LOC.ROW$Y.ROW==Y[i]])| LOC.ROW$Y.ROW!=Y[i],]
      }
    }
    
    X.TM.LBL <- rep(LOC.ROW$X.ROW, each=length(tickmark))
    
    Y.TM.LBL <- rep(LOC.ROW$Y.ROW, each=length(tickmark)) + rep((tickmark - VALUE.MIN) * VALUE.SIZE, nrow(LOC.ROW)) + 0.05 * H.SIZE
    
    X.TM <- rep(X2, each=length(tickmark))
    
    Y.TM <- rep(Y, each=length(tickmark)) + rep((tickmark - VALUE.MIN) * VALUE.SIZE, nrow(RIVER.DATA)) + 0.05 * H.SIZE
    
    segments(X.TM, 
             Y.TM, 
             X.TM + pos * L.TM, 
             Y.TM, col = tm.col)
    
    if (lbl.shw){
      text(X.TM.LBL + L.TM,
           Y.TM.LBL,
           labels = label,
           cex = lbl.cex,
           adj = c(0,0.5),
           col = lbl.col)
    }
    
  }
  
  
}

# RiverAxisLabel ====================

RiverAxisLabel <- function(label, riverlayout,
                           cex = 0.7,
                           adj = c(0.5, -2),
                           srt = 90,
                           col = "black",
                           pos = NULL,
                           offset = 0.5,
                           side = "L",
                           mainonly = TRUE){
  
  # Data transfer
  RIVER.DATA <- riverlayout[[1]]
  H.MAX <- riverlayout[[2]]
  H.SIZE <- riverlayout[[3]]
  W.MAX <- riverlayout[[4]]
  W.SIZE <- riverlayout[[5]]
  X1 <- riverlayout[[6]]
  X2 <- riverlayout[[7]]
  Y <- riverlayout[[8]]
  DIRECTION <- riverlayout[[9]]
  
  # left or right?
  
  if (DIRECTION == -1){
    X <- X1
    X1 <- X2
    X2 <- X 
  }
  
  
  # River Axis Title
  if (side == "L"){
    
    LOC.ROW <- data.frame(X.ROW = X1, Y.ROW = Y, ROW = RIVER.DATA$row)
    
    for (i in (1:length(Y))[duplicated(Y)]){
      LOC.ROW <- LOC.ROW[LOC.ROW$X.ROW != max(LOC.ROW$X.ROW[LOC.ROW$Y.ROW==Y[i]]) | LOC.ROW$Y.ROW!=Y[i],]
    }
    
    if (mainonly){
      text(LOC.ROW$X.ROW[LOC.ROW$ROW==0],
           LOC.ROW$Y.ROW[LOC.ROW$ROW==0] + 0.5 * H.SIZE,
           labels = label,
           cex = cex,
           adj = adj,
           col = col,
           srt = srt)
    }else{
      text(LOC.ROW$X.ROW,
           LOC.ROW$Y.ROW + 0.5 * H.SIZE,
           labels = label,
           cex = cex,
           adj = adj,
           col = col,
           srt = srt)
    }
    
  }
  
  if (side == "R"){
    
    LOC.ROW <- data.frame(X.ROW = X2, Y.ROW = Y, ROW = RIVER.DATA$row)
    
    
    for (i in (1:length(Y))[duplicated(Y)]){
      LOC.ROW <- LOC.ROW[LOC.ROW$X.ROW != max(LOC.ROW$X.ROW[LOC.ROW$Y.ROW==Y[i]]) | LOC.ROW$Y.ROW!=Y[i],]
    }
    
    if (mainonly){
      text(LOC.ROW$X.ROW[LOC.ROW$ROW==0],
           LOC.ROW$Y.ROW[LOC.ROW$ROW==0] + 0.5 * H.SIZE,
           labels = label,
           cex = cex,
           adj = adj,
           col = col,
           srt = srt)
    }else{
      text(LOC.ROW$X.ROW,
           LOC.ROW$Y.ROW + 0.5 * H.SIZE,
           labels = label,
           cex = cex,
           adj = adj,
           col = col,
           srt = srt)
    }
    
  }
}


# RiverSite - plot sites of interest=======================

RiverSite <- function(site, river, distance, group, riverlayout,
                      pt.pch = 21,
                      pt.col = NA,
                      pt.bg = "red",
                      pt.cex = 1,
                      lbl.cex = 0.5,
                      lbl.srt = 0,
                      lbl.adj = c(0.5,2),
                      lbl.col = "black",
                      lbl.pos = 1,
                      lbl.ofs = 0.5,
                      lbl.shw = TRUE){
  
  # Data transfer
  RIVER.DATA <- riverlayout[[1]]
  H.MAX <- riverlayout[[2]]
  H.SIZE <- riverlayout[[3]]
  W.MAX <- riverlayout[[4]]
  W.SIZE <- riverlayout[[5]]
  X1 <- riverlayout[[6]]
  X2 <- riverlayout[[7]]
  Y <- riverlayout[[8]]
  DIRECTION <- riverlayout[[9]]
  
  if (DIRECTION == -1){
    X.RIVER <- X2
    length  <- RIVER.DATA$length[match(river, RIVER.DATA$river)]
    distance <- length - distance
  }else{
    X.RIVER <- X1
  }
  
  
  # Location coordinates
  X.LOC <- X.RIVER[match(river, RIVER.DATA$river)] + distance * W.SIZE # X of locations
  
  Y.LOC <- Y[match(river, RIVER.DATA$river)] # Y of locations
  
  if(length(group)==1){
    points(X.LOC, Y.LOC, type = "p", pch = pt.pch[1], col = pt.col[1], bg = pt.bg[1], cex = pt.cex[1])
    if(lbl.shw){
      text(X.LOC, Y.LOC, labels = site, cex = lbl.cex, adj = lbl.adj, srt = lbl.srt, offset = lbl.ofs, col = lbl.col, pos = lbl.pos)
    }
  }else{
    points(X.LOC, Y.LOC, 
           type = "p", 
           pch = rep(pt.pch, length(group))[c(group)], 
           col = rep(pt.col,length(group))[c(group)], 
           bg = rep(pt.bg,length(group))[c(group)], 
           cex = rep(pt.cex, length(group))[c(group)])
    
    if (lbl.shw){
      text(X.LOC, Y.LOC, labels = site, 
           cex = rep(lbl.cex, length(group))[c(group)], 
           adj = rep(lbl.adj, length(group))[c(group)], 
           col = rep(lbl.col, length(group))[c(group)], 
           srt = lbl.srt,
           pos = lbl.pos,
           offset = lbl.ofs)
    }
  }
  
  
  
  
}




# RiverReach - draw river reaches =======================

RiverReach <- function(reach, river, from, to, group, style, riverlayout,
                       rea.pos = NA,       # absolute positions of lines
                       rea.col = "lightblue",
                       rea.lty = 1,
                       rea.lwd = 1,
                       rea.den = NULL,
                       bd.col = "black",
                       ln.shw = T,
                       ln.col = "grey40",
                       ln.lty = 3,
                       ln.lwd = 1,
                       pt.shw = T,
                       pt.col = "black",
                       pt.pch = 20,
                       pt.bg = "black",
                       pt.cex = 1,
                       pt.lwd = 1){
  
  # Data transfer
  RIVER.DATA <- riverlayout[[1]]
  H.MAX <- riverlayout[[2]]
  H.SIZE <- riverlayout[[3]]
  W.MAX <- riverlayout[[4]]
  W.SIZE <- riverlayout[[5]]
  X1 <- riverlayout[[6]]
  X2 <- riverlayout[[7]]
  Y <- riverlayout[[8]]
  DIRECTION <- riverlayout[[9]]
  
  # Lines or Bands
  Style.Line <- style
  Style.Line[Style.Line == 99] <- NA
  
  Style.Band <- style
  Style.Band[Style.Band != 99] <- NA
  
  if (length(style[style<99]) > 0){
    N.LINE <- max(style[style<99])
  }else{
    N.LINE <- 0
  }
  N.OBS <- nrow(style)
  
  
  # Direction converting
  
  if (DIRECTION == -1){
    length  <- RIVER.DATA$length[match(river, RIVER.DATA$river)]
    from <- length - from
    to <- length - to
    X.FROM <- X2[match(river, RIVER.DATA$river)] + from * W.SIZE # X of reach
    X.TO <- X2[match(river, RIVER.DATA$river)] + to * W.SIZE # X of reach
  }else{  
    # Calculate X and Y
    X.FROM <- X1[match(river, RIVER.DATA$river)] + from * W.SIZE # X of reach
    X.TO <- X1[match(river, RIVER.DATA$river)] + to * W.SIZE # X of reach
    
  }
  
  Y.REACH.BASE <- Y[match(river, RIVER.DATA$river)] # Y base of locations
  
  Y.REACH.BAND <- Y.REACH.BASE + Style.Band-99
  
  # draw bands and lines
  
  rect(X.FROM, Y.REACH.BAND, X.TO, (Y.REACH.BAND + H.SIZE),
       border = bd.col, density = rea.den,
       col = rep(rea.col, nlevels(group))[c(group)]
  )
  
  if (all(is.na(rea.pos))){
    
    Y.REACH.LINE <- Y.REACH.BASE + H.SIZE/(N.LINE+1)*Style.Line
    
    
  }else{
    
    rea.pos[which(style == 99)] <- NA
    
    Y.REACH.LINE <- Y.REACH.BASE + rea.pos * H.SIZE
    
  }
  
  segments(X.FROM, Y.REACH.LINE, X.TO, Y.REACH.LINE, 
           col = rep(rea.col, nlevels(group))[c(group)],
           lty = rep(rea.lty, nlevels(group))[c(group)], 
           lwd = rep(rea.lwd,nlevels(group))[c(group)])
  
  # Plot lead line
  Y.PARENT <- Y[match(RIVER.DATA$parent, RIVER.DATA$river)] # Y of Parent of each river, using dictionary technique
  
  if (ln.shw){
    segments(X1, Y, X1, Y.PARENT, col = ln.col, lty = ln.lty, lwd = ln.lwd)    
  }
  
  # Plot anchor points
  if (pt.shw){
    points(X1, Y.PARENT, pch=pt.pch, col = pt.col, bg = pt.bg, cex = pt.cex, lwd = pt.lwd) # plot anchor points
  }
  
}


# TitlePaste =========================

TitlePaste <- function(x){
  
  if(length(x)==1){
    x[[1]]
  }else{
    y <- NA
    
    for(i in 1:length(x)){
      y[i] <- paste(x[[i]], collapse = "/")
      
    }
    
    y
  }
  
  
}




# River Block Chart ====================================


RiverBlockChart <- function(site, river, distance, value, arrangement,
                            h.gap = 0.1, # by default, H.GAP is H.BLOCK * 0.5
                            w.gap = 0.25, # by default, W.GAP is W.BLOCk * 0.25
                            w.gap.s = 0.1, # by default, gap between small blocks is W.BLOCK *0.1
                            r.gap = 0.25, # by default, gap between rivers is W.BLOCK * 0.25
                            block.col = NA,
                            block.lwd = 1,
                            border.col = "grey20",
                            bg.col = "lightgrey",
                            mar = 0.1, # smallest margin
                            hw.rat = 1.5,
                            site.shw = TRUE,
                            site.pos = 1,
                            site.ofs = 1.5,
                            site.cex = 0.5,
                            site.col = "black",
                            site.order = "A", # alphabetical order ("A"), river flow left ("L"), river flow right ("R")
                            rvr.shw = TRUE,
                            rvr.ofs = 1.5,
                            rvr.cex = 0.7,
                            rvr.col = "black",
                            rvr.t.b = "b",
                            rvr.order = NA, # alphabetical order (NA) or custom order (a vector)
                            par.shw = TRUE,
                            par.pos = 2,
                            par.ofs = 1,
                            par.cex =0.6,
                            par.adj = c(1,0.5),
                            par.col = "black"){ 
  
  # plot new
  par(mar=c(0,0,0,0))
  plot.new()
  
  # layout
  N.SITE <- nrow(value) # site/observation number
  N.LINE <- length(arrangement)
  N.RIVER <- length(unique(river))
  
  if((N.SITE+N.RIVER-1)>N.LINE){
    MARGIN.LEFT <- mar
    MARGIN.RIGHT <- 1- mar
    WIDTH <- 1-mar*2
    W.BLOCK <- WIDTH/((1+w.gap)*N.SITE + N.RIVER*w.gap + (N.RIVER-1)*r.gap)
    H.BLOCK <- W.BLOCK * hw.rat
    HEIGHT <- H.BLOCK * (N.LINE + h.gap * N.LINE +h.gap)
    MARGIN.TOP <- (1-HEIGHT)/2
    MARGIN.BOTTOM <- MARGIN.TOP
  }else{
    MARGIN.TOP <- 1- mar
    MARGIN.BOTTOM <- mar
    HEIGHT <- 1-mar*2
    H.BLOCK <- HEIGHT/(N.LINE + h.gap * N.LINE +h.gap)
    W.BLOCK <- H.BLOCK/hw.rat
    WIDTH <- (W.BLOCK + w.gap*W.BLOCK) * N.SITE + N.RIVER*w.gap*W.BLOCK + (N.RIVER-1)*r.gap*W.BLOCK
    MARGIN.LEFT <- (1-WIDTH)/2
    MARGIN.RIGHT <- 1- MARGIN.LEFT
    
  }
  
  H.GAP <- H.BLOCK * h.gap # "GAP" is the distance between blocks
  W.GAP <- W.BLOCK * w.gap
  W.GAP.S <- W.BLOCK * w.gap.s
  R.GAP <- W.BLOCK * r.gap
  
  
  if(all(is.na(rvr.order))){
    rvr.order <- sort(unique(river))
  }
  
  # Site coordinates
  
  LOC.RIVER.CAL <- data.frame(river = sort(unique(river)), number = matrix(table(river)))
  
  LOC.RIVER.CAL$river <- factor(LOC.RIVER.CAL$river, levels = rvr.order)
  
  LOC.RIVER.CAL <- LOC.RIVER.CAL[order(LOC.RIVER.CAL$river),]
  
  for(i in N.RIVER:1){
    LOC.RIVER.CAL$accnumber[i] <- sum(LOC.RIVER.CAL$number[1:i-1])
  }
  
  LOC.RIVER.CAL$gapnum <- 0:(N.RIVER-1)
  LOC.RIVER.CAL$x <- LOC.RIVER.CAL$gapnum*(W.GAP + R.GAP) + LOC.RIVER.CAL$accnumber * (W.BLOCK + W.GAP)
  
  A <- data.frame(site, river, distance, NumberinGroup = rep(0,N.SITE), value)
  A$river <- factor(A$river, levels = rvr.order)
  
  if (site.order == "A"){
    A <- A[order(A$river, A$site),]
  }else if (site.order == "L"){
    A <- A[order(A$river, A$distance),]
  }else if (site.order == "R"){
    A <- A[order(A$river, -A$distance),]
  }else{
    
  }
  
  for (i in river){
    k = 1
    for (j in 1:N.SITE){
      if (A$river[j]==i){
        A[j, "NumberinGroup"] <- k
        k <- k + 1
      }
    }
  }
  
  X.RIVER <- LOC.RIVER.CAL$x[match(A$river,LOC.RIVER.CAL$river)]
  
  NumberinGroup <- A$NumberinGroup
  
  
  X.SITE <-  MARGIN.LEFT + X.RIVER + (NumberinGroup-1)*(W.BLOCK + W.GAP) + W.GAP
  
  Y.SITE <- rep((H.GAP + MARGIN.BOTTOM), N.SITE) # Y of sampling sites
  
  
  # Background
  
  X1.RIVER.BG <- MARGIN.LEFT + LOC.RIVER.CAL$x
  Y1.RIVER.BG <- MARGIN.BOTTOM + rep(0,N.RIVER) - 3 * H.GAP
  X2.RIVER.BG <- MARGIN.LEFT + LOC.RIVER.CAL$x + LOC.RIVER.CAL$number*(W.BLOCK+W.GAP)+W.GAP
  Y2.RIVER.BG <- MARGIN.BOTTOM +  rep(HEIGHT,N.RIVER) + 3 * H.GAP
  
  rect(X1.RIVER.BG, Y1.RIVER.BG, X2.RIVER.BG, Y2.RIVER.BG, border = NA, col = bg.col)
  
  
  # Block coordinates
  Y.BLOCK <- rep(N.LINE:1, arrangement) # repeat the sequence VStr:1, repeating times is provided by VStr
  Y.BLOCK <- (Y.BLOCK-1)*(H.GAP+H.BLOCK)
  X.BLOCK <- sequence(arrangement)-1 # VStr provides the "to"s of the sequence
  
  N.PERLINE <- rep(arrangement, arrangement) # number of blocks per line
  
  VALUE.BLOCK <- A[(ncol(A)-ncol(value)+1):ncol(A)]
  
  for (i in 1:N.SITE){   # draw small blocks site by site
    
    rect((X.SITE[i]+X.BLOCK*(W.BLOCK+W.GAP.S)/N.PERLINE),
         (Y.SITE[i]+Y.BLOCK),
         (X.SITE[i]+X.BLOCK*(W.BLOCK+W.GAP.S)/N.PERLINE+(W.BLOCK-(N.PERLINE-1)*W.GAP.S)/N.PERLINE),
         (Y.SITE[i]+Y.BLOCK + H.BLOCK), 
         col = block.col[as.numeric(VALUE.BLOCK[i,])], border = border.col, lwd = block.lwd)
    
  }
  
  # Site name
  if (site.shw){
    text(X.SITE+W.BLOCK/2, Y.SITE, A$site, pos= site.pos, offset = site.ofs, cex = site.cex, col = site.col)
  }
  
  # River name
  if (rvr.shw){
    if(!all(is.na(rvr.order))){
      RIVER.NAME <- rvr.order
    } else{
      RIVER.NAME <- sort(unique(river))
    }
    
    if (rvr.t.b == "b"){
      text((X1.RIVER.BG + X2.RIVER.BG)/2, Y1.RIVER.BG, RIVER.NAME, pos = 1, offset = rvr.ofs, cex = rvr.cex, col = rvr.col)
    } 
    if (rvr.t.b == "t"){
      text((X1.RIVER.BG + X2.RIVER.BG)/2, Y2.RIVER.BG, RIVER.NAME, pos = 3, offset = rvr.ofs, cex = rvr.cex, col = rvr.col)
    }
    
  }
  
  # Parameter names
  PAR.NAME.LIST <- split(colnames(value), rep(1:length(arrangement), arrangement))
  
  PAR.NAMES <- TitlePaste(PAR.NAME.LIST)
  
  X.PAR <- rep(min(X.SITE), N.LINE)
  
  Y.PAR <- sort(unique(min(Y.SITE)+Y.BLOCK + H.BLOCK/2), decreasing = TRUE)
  
  if (par.shw){
    text(X.PAR, Y.PAR, PAR.NAMES, pos = par.pos, offset = par.ofs, cex = par.cex, adj = par.adj, col = par.col)
  }
  
  
}


# RiverDirection ==================

RiverDirection <- function(riverlayout, 
                           loc = NA, # or loc = c(mar.l, mar.b)
                           arw.length = 0.05,
                           arw.lty = 1,
                           arw.lwd = 1,
                           arw.angle = 30,
                           arw.col = "black",
                           label = "Flow direction",
                           lbl.cex = 0.5,
                           lbl.pos = 4,
                           lbl.ofs = 0.5){
  
  # Data transfer
  
  DIRECTION <- riverlayout[[9]]
  
  # locator
  
  if (all(is.na(loc))){
    loc <- locator()
  }
  
  arrows(loc[[1]], loc[[2]], loc[[1]] + arw.length, loc[[2]], 
         length = arw.length * 2, code = 1.5-DIRECTION/2, lty = arw.lty, lwd = arw.lwd, angle = arw.angle, col = arw.col)
  
  text(max(loc[[1]],loc[[1]] + arw.length), loc[[2]], label,
       adj = 0, cex = lbl.cex, pos = lbl.pos, offset = lbl.ofs)
}

# RiverScale ======================================

RiverScale <- function(length, label, riverlayout,
                       loc = NA, # or loc = c(mar.l, mar.b)
                       scl.col = "black",
                       scl.lwd = 1,
                       lbl.cex = 0.5,
                       lbl.pos = 4,
                       lbl.ofs = 0.5){
  
  # Data transfer
  RIVER.DATA <- riverlayout[[1]]
  H.MAX <- riverlayout[[2]]
  H.SIZE <- riverlayout[[3]]
  W.MAX <- riverlayout[[4]]
  W.SIZE <- riverlayout[[5]]
  X1 <- riverlayout[[6]]
  X2 <- riverlayout[[7]]
  Y <- riverlayout[[8]]
  DIRECTION <- riverlayout[[9]]
  
  # locator
  
  if (all(is.na(loc))){
    loc <- locator()
  }
  
  segments(loc[[1]],loc[[2]], loc[[1]]+length*W.SIZE, loc[[2]], 
           col = scl.col, lty = 1, lwd = scl.lwd)
  segments(loc[[1]], loc[[2]], loc[[1]], loc[[2]] + length*W.SIZE*0.1, 
           col = scl.col, lty = 1, lwd = scl.lwd)
  segments(loc[[1]]+length*W.SIZE, loc[[2]], loc[[1]]+length*W.SIZE, loc[[2]] + length*W.SIZE*0.1, 
           col = scl.col, lty = 1, lwd = scl.lwd)
  
  text(max(loc[[1]],loc[[1]] + length*W.SIZE), loc[[2]], label,
       adj = 0, cex = lbl.cex, pos = lbl.pos, offset = lbl.ofs)
  
}