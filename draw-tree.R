setwd("~/Documents/Work/Sandbox/Miraine/")

# Load tree data
treedata <- read.csv("treedata.csv", header = TRUE)
treedata

# Check validity of the data
check.tree.data <- function(tab){
  # Check that no individual infected itself
  # (first individual has by=NA)
  autoinf <- any(tab$id == tab$by, na.rm = TRUE)
  if(autoinf) stop("Check infection data: some individuals infected themselves")
  
  # Check that infector was alive when infected a focal individual
  infb <- tab$tbirth > tab[tab$by,]$tbirth
  infd <- tab$tbirth < tab[tab$by,]$tdeath
  if(any(infb==FALSE, na.rm = TRUE)) stop("Check infection times: some individuals were infected by individuals that were not already infected")
  if(any(infd==FALSE, na.rm = TRUE)) stop("Check infection times: some individuals were infected by individuals that were already dead")
  
  return(TRUE)
}

check.tree.data(treedata)

nindiv <- length(treedata$id)

#######################################################################
# PLOT TREE
#######################################################################
# Initialize the plot window
plot(0, 0, xlim=c(0, nindiv+1), ylim=c(0, max(treedata$tdeath)), #
     type="n", xlab = "", ylab = "Time", axes = FALSE)
axis(2)

xvals <- 1:nindiv # x positions of the individuals

# Plot lives
segments(xvals, treedata$tbirth, xvals, treedata$tdeath, col = gray(0.2), lwd=2)

# Plot infections
segments(xvals, treedata$tbirth, treedata$by, treedata$tbirth, col = gray(0), lwd = 1, lty = 2)

# CONTOUR
# Find x positions of the individuals
xcont <- rep(0, nindiv)
xcont[1] <- 0
for(i in 2:nindiv){
  xcont[i] <- xcont[i-1] + treedata$tdeath[i-1] - treedata$tbirth[i]
}



