rm(list=ls()) # Clear the memory

setwd("~/Documents/Work/Sandbox/Miraine/") # Change on your machine!
filename <- "treedata_mix.csv"

# General graphical parameters
par(las=1) # Horizontal labels

# LOAD TREE DATA
# Tree data is stored in a csv file with 4 columns:
#   id: identity of the individual
#   tbirth: time at which the individual was infected/born
#   tdeath: time at which the individual cleared the infection/died
#   by: id of the individual who infected the focal individual
#       (NA for the individual with the earliest birth/infection date)
treedata <- read.csv(filename, header = TRUE)

nindiv <- length(treedata$id) # Number of individuals

# Compute vector of indices of the infectors
# (necessary if the IDs != indives on the table)
bywhomindex <- rep(0, nindiv) # initialize
for(i in 1:nindiv){
  if(is.na(treedata$by[i])) # if no infector (first individual)
    bywhomindex[i] <- NA 
  else 
    bywhomindex[i] <- which(treedata$id==treedata$by[i]) # index of the infector
}

#######################################################################
# CHECK VALIDITY OF THE TREE DATA
#######################################################################
# Do a series of tests on the data:

# - Check that no individual infected itself
#   (first individual has by=NA)
autoinf <- any(treedata$id == treedata$by, na.rm = TRUE)
if(autoinf){
  errlines <- which(treedata$id == treedata$by)
  errmsg <- paste("Check lines", errlines, ". Individuals infected themselves.")
  stop(errmsg)
}

# - Check that infector was alive when infected a focal individual
# Infector was born
infb <- treedata$tbirth > treedata[bywhomindex,]$tbirth
if(any(infb==FALSE, na.rm = TRUE)){
  errlines <- which(infb==FALSE)
  errmsg <- paste("Check lines ", errlines, ". individuals were infected by individuals that were not already infected.")
  stop(errmsg)
}
# Infector was not dead
infd <- treedata$tbirth < treedata[bywhomindex,]$tdeath
if(any(infd==FALSE, na.rm = TRUE)){
  errlines <- which(infd==FALSE)
  errmsg <- paste("Check lines ", errlines, ". individuals were infected by individuals that were already dead")
  stop(errmsg)
}  



#######################################################################
# SORT TREE FOR CONTOUR
#######################################################################
sort.tree <- function(treedata){
  # Initialize the vector of indices of sorted individuals
  sorted <- rep(0, nindiv)
  # We have to be careful not to mix individual IDs
  # and individual indices on the tree (line in table). 
  
  # First one is the index on the tree of the individual with the earliest birth date
  sorted[1] <- which.min(treedata$tbirth)
  # currentbranch is the ID of the branch we are currently considering
  currentbranch <- treedata$id[sorted[1]] 
  # Remove the sorted lines from the tree
  tmptree <- treedata[-sorted[1], ] 

  # Initialize vector of branches (need to keep track)
  branches <- currentbranch

  for (i in 2:nindiv){
    # First, check whether we have some infections left on the branch. 
    while(all(tmptree$by != currentbranch)){ 
      # If there are no infections left, we have to go back on previous branches
      for(j in seq(length(branches)-1, 1, by=-1)){
        currentbranch <- branches[j]
        # Stop on the first branch back where there are infections left
        if(any(tmptree$by == currentbranch)) break 
      }
    } # end of the while loop
    
    # There are infections left on the current branch
    # so find the latest remaining one
    
    # Temporary tree of just the infections from the current branch
    infecttree <- tmptree[tmptree$by == currentbranch, ] 
    # Find index of latest infection on this infecttree
    whichmaxoninfect <- which.max(infecttree$tbirth)
    # ID of the individual
    itmpmax <- infecttree$id[whichmaxoninfect]
    # Its position on the tmptree
    jtmp <- which(tmptree$id == itmpmax)

    # Update the current branch with the ID of the new individual
    currentbranch <- itmpmax # ID
    # Update the sorted vector with the index in the original treedata of the new individual
    sorted[i] <- which(treedata$id==itmpmax) 
    # Remove the individual from the list of remaining individuals
    tmptree <- tmptree[- jtmp, ] # jmpt is its index on tmptree
    # Keep a record of the branches
    branches <- c(branches, currentbranch)
  }
  # Return sorted tree
  return(treedata[sorted, ])
}

# Save the original tree
oldtreedata <- treedata
# Sort it for the contour
treedata <- sort.tree(treedata)
treedata

# Erase bywhomindex because it has changed since we sorted the tree
rm(bywhomindex) 

#######################################################################
# PLOT TREE
#######################################################################
# Function to initialize the plot window
initplot <- function(xmax = nindiv+1){
  plot(0, 0, xlim=c(0, xmax), ylim=c(0, max(treedata$tdeath)), #
       type="n", xlab = "", ylab = "Time", axes = FALSE)
  # Draw vertical axis
  axis(2)
}

# x Positions of the individuals
xvals <- 1:nindiv 

# Function to plot the tree
plot.tree <- function(treedata, xvals, coltree = gray(0.2)){
  # treedata: table of tree data values
  # xvals: x positions of the individuals (typically, 1:nindiv)
  # coltree: color of the tree 

  # Plot "lives" (segments from birth to death)
  segments(xvals, treedata$tbirth, xvals, treedata$tdeath, col = coltree, lwd = 3)
  
  # Compute vector of indices of the infectors
  # (necessary if the IDs != indives on the table)
  # (it's inside of the function so that the function works for different trees)
  bywhomindex <- rep(0, nindiv) # initialize
  for(i in 1:nindiv){
    if(is.na(treedata$by[i])) # if no infector (first individual)
      bywhomindex[i] <- NA 
    else 
      bywhomindex[i] <- which(treedata$id==treedata$by[i]) # index of the infector
  }
  
  # Plot infections
  segments(xvals, treedata$tbirth, xvals[bywhomindex], treedata$tbirth, col = coltree, lwd = 1, lty = 4)
}

# PLOTTING !

# Plot original tree
initplot()
plot.tree(oldtreedata, xvals)
title("original tree")

prt <- readline("press enter to see the tree sorted by birth dates")

# Sorting the tree by birth events
initplot()
isortbirth <- sort(oldtreedata$tbirth, index.return = TRUE)$ix
plot.tree(oldtreedata[isortbirth, ], xvals)
title("tree sorted by birth events")

prt <- readline("press enter to see the tree sorted for the contour")

initplot()
plot.tree(treedata, xvals)
title("sorted tree (for the contour)")

prt <- readline("press enter to continue and see the contour")



#######################################################################
# PLOT CONTOUR
#######################################################################

# Find x positions of the individuals 
# (they depend on the tree and are such that the contour lines all have the same slope)
xcont <- rep(0, nindiv+1) # initialize the vector
xcont[1] <- 0 # First individual at time t=0
for(i in 2:nindiv){
  xcont[i] <- xcont[i-1] + treedata$tdeath[i-1] - treedata$tbirth[i]
}
# then the contour goes to 0 after the last individual
xcont[nindiv+1] <- xcont[nindiv] + treedata$tdeath[nindiv]

# Function to plot the contour
plot.contour <- function(treedata, xcont, contcol = rgb(1, 0.5, 0)){
  # treedata: table of tree data
  # xcont: x positions of the individuals
  # contcol: color of the contour
  segments(xcont[1:nindiv], treedata$tdeath, xcont[2:(nindiv+1)], c(treedata$tbirth[2:nindiv], 0), col = contcol, lwd = 2)
  segments(xcont[1:nindiv], treedata$tbirth, xcont[1:nindiv], treedata$tdeath, col = contcol, lwd = 1, lty = 2)  
}

# PLOTTING !
initplot(xmax = max(30, xcont))
axis(1, pos=0)
title("Tree and its contour")
# Plot the tree, in grey
plot.tree(treedata, xcont[1:nindiv], gray(0.4))

prt <- readline("press enter to add the contour")

plot.contour(treedata, xcont)

