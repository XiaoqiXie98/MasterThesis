library(terra)
library(here)
getwd()#check my current directory
#find the interested locations
#We are interested in the locations with an active management status (which means the cell value is 1 in the 'map' map)
map=rast(here("study_area.tif"))#load the raster map 
plot(map)#see how the map looks like
xy=xyFromCell(map,1:ncell(map))#extract all xy coordinates from this map. Coordinates contain for cell value 0 and 1.
v=extract(map,xy)#extract the cell value using the extracted xy coordinates

xy=xy[which(v==1),]#only leave cells with a value 1
print(dim(xy)[1])#how many cells (rows) have a cell value of 1
#Overall, we have 622286 cells

#Now, we look at the covariate - current year's infestations (response variable).

mz2021=rast(here("mz2021"))#load the raster map with values representing management status in 2021. We will use this map to align the position of cells in all other maps. We want cell positions in different maps to be same
mz2021=disagg(mz2021,4)#The orignal cell size is 2000 m. We disaggregate it into 500 m by separate each cell to four cells.

#Load maps with detection of red-top tree (last year's infested tree) and green-dead tree (current year's infested tree).
#We have such data from 2007 to 2021 for red-top trees. Green-dead trees only observed in some years.
#I already separated maps of red-top trees and green-dead tree.
#In the following few lines, we use pattern="Red" and "Green" to call raster maps for red-top tree and green-dead tree.
#namR and namG give the directory of maps for red-top tree and green-dead tree.
namR=list.files(path=here("Infest Green Red"),pattern="Red",full.names=T)#red-top tree map
namR=namR[c(-1,-17)]#remove the unused directories
namG=list.files(path=here("Infest Green Red"),pattern="Green",full.names=T)#green-dead tree map
namG=namG[c(-12,-13)]#remove the unused directories
#
rastR=lapply(namR,rast)#load maps
rastG=lapply(namG,rast)#load maps
#align the cells in different maps by using management zone 2021 map

rastR_new=lapply(rastR,function(x) resample(x,mz2021,"near"))#resample cells by nearest neighborhood method
rastG_new=lapply(rastG,function(x) resample(x,mz2021,"near"))#resample cells by nearest neighborhood method

#estimate total number of green-dead trees for each year
GreenFromRed=as.numeric()
for(i in 2:15){#2008-2021 red-top trees=2007-2020 green-dead trees
  e=as.matrix(extract(rastR_new[[i]],xy))
  e[is.na(e)]=0
  GreenFromRed=cbind(GreenFromRed,e)
}
for(i in 2:15){#name columns to be 2007-2020 green-dead trees
  colnames(GreenFromRed)[i-1]=paste("GreenFromRed",2005+i)
}
#detected green
#Now I work on the green-dead tree's maps. We did not have green-dead trees every year. 
#In the dataset detectedGreen, each column represents one year's green-dead trees. I gradually added year's green-dead trees.
detectedGreen=as.numeric()
for(i in 1:3){#add 2007,2008,2009 green-dead trees
  e=as.matrix(extract(rastG_new[[i]],xy))#extract cell values
  e[is.na(e)]=0#For those cells with NA, we set NA to be 0
  detectedGreen=cbind(detectedGreen,e)#add a column
}
detectedGreen=cbind(detectedGreen,0)# 2010, no observed green-dead trees
for(i in 4:9){#add 2011, 2012, 2013, 2014, 2015, 2016 green-dead trees
  e=extract(rastG_new[[i]],xy)#extract cell values
  e[is.na(e)]=0#For those cells with NA, we set NA to be 0
  detectedGreen=cbind(detectedGreen,e)#add a column
}
detectedGreen=cbind(detectedGreen,0,0)#add 2017,2018, no detected green-dead trees
e=extract(rastG_new[[10]],xy)#2019, extract green-dead tree numbers from the map
e[is.na(e)]=0#set locations with NA to be 0
detectedGreen=cbind(detectedGreen,e)#add 2019, add 2019 green-dead trees into our data
detectedGreen=cbind(detectedGreen,0)#add 2020, no observed green-dead trees. I added 0.
for(i in 2:15){#add column names
  colnames(detectedGreen)[i-1]=paste("DetectedGreen",i+2005)
}

#count total number of green-dead trees from red-top tree maps and green-dead tree maps
greenStayToNextYear=detectedGreen+GreenFromRed#these green trees were not cut down through human operations

#adds controlled trees, which are green-dead trees cut down in the detected year so that they can not cause beetle emergence in the following year.
#But these trees still represent pines killed in the current year by the mountain pine beetle.
#Thus, the total number of current year's infestations are cut green+ green from red +detected green
namCtr=list.files(path=here("ctrol"),pattern="ctr",full.names=T)#use pattern "ctr" to find directories of all controlled trees.
rastCtr=lapply(namCtr,rast)#load raster maps
rastCtr_new=lapply(rastCtr, function(x) resample(x,mz2021,"near"))#resample by using the map management status 2021
#In total, we have the control maps from 2006-2020
#extract each year cut trees
#each column in data cutGreen represents one year's cut trees.
cutGreen=as.numeric()
for(i in 2:15){#2007-2020
  e=as.matrix(extract(rastCtr[[i]],xy))#extract value from cell center
  e[is.na(e)]=0#set cells with a value NA to be cells with a value 0
  cutGreen=cbind(cutGreen,e)#add column
}
for(i in 2:15){#add column names
  colnames(cutGreen)[i-1]=paste("cutGreen",i+2005)
}
#add cut trees 
CurrentYearInfest=greenStayToNextYear+cutGreen
save(CurrentYearInfest,file="CurrentYearInfest.RData")#save data
