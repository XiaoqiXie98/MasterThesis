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

mz2021=rast(here("mz2021"))#load the raster map with values representing management status in 2021. We will use this map to align the position of cells in all other maps. We want cell positions in different maps to be same
mz2021=disagg(mz2021,4)#The orignal cell size is 2000 m. We disaggregate it into 500 m by separate each cell to four cells.

#Load maps with detection of red-top tree (last year's infested tree) and green-dead tree (current year's infested tree).
#We did not consider cut and burn data because we only care about the influence of last year's infestations that will affect cuurent year.
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

#
nei3Red=list()#neighbor's values within 0.75 km. This is similar to look for the cells adjacent tto the target cell.
#Now modify the value in each cell in rastR and rastG by the sum of their adjacent cells.
for(i in 1:14){
  nei3Red[[i]]=focal(rastR_new[[i]],w=matrix(c(1,1,1,1,0,1,1,1,1),nrow=3,byrow=T),fun=sum,na.rm=T)
}
nei3Green=list()
for(i in 1:11){
  nei3Green[[i]]=focal(rastG_new[[i]],w=matrix(c(1,1,1,1,0,1,1,1,1),nrow=3,byrow=T),fun=sum,na.rm=T)
}
#make table from raster
nei3=as.numeric()
#Now, I look at 2006-2019
for(i in 1:14){
  e=as.matrix(extract(nei3Red_new[[i]],xy))#extract the sum of infestations from adjacent cells for all interested coordinates.
  e[is.na(e)]=0#set NA to be 0
  nei3=cbind(nei3,e)#One column represents one year's value. I add data for different years into table.
}

Green=as.numeric()#last year's green-dead trees
for(i in 1:3){
  e=as.matrix(extract(nei3Green_new[[i]],xy))#extract the sum of infestations from adjacent cells for all interested coordinates.
  e[is.na(e)]=0#set NA to be 0
  Green=cbind(Green,e)#One column represents one year's value. I add data for different years into table.
}
Green=cbind(0,Green,0)#add 2006, 2006 has no green trees, similar for 2010
for(i in 4:9){
  e=as.matrix(extract(nei3Green_new[[i]],xy))#extract the sum of infestations from adjacent cells for all interested coordinates.
  e[is.na(e)]=0#set NA to be 0
  Green=cbind(Green,e)#One column represents one year's value. I add data for different years into table.
}
Green=cbind(Green,0,0)#2017,2018, no detected green
e=as.matrix(extract(nei3Green_new[[10]],xy))#extract the sum of infestations from adjacent cells for all interested coordinates.
e[is.na(e)]=0#set NA to be 0
Green=cbind(Green,e)#One column represents one year's value. I add data for different years into table.
#
nei3=nei3+Green
save(nei3,file="Climate Change/last year neighborhood infestation 0.75km.RData")









