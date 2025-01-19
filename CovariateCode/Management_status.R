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

loc=xy#change name of the interested cells
nam=list.files(path=here(),pattern="mz",full.names=T)#look for the directories of maps with common pattern mz
nam=nam[c(1:14)]#We only make use of the first 14 names.
rastm=lapply(nam,rast)#load maps
for(i in 1:14){#align the cell positions to be consistent with positions in map mz2021
  rastm[[i]]=resample(rastm[[i]],mz2021,"near")
}

#The number used to represent management status can be different for each years. Thus, in the following, I modify the numbers used in different years to be consistent.
#I set 0 means inactive management and 1 means active management.
i=1#2007
l=extract(rastm[[i]],loc)#extract values for interested cells
l[is.na(l),1]=0#For cells with NA, I replace NA by 0.
mz=l#The first column values
for(i in 2:14){#add the rest of columns (2008-2020)
  l=extract(rastm[[i]],loc)
  l[is.na(l),1]=0
  mz=cbind(mz,l)
}

#I set leading 1, all others 0
#For different years, the number represent active and inactive are different. In the following code, I modify them so that 1 means active and 0 means inactive.
l1<-which(mz[,1]==2)
l2<-which(mz[,1]!=2)
mz[l1,1]<-1
mz[l2,1]<-0
#####
######
l1<-which(mz[,2]==2)
l2<-which(mz[,2]!=2)
mz[l1,2]<-1
mz[l2,2]<-0
###09 1 leading
mz[,3][which(mz[,3]!=1)]=0

l1<-which(mz[,4]==2)
l2<-which(mz[,4]!=2)
mz[l1,4]<-1
mz[l2,4]<-0

l1<-which(mz[,5]==2)
l2<-which(mz[,5]!=2)
mz[l1,5]<-1
mz[l2,5]<-0

l1<-which(mz[,6]==2)
l2<-which(mz[,6]!=2)
mz[l1,6]<-1
mz[l2,6]<-0

l1<-which(mz[,7]==2)
l2<-which(mz[,7]!=2)
mz[l1,7]<-1
mz[l2,7]<-0

l1<-which(mz[,8]==2)
l2<-which(mz[,8]!=2)
mz[l1,8]<-1
mz[l2,8]<-0

#15,16,17,18, 1 leading
for(i in 9:12){
  mz[,i][which(mz[,i]!=1)]=0
}
###
l1<-which(mz[,13]==2)
l2<-which(mz[,13]!=2)
mz[l1,13]<-1
mz[l2,13]<-0


l1<-which(mz[,14]==2)
l2<-which(mz[,14]!=2)
mz[l1,14]<-1
mz[l2,14]<-0
#############################################################################################
#After correcting cell values to use the same rule, I save this data.
save(mz,file="Climate Change/mz.RData")
