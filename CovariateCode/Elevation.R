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

#elevation
elevation=rast(here('elevationaep'))#load elevation map
elevation=aggregate(elevation,fact=5,fun=mean)#The elevation map has a cell size of 100 m. I change the cell size to be 500 m by taking mean of the covered cell's values.
elevation=resample(elevation,mz2021,'near')#I align the cell position to be consistent with positions in map mz2021.
elevation=focal(elevation,matrix(1,5,5),fun=mean,na.policy="only")#For those cells with NA, I compute cell values by mean of neighborhoods.
e=as.matrix(extract(elevation,xy))#extract elevation values for interested xy.
length(which(is.na(e)))#if there is missing value
elevation=as.numeric(replicate(14,as.matrix(e)))#Elevation won't change over years. Thus, I created a table with 14 columns representing years with same values.
save(elevation,file="Climate Change/elevation.RData")#save elevation data
#