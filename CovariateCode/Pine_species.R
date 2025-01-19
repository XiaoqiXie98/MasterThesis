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

#
hybrid=rast(here("hybridaep1000"))#load map
hybrid=focal(hybrid,matrix(1,3,3),fun=mean,na.policy="only")#For those locations with NA, we replace the cell value by mean of the nearest 9 cells' values. 
hybrid=project(hybrid,method='near',crs(mz2021),res=500)#adjuat the coordinate system to be consistent with the map regrading management status in 2021.
hybrid=resample(hybrid,mz2021,'near')#align the cell position to be consistent with mz2021 map.
hybrid=disagg(hybrid,fact=2)#The orginal map has a cell value of 1000 m. I changed it into 500 m.
q=as.matrix(extract(hybrid,xy))#extract cell values for the interested coordinates xy.

save(q,file="Climate Change/Q.RData")
