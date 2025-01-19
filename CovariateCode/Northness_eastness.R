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

#Northness and eastness are estimated from aspect
aspect=rast(here("aspectaep"))#load aspect map. This map uses a cell size of 100 m.
aspect=aggregate(aspect,fact=5,fum=mean)#modify the cell size to be 500 m by taking mean of the covered cells.
aspect=resample(aspect,mz2021,'near')#align cell positions to be consistent with map ma2021.
aspect=focal(aspect,matrix(1,5,5),fun=mean,na.policy="only")#For locals with NA, I compute them by the mean of neighborhoods.
a=as.matrix(extract(aspect,xy))#extract aspect values for interested cells
rad=a/180*pi#degree to radius
northness=as.numeric(replicate(14,cos(rad)))#Facing north is estimated by cos(radius). These values won't change over years. Thus, the table contains 14 columns representing 14 years with same value every year.
easterness=as.numeric(replicate(14,sin(rad)))#Facing east is estimated by sin(radius). These values won't change over years. Thus, the table contains 14 columns representing 14 years with same value every year.
save(northness,file="Climate Change/northness.RData")#save the computed values
save(easterness,file="Climate Change/easterness.RData")#save the computed values
