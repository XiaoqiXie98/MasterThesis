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

nam=list.files(path=here("NFI layers 2011"),pattern="Age",full.names=T)#use the pattern "Age" to find the directory of age map.
age=rast(nam)#load this map
age=aggregate(age,2,mean)#This map has a cell size 250 m. I aggregate this map to have a cell size of 500 m. The method for aggregation is mean.
age=project(age,crs(mz2021),res=500,method="near")#change the coordinate system to be similar to map mz2021.
age=resample(age,mz2021,'near')#align the cell position in the age map to be similar as in the map mz2021.
Age2011=round(as.matrix(extract(age,xy)),0)#only keep integer values.
length(which(is.na(Age2011)))#if there is missing values in this map
Age1=Age2011#create a table, which starts with one column. In the following loop, I will add up columns representing different years into this table.
#We only have one year average tree age. Thus, I estimate other years by subtracting and add 1.
Age=round(Age2011,0)#keep integer values
for(i in 1:4){#years before 2011
  e=Age2011-i#2007-2010
  e[which(e<0)]=0
  Age=cbind(e,Age)
}
for(i in 1:9){#years after 2011
  e=Age2011+i#2012-2020
  Age=cbind(Age,e)
}
colnames(Age)=2007:2020#rename column
save(Age,file="Climate Change/Age.RData")
