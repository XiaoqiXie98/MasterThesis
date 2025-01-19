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

#In this script, we estimate the number of lodgepole and jack pines per cell.

#total volume *(percentage abundance of lodgepole and jack pines)=volume of lodgepole and jack pines
total_volume=rast(here("NFI layers 2011/Structure_Volume_Total.tif"))#load the total tree volume map
total_volume=aggregate(total_volume,fact=2,method="mean")#This map uses a cell of size 250 m. I aggregate cells so that the cell size becomes 500 m. Aggregation method is average of four cells.
jack_percent=rast(here("NFI layers 2011/Species_Pinu_Ban.tif"))#load the map containing  percentage abundance of jack pine.
jack_percent=aggregate(jack_percent,fact=2,method="mean")#The cell size in this map is 250 m. We aggregate cells so that cell size becomes 500 m. Aggregation method is average of four cell.
lodge_percent=rast(here("NFI layers 2011/Species_Pinu_Con.tif"))#load the map containing percentage abundance of lodgepole pine.
lodge_percent=aggregate(lodge_percent,fact=2,method="mean")#The cell size in this map is 250 m. I aggregate cells so that the cell size becomes 500 m.

total_volume1=project(total_volume,crs(mz2021),res=500)#change the coordinate system of the total-volume map
total_volume2=resample(total_volume1,mz2021,method="bilinear")#set the cell position in this changed map to be consistent with the map of management zone in 2021.
total_volume2=crop(total_volume2,mz2021)#The area covered in the total_volume2 is larger than we needed. I crop the map to roughly cover Alberta.

#I did the similar operations for the percentage abundance lodgepole and jack pine maps.
jack_percent1=project(jack_percent,crs(mz2021),res=500)
jack_percent2=resample(jack_percent1,mz2021,method="bilinear")
jack_percent2=crop(jack_percent2,mz2021)

lodge_percent1=project(lodge_percent,crs(mz2021),res=500)
lodge_percent2=resample(lodge_percent1,mz2021,method="bilinear")
lodge_percent2=crop(lodge_percent2,mz2021)
##################################################################################
#After preparing the maps, we start computing.

xy=xyFromCell(total_volume2,1:ncell(total_volume2))#extract coordinates
totalVolume=extract(total_volume2,xy)#find cell value for each coordinate
jackPercent=extract(jack_percent2,xy)/100#transform percentage to decimals
lodgePercent=extract(lodge_percent2,xy)/100#transform percentage to decimals
#Notice that the total volume and percentage abundances have a unit value per hectare (value/100x100 m^2)
pine_percent_v=totalVolume*(jackPercent+lodgePercent)#estimate volume per hectare

#This is the function derived in Goodsman 2016 to transform volume to trees per hectare.
E=function(volume_hectare){
  a=17.6
  delta=0.00527
  f=a*volume_hectare*exp(-delta*volume_hectare)
  return(f)
  
}
###################################################################
#In this section, I computed trees per hectare value for each cell and made a raster map for these values.
expected_volume=E(pine_percent_v)
d=cbind(xy,expected_volume)
library(raster)
r=rasterFromXYZ(d)
crs(r)=crs(mz2021)
writeRaster(r, filename = "pine_density_numbers_hectare_beaudoin_11.tif", format = "GTiff",overwrite=TRUE)
#numbers per hectare
#############################################################
#The made map contains pine numbers per hectare in 2011. In the following part, I will adjust the numbers per cell by adding and subtracting the current year's infestations.

density=rast(here("pine_density_numbers_hectare_beaudoin_11.tif"))#load the made map
load(here("CurrentYearInfest.RData"))#load current year's infestations
#current pine density before attack
density2011=as.matrix(extract(density,xy)*500^2/100^2)#transform numbers per hectare to numbers per cell
density2011=round(density2011)#round values to integers
#Each column in pine density represents one year's pine trees per cell (I only considered lodgepole and jack pines).
pinedensity=as.numeric()
for(i in 1:3){#add 2007-2009 pines per cell
  #year before 2011 by adding infested trees
  v=apply(CurrentYearInfest[,c(i:4)],1,sum)
  pinedensity=cbind(pinedensity,density2011+v)
}
pinedensity=cbind(pinedensity,density2011+CurrentYearInfest[,4])#add 2010 pines per cell
pinedensity=cbind(pinedensity,density2011)#add 2011
v=CurrentYearInfest[,5]#2011 infestations
v1=density2011-v#2011 pines per cell - 2011 infestations
v1[which(v1<0)]=0 #cells with value less than 0 be come 0
pinedensity=cbind(pinedensity,v1)#add 2012

for(i in 6:13){#do for the rest of years
  v=apply(CurrentYearInfest[,c(5:i)],1,sum)#sum the value of rows from 2011 to the year we estimate for.
  v1= density2011-v#subtract those killed trees
  v1[which(v1<0)]=0#We set the cells with value less than 0 after subtraction being 0.
  pinedensity=cbind(pinedensity,v1)#add the counted year's value
}
colnames(pinedensity)=2007:2020#rename columns

save(pinedensity,file="Climate Change/pinedensity.RData")

