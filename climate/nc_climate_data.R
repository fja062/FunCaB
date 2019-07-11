# load packages
library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)


# read in file
ncpath <- "~/OneDrive - University of Bergen/Research/FunCaB/Data/climate_data/"
ncname <- "clt_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_mon_197001-197012_box"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "clt" # cloud_area_fraction

ncin <- nc_open(ncfname)
print(ncin)

# get longitude and latitude
lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon)
lat <- ncvar_get(ncin,"lat")
nlat <- dim(lat)

# get time
time <- ncvar_get(ncin,"time")
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)

# get temperature variable and its attributes
tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)

# convert time -- split the time units string into fields
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth <- as.integer(unlist(tdstr)[2])
tday <- as.integer(unlist(tdstr)[3])
tyear <- as.integer(unlist(tdstr)[1])
chron(time,origin=c(tmonth, tday, tyear))

# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA

# get a single slice or layer (January)
m <- 1
tmp_slice <- tmp_array[,,m]

# quick map
ggplot(aes(x = lon, y = lat)) +
  geom_map()

image(x = lon, y = lat, z = as.matrix(tmp_slice), col=rev(brewer.pal(10,"RdBu")))


ggplot(grid, aes(x= lon, y = lat)) +
  geom_raster()

grid <- expand.grid(lon=lon, lat=lat)
cutpts <- c(0,10,20,30,40,50, 60, 70, 80,90)
levelplot(tmp_slice ~ lon * lat, data=grid, pretty=T, 
          col.regions=(rev(brewer.pal(10,"RdBu"))))
