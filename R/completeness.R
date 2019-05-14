#' @title completeness
#' @description A function to plot biodiversity data.
#' @param input  df or raster
#' @param crs projection of data provided ('longlat'/'cea'/'auto')
#' @param zoom_out zoom level for the map
#' @examples
#' crop_map_world(df)
#' crop_map_world(df,crs='longlat')
#' @export
#'

completeness=function(gbif_df,raster_iucn,crs='longlat'){
  df=gbif_df
  df_coord=df[,3:4]
  if(crs=='cea'){
    spdf <- SpatialPointsDataFrame(df_coord, df,proj4string =
                                     CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

    r=raster(extent(raster_iucn), resolution=100000,crs=CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

  }
  else
  {
    spdf <- SpatialPointsDataFrame(df_coord, df,proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

    spdf=spTransform(spdf,CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

    r=raster(extent(raster_iucn), resolution=100000,crs=CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

  }

  print('rasterizing gbif')

  gbif_temp = rasterize(spdf,r,fun=function(x, ...) {length(unique(na.omit(x)))},field=spdf$species)


  #shape_iucn=gBuffer(shape_iucn,0.5,byid=T,id=shape_iucn$layer)

  #iucn_temp=rasterize(shape_iucn,r,fun=function(x, ...) {length(unique(na.omit(x)))},field=shape_iucn$layer)

  #iucn_temp_cea= spTransform(shape_iucn,CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  #r_cea=raster(ncol=84, nrow=77,extent(raster_iucn), resolution=100000,crs=CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  #print('rasterizing iucn cea')
  #iucn_temp_cea=rasterize(iucn_temp_cea,r_cea,fun=function(x, ...) {length(unique(na.omit(x)))},field=iucn_temp_cea$ID)
  iucn_temp=raster_iucn


  #gbif_temp = rasterize(spdf,r,fun=function(x, ...) {length(unique(na.omit(x)))},field=spdf$species)

  r_final=gbif_temp/raster_iucn
  plot(r_final)
  r_final_list=list(r_final,gbif_temp)
  #final_object=list(r_final,gbif_temp,iucn_temp,iucn_temp_cea)
  return(r_final_list)

}
