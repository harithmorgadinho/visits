#' @title clean_gbif_completeness
#' @description A function to plot biodiversity data.
#' @param input  df or raster
#' @param crs projection of data provided ('longlat'/'cea'/'auto')
#' @param zoom_out zoom level for the map
#' @examples
#' crop_map_world(df)
#' crop_map_world(df,crs='longlat')
#' @export
#'@importFrom CoordinateCleaner cc_cap
#'@importFrom CoordinateCleaner cc_cen
#'@importFrom CoordinateCleaner cc_inst
#'@importFrom CoordinateCleaner cc_dupl
#'@importFrom CoordinateCleaner cc_equ
#'
#'
clean_gbif_completeness=function(gbif_df,iucn_shp,taxa_list,method='zizka',crs=NULL){
  #only use names in gbif that exist in iucn_shp
  new_names_list=intersect(iucn_shp@data$BINOMIAL,gbif_df$species)
  gbif_df=gbif_df %>% filter(species%in%new_names_list)

  #subset iucn to taxa of interest
  iucn_shp=iucn_shp[iucn_shp$BINOMIAL %in% taxa_list, ]


  if(method=='zizka'){
    #clean species outside polygons method 1
    x=gbif_df
    colnames(x)=c('BINOMIAL',"year","decimallongitude","decimallatitude")
    range=iucn_shp
    gbif_df=cc_iucn(x, range, lon = "decimallongitude", lat = "decimallatitude",
                    species = "BINOMIAL", buffer = 0, value = "clean", verbose = TRUE)
  }


  if(method=='josue'){
    #clean species outside polygons method 2

    if(is.null(crs=='cea')){
      spdf <- SpatialPointsDataFrame(df_coord, df,proj4string =
                                       CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    }
    else
    {
      spdf <- SpatialPointsDataFrame(df_coord, df,proj4string =
                                       CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    }

    i=1
    spp_to_intersect=unique(spdf$species)
    gbif_new_shp_pts=spdf
    dat_2<- raster::intersect(gbif_new_shp_pts[gbif_new_shp_pts$species==spp_to_intersect[i],],buffer(iucn_shp[iucn_shp@data$BINOMIAL==spp_to_intersect[i],],0))
    #Set a progress bar
    pb <- txtProgressBar(min = 0, max = length(spp_to_intersect), style = 3)
    # Intersect each
    for (i in 2: length(spp_to_intersect)){
      setTxtProgressBar(pb, i)
      spp_now <- spp_to_intersect[i]
      dat_1<- raster::intersect(gbif_new_shp_pts[gbif_new_shp_pts$species==spp_now,],buffer(iucn_shp[iucn_shp@data$BINOMIAL==spp_now,],0))
      if (length(dat_1)==0){
        print(spp_now)
      } else{
        dat_2<- rbind(dat_2,dat_1)
      }
    }

    gbif_df=as.data.frame(dat_2)
    gbif_df=gbif_df[,c(1,2,5,6)]
    max(gbif_df$decimalLongitude.1)
  }
  #use ccleaner
  colnames(gbif_df)=c('species',"year","decimalLongitude","decimalLatitude")
  #gbif_df=cc_cap(gbif_df, lon = 'decimalLongitude', lat =  'decimalLatitude')
  gbif_df=cc_cen(gbif_df, lon = 'decimalLongitude', lat =  'decimalLatitude',test = "country")
  #gbif_df=cc_inst(gbif_df, lon = 'decimalLongitude', lat =  'decimalLatitude')
  #cc_dupl(input3, lon = 'decimalLongitude', lat =  'decimalLatitude')

  gbif_df=cc_equ(gbif_df, lon = 'decimalLongitude', lat =  'decimalLatitude')

  gbif_iucn=list(gbif_df,iucn_shp)
  return(gbif_iucn)

}
