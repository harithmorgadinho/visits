#' @title visits_per_cell
#' @description A function to plot biodiversity data.
#' @param input  df or raster
#' @param crs projection of data provided ('longlat'/'cea'/'auto')
#' @param zoom_out zoom level for the map
#' @examples
#' crop_map_world(df)
#' crop_map_world(df,crs='longlat')
#' @export
#'

visits_per_cell=function(input,min_n_visits=5,extent_object){


  df_temp=input %>% group_by(cell_id) %>% summarise(a_p_sum = sum(a_p))
  cell_list_fixed=df_temp$cell_id[which(df_temp$a_p_sum>=min_n_visits)]
  input= input %>% filter(cell_id%in%cell_list_fixed)
  input=input[input$cell_id%in%cell_list_fixed,]
  df_temp=input %>% group_by(cell_id) %>% summarise(a_p_sum = sum(a_p))

  r=raster(ncol=84, nrow=77,extent(extent_object), resolution=100000,crs=CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

  df=as.data.frame(1:ncell(r))
  df$values=NA
  colnames(df)=c('cell_id_temp','values')
  colnames(df_temp)=c('cell_id_temp','values')
  merged=merge(df,df_temp,by='cell_id_temp',all=TRUE)
  values(r)=merged$values.y
  return(r)
}
