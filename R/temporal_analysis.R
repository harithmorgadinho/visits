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
temporal_analysis=function(input,raster_input){
  #df = species / date of record / long/ lat
  print('step 1')

  df=input
  df_coord=df[,3:4]
  spdf <- SpatialPointsDataFrame(df_coord, df,proj4string =
                                   CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

  year_list=unique(input$year)
  year_list <- sort(year_list, decreasing = FALSE, na.last = NA)

  r=raster(ncol=84, nrow=77,extent(shp_behrmann), resolution=100000,crs=CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

  r2=r
  values(r2) <- seq(1:length(r2))
  rr <- rasterToPolygons(r2)
  rr_sf=st_as_sf(rr)


  spdf_sf=st_as_sf(spdf)
  spdf_sf_st <- st_transform(spdf_sf, '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')
  list_cells_df=st_intersection(spdf_sf_st,rr_sf)

  list_cells_df=unique(list_cells_df$layer)

  ##calculated in script loop_find_average

  print('step 2')
  #start here
  df_big=matrix(0, nrow = length(unique(input$species)), ncol = length(list_cells_df))
  rownames(df_big)=unique(input$species)
  colnames(df_big)=list_cells_df


  #for (i in seq_along(year_list)){


  #year=year_list[[i]]

  # subset_df_by_year=spdf[spdf$year==year,]
  # subset_df_by_year=st_as_sf(subset_df_by_year)

  # subset_df_by_year <- st_transform(subset_df_by_year, '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')

  #subset_df_by_year_inters=st_intersection(subset_df_by_year,rr_sf)
  # subset_df_by_year_inters_table<-table(subset_df_by_year_inters$layer,subset_df_by_year_inters$species)

  #df_small=t(subset_df_by_year_inters_table)

  # df_big[rownames(df_small), colnames(df_small)] <-
  # df_big[rownames(df_small), colnames(df_small)] +
  # df_small

  #}

  print('step 3')

  DF_BY_YEAR = function(year_list) {

    year=year_list

    subset_df_by_year=spdf[spdf$year==year,]
    subset_df_by_year=st_as_sf(subset_df_by_year)

    subset_df_by_year <- st_transform(subset_df_by_year, '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')

    subset_df_by_year_inters=st_intersection(subset_df_by_year,rr_sf)
    subset_df_by_year_inters_table<-table(subset_df_by_year_inters$layer,subset_df_by_year_inters$species)

    df_small=t(subset_df_by_year_inters_table)

    df_big[rownames(df_small), colnames(df_small)] <-
      df_big[rownames(df_small), colnames(df_small)] +
      df_small



    return(df_big)
  }

  list_of_big_df=pblapply(year_list,DF_BY_YEAR)
  names(list_of_big_df)=year_list



  # loop each year

  print('step 4')

  k=1

  vec_final=vector("list",length = length(year_list))

  df_final=list_of_big_df[[year_list[k]]]

  #loop calculates ratio in cells after every year

  for (k in 1:length(year_list)){
    cat('.... .... .... year [', year_list[k], '/', year_list[k], ']\n', sep = '')


    df_final_richness=colSums(df_final)
    #problem
    #year1, sp = 1,2,3, rich=3
    #year2, sp = 1, rich = 4
    ratio_df=df_final_richness/raster_input@data@values[list_cells_df]

    vec_final[[k]]=ratio_df

    df_final=df_final+list_of_big_df[[year_list[k+1]]]
    df_final[df_final>1]=1
  }

  names(vec_final)=year_list

  x <- vec_final
  na_to_zero <- function(x) {
    if (length(x[[1]]) > 1) {
      res <- vector(mode = 'list', length = length(x))
      for (i in seq_along(x)) {
        res[[i]] <- na_to_zero(x[[i]])
      }
    } else {
      x[is.na(x)] <- 0
      res <- x
    }
    res
  }
  res <- na_to_zero(x)
  vec_final=res
  names(vec_final)=year_list




  #convert to df
  vec_final_df=lapply(vec_final, function(x) {if(any(class(x)=="matrix")) as.data.frame(x) else x})

  print('step 5')

  vec_final_df2=vec_final_df

  list_of_cells=names(vec_final_df2[[year_list[k]]])


  for (k in seq_along(year_list)){
    cat('.... .... .... year [', year_list[k], '/', tail(n=1,year_list), ']\n', sep = '')

    #make year column for df
    year_df=c(rep(year_list[k],length(vec_final_df[[year_list[k]]])))

    vec_final_df2[[year_list[k]]]=cbind(vec_final_df[[year_list[k]]],year_df)

    year=year_list[k]

    a_p_final=data.frame()

    for (i in seq_along(list_of_cells)){
      #cat('.... .... .... year [', list_of_cells[i], '/', list_of_cells[i], ']\n', sep = '')

      #add a_p column

      if (as.numeric(list_of_cells[i]) %in% num_list_visits[[year]]==T){
        a_p_temp=cbind(rownames(vec_final_df2[[year_list[k]]])[i],1)}

      else{
        a_p_temp=cbind(rownames(vec_final_df2[[year_list[k]]])[i],0)
      }

      a_p_final=rbind(a_p_final,a_p_temp, stringsAsFactors=FALSE)
    }

    print('step 6')

    length(a_p_final[,1])
    #go to next cell

    #finish dataframe
    a_p_final_temp=c(rep(0,length(list_of_cells)))
    vec_final_df2[[year_list[k]]]=cbind(vec_final_df2[[year_list[k]]],a_p_final_temp)
    vec_final_df2[[year_list[k]]][,3]=as.numeric(vec_final_df2[[year_list[k]]][,3]) + as.numeric(a_p_final$V2)

    #go to another year and repeat
  }


  #converting list of matrix into mother df
  df_final=do.call(rbind, vec_final_df2)

  df_final_backup=df_final

  df_final_test=cbind(df_final_backup,rownames(df_final_backup))

  df_final2=as.data.frame(df_final_test)

  colnames(df_final2)=c('ratio','year','a_p','cell_id')

  df_final2_back=df_final2
  df_final2$ratio=as.numeric(as.character(df_final2$ratio))
  df_final2$year=as.numeric(as.character(df_final2$year))
  df_final2$ratio[is.nan(df_final2$ratio)]=0
  df_final2$ratio[is.infinite(df_final2$ratio)]=0
  df_final2$cell_id=as.numeric(as.character(df_final2$cell_id))
  df_final2$a_p=as.numeric(as.character(df_final2$a_p))
  return(df_final2)
}
