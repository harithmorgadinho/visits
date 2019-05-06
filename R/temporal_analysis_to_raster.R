#' @title temporal_analysis_to_raster
#' @description A function to plot biodiversity data.
#' @param input  df or raster
#' @param crs projection of data provided ('longlat'/'cea'/'auto')
#' @param zoom_out zoom level for the map
#' @examples
#' crop_map_world(df)
#' crop_map_world(df,crs='longlat')
#' @export
#'@importFrom BBmisc normalize
#'@importFrom raster raster

temporal_analysis_to_raster= function(input,from=NULL,to=NULL,ratio_low=NULL,ratio_high=NULL,var_2=NULL,remove_ratio_low=NULL,remove_ratio_high=NULL) {

  #input$ratio[which(is.na(input$ratio))]=0

  df_temp=input %>% group_by(cell_id) %>% summarise(a_p_sum = sum(a_p))
  cell_list_fixed=df_temp$cell_id[which(df_temp$a_p_sum>=5)]
  input= input %>% filter(cell_id%in%cell_list_fixed)
  input=input[input$cell_id%in%cell_list_fixed,]


  if(is.null(from) & is.null(to)){
    input=input}

  if(!is.null(from) & is.null(to)){

    input=input %>% filter(year>=from)
    df_temp=input %>% group_by(cell_id) %>% summarise(a_p_sum = sum(a_p))
    cell_list_fixed=df_temp$cell_id[which(df_temp$a_p_sum>=5)]
    input= input %>% filter(cell_id%in%cell_list_fixed)
    input=input[input$cell_id%in%cell_list_fixed,]

  }

  if(is.null(from) & !is.null(to)){
    input=input %>% filter(year<=to)
    df_temp=input %>% group_by(cell_id) %>% summarise(a_p_sum = sum(a_p))
    cell_list_fixed=df_temp$cell_id[which(df_temp$a_p_sum>=5)]
    input= input %>% filter(cell_id%in%cell_list_fixed)
    input=input[input$cell_id%in%cell_list_fixed,]
  }

  if(!is.null(from) & !is.null(to)){
    input=input %>% filter(year<=to) %>% filter(year>=from)
    df_temp=input %>% group_by(cell_id) %>% summarise(a_p_sum = sum(a_p))
    cell_list_fixed=df_temp$cell_id[which(df_temp$a_p_sum>=5)]
    input= input %>% filter(cell_id%in%cell_list_fixed)
    input=input[input$cell_id%in%cell_list_fixed,]
  }

  if(is.null(ratio_low) & is.null(ratio_high)){
    input=input}

  if(!is.null(ratio_low) & is.null(ratio_high)){
    keep=unique(input$cell_id[which(input$ratio>=ratio_low)])
    input=input %>%filter(cell_id%in%keep)
  }
  if(is.null(ratio_low) & !is.null(ratio_high)){

    keep=unique(input$cell_id[which(input$ratio<=ratio_high)])
    input=input %>%filter(cell_id%in%keep)
  }
  if(!is.null(ratio_low) & !is.null(ratio_high)){

    keep=unique(subset(input$cell_id, input$ratio>=ratio_low & input$ratio<=ratio_high))
    input=input %>%filter(cell_id%in%keep)
  }



  if(!is.null(remove_ratio_low) & !is.null(remove_ratio_high)){
    df=input %>% group_by(cell_id) %>% summarise(ratio = max(ratio))
    keep=(df$cell_id[which(df$ratio>=remove_ratio_low & df$ratio<=remove_ratio_high)])
    input=input %>%filter(cell_id%in%keep)
  }

  if(is.null(remove_ratio_low) & is.null(remove_ratio_high)){
    input=input}

  if(!is.null(remove_ratio_low) & is.null(remove_ratio_high)){
    df=input %>% group_by(cell_id) %>% summarise(ratio = max(ratio))
    keep=df$cell_id[which(df$ratio>=remove_ratio_low)]
    input=input %>%filter(cell_id%in%keep)

  }

  if(is.null(remove_ratio_low) & !is.null(remove_ratio_high)){
    df=input %>% group_by(cell_id) %>% summarise(ratio = max(ratio))
    keep=df$cell_id[which(df$ratio<=remove_ratio_high)]
    input=input %>%filter(cell_id%in%keep)
  }



  if (nrow(input)==0){

    r=raster(ncol=84, nrow=77,extent(shp_behrmann), resolution=100000,crs=CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    values(r)=0
    return(r)
  }
  else{


    df_final=data.frame()

    cell_id_list=unique(input$cell_id)
    print('step 7')

    input_stand=normalize(input[,c(1,2)], method = "standardize", range = c(-1, 1), margin = 1L, on.constant = "quiet")
    input_stand$cell_id=input$cell_id
    input_stand$a_p=input$a_p

    for (j in seq_along(cell_id_list)){

      cell_id_temp=cell_id_list[[j]]
      print(cell_id_temp)
      subset=input_stand %>% filter(cell_id==cell_id_temp)

      if(sum(input_stand$a_p)==0){
        print(cell_id_temp)
      }else{

        # lm model----

        #lm<-glm(a_p~ratio+year,data=subset,family='binomial')
        lm<-tryCatch(glm(a_p~ratio+year,data=subset,family='binomial'),error=function(e) e, warning=function(w) w)

        a=as.data.frame(lm$coefficients)

        if(!is.null(var_2)){

          lm_value=a$`lm$coefficients`[2]
          lm_value_b=a$`lm$coefficients`[3]
          if(is(lm,"warning")) lm_value=NA
          if(is(lm,"warning")) lm_value_b=NA
          df_temp=cbind(cell_id_temp,lm_value,lm_value_b)
          df_final=rbind(df_final,df_temp)


        }else{

          lm_value=a$`lm$coefficients`[2]
          if(is(lm,"warning")) lm_value=NA
          df_temp=cbind(cell_id_temp,lm_value)
          df_final=rbind(df_final,df_temp)
        }

      }
    }

    r=raster(ncol=84, nrow=77,extent(shp_behrmann), resolution=100000,crs=CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

    df=as.data.frame(1:ncell(r))
    df$values=NA
    colnames(df)=c('cell_id_temp','values')
    merged=merge(df,df_final,by='cell_id_temp',all=TRUE)
    values(r)=merged$lm_value

    if(!is.null(var_2)){
      r2=r
      df=as.data.frame(1:ncell(r))
      df$values=NA
      colnames(df)=c('cell_id_temp','values')
      merged=merge(df,df_final,by='cell_id_temp',all=TRUE)
      values(r2)=merged$lm_value_b
      #r_list=as.list(r,r2)
      r_list=list(r,r2)
      return(r_list)

    }else{

      #plot(r)
      return(r)}}
}
