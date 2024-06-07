#' @title Draw spatial maps of various quantities over regions in the US
#'
#' @description
#' Creates a map of any given quantity (at a selected time or averaged over time) for regions in the US specified by state and county
#'
#' @usage USDmapCount(state.sel,dat,scol,tcol=NULL,tsel=NULL,cname,uplim=NULL)
#'
#' @param state.sel character vector giving the selected states
#' @param dat data frame having named components: y - the necessary quantity (numeric), sid - the region indices, tid - the time indices
#' @param scol column index of the region to choose
#' @param tcol (optional) column index of the time point to choose
#' @param tsel (optional) selected time point
#' @param cname character vector of county names, must match those in USAcities
#' @param uplim (optional) numeric, upper limit for the given quantity
#'
#' @import ggplot2
#' @import maps
#' @import dplyr
#' @import viridis
#' @import datasets
#' @import utils
#'
#' @return spatial map of the required quantity over the specified region
#' @export
#'
USDmapCount = function(state.sel,dat,scol,tcol=NULL,tsel=NULL,cname,uplim=NULL){

  if(!is.character(state.sel)){stop("state.sel must be character vector")}
  if(!is.numeric(scol)){stop("scol must be numeric")}
  if(!is.null(tcol)){if(!is.numeric(tcol)){stop("tcol must be NULL or numeric")}}
  if(!is.null(tsel)){if(!is.numeric(tsel)){stop("tsel must be NULL or numeric")}}
  if(!is.null(tsel)){if(is.null(tcol)){stop("tcol must be supplied if tsel is not null")}}
  if(!is.character(cname)){stop("cname must be character vector")}
  if(!is.null(uplim)){if(!is.numeric(uplim)){stop("uplim must be NULL or numeric")}}
  if(!is.data.frame(dat)){stop("dat must be a data frame")}
  if(is.null(dat$y)){stop("dat must have an entry named y")}
  if(is.null(dat$sid)){stop("dat must have an entry named sid")}
  if(is.null(dat$tid)){stop("dat must have an entry named tid")}
  if(prod(toupper(state.sel) %in% datasets::state.abb) != 1){stop("State abbreviation does not match")}
  USAcities <- BSTZINB::USAcities
  if(prod(toupper(cname) %in% toupper(USAcities$county_name)) != 1){stop("County names do not match")}

  if(is.null(uplim)){
    ulim = max(dat$y)
  }else{
    ulim = uplim
  }
  stateabb <- toupper(state.sel)
  statefull <- tolower(datasets::state.name[which(datasets::state.abb == stateabb)])
  state.map <- map_data("county",region=statefull)
  sid = names(dat) %>% .[scol]
  tid = names(dat) %>% .[tcol]
  if(is.null(tsel)){
    zinb.summary     = dat %>% dplyr::group_by(sid) %>% dplyr::summarise("cnt"=mean(dat$y))
  }else{
    zinb.summary     = dat %>% dplyr::group_by(sid) %>% dplyr::filter(tid==tsel)
    zinb.summary$cnt = zinb.summary$y
    zinb.summary = zinb.summary %>% dplyr::select(zinb.summary$sid,zinb.summary$cnt)
  }

  zinb.summary$cid = tolower(cname[zinb.summary$sid])
  map.df  <- dplyr::left_join(state.map, zinb.summary, by = c("subregion"="cid"))
  p <- ggplot2::ggplot(data = map.df, aes(x = map.df$long, y = map.df$lat, group=map.df$group)) +
    ggplot2::geom_polygon(aes(fill = map.df$cnt), color="black") +
    ggplot2::scale_fill_gradientn(colours = viridis::viridis(25))+
    ggplot2::theme_bw() +
    ggplot2::labs(fill= "",
         title = "", x="", y="")
  return(p)

}

#' @title Bar plot for time-averaged log-q estimates over quantile-representative counties (descending order)
#'
#' @description
#' Produce a descending order of bar plot for time-averaged log-q estimates over quantile-representative counties
#'
#' @usage qRankPar(state.set,ns,nt,cname,stfit,vn=12)
#'
#' @param state.set character vector of set of states on which the the graphics is to be made
#' @param ns positive integer, the number of counties
#' @param nt positive integer, the number of timepoints
#' @param cname character vector of the names of the counties
#' @param stfit the fitted data for BSTP, BSTNB or BSTZINB
#' @param vn positive integer, number of sample counties to display
#'
#' @import dplyr
#' @import boot
#' @import ggplot2
#' @import datasets
#' @import graphics
#' @import utils
#'
#' @return bar graph
#' @export
qRankPar = function(state.set,ns,nt,cname,stfit,vn=12){

  if(!is.character(state.set)){stop("state.set must be character vector")}
  if(!is.numeric(ns)){stop("ns must be numeric")}
  if(!is.numeric(nt)){stop("nt must be numeric")}
  if(!is.character(cname)){stop("cname must be character vector")}
  if(is.null(stfit$Eta1)){stop("stfit must have an entry named Eta1")}
  if(!is.numeric(vn)){stop("vn must be a positive integer")}
  if(vn <= 0){stop("vn must be a positive integer")}
  if(prod(toupper(state.set) %in% datasets::state.abb) != 1){stop("State abbreviation does not match")}
  USAcities <- BSTZINB::USAcities
  if(prod(toupper(cname) %in% toupper(USAcities$county_name)) != 1){stop("County names do not match")}



  qij.mat <- matrix(boot::inv.logit(apply(stfit$Eta1,2,mean)),nrow=ns)
  zinb.summary <- data.frame("County"=cname, "m" = rowMeans(qij.mat))
  zinb.summary <- zinb.summary[order(zinb.summary$m),]
  zinb.summary$County = factor(zinb.summary$County)
  zinb.summary.sample = zinb.summary[floor(seq(1,ns,length.out=vn)),]

  par(mfrow=c(1,1),mar=c(3,5,1,1))
  p <- ggplot2::ggplot(zinb.summary.sample,aes(x=reorder(zinb.summary$County, zinb.summary$m), y=zinb.summary$m, fill=zinb.summary$County)) + ggplot2::geom_bar(alpha=0.8,stat="identity") +
    ggplot2::xlab("") + ggplot2::ylab("Probability at risk") + ggplot2::ylim(0,1) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "")
  return(p + ggplot2::coord_flip()+ggplot2::theme(axis.title.x = element_text(size=24),
                                axis.text.x = element_text(size=24),
                                axis.text.y = element_text(size=24)))
}

#' @title Bar plot for time-averaged log-q estimates over top ranking counties (descending order)
#'
#' @description
#' Produce a descending order of bar plot for time-averaged log-q estimates over top ranking counties
#'
#' @usage qRankParTop(state.set,ns,nt,cname,stfit,vn=12)
#'
#' @param state.set character vector of set of states on which the the graphics is to be made
#' @param ns positive integer, the number of counties
#' @param nt positive integer, the number of timepoints
#' @param cname character vector of the names of the counties
#' @param stfit the fitted data for BSTP, BSTNB or BSTZINB
#' @param vn positive integer, number of sample counties to display
#'
#' @import dplyr
#' @import BayesLogit
#' @import ggplot2
#' @import datasets
#' @import graphics
#' @import utils
#'
#' @return bar graph
#' @export
qRankParTop = function(state.set,ns,nt,cname,stfit,vn=12){

  if(!is.character(state.set)){stop("state.set must be character vector")}
  if(!is.numeric(ns)){stop("ns must be numeric")}
  if(!is.numeric(nt)){stop("nt must be numeric")}
  if(!is.character(cname)){stop("cname must be character vector")}
  if(is.null(stfit$Eta1)){stop("stfit must have an entry named Eta1")}
  if(!is.numeric(vn)){stop("vn must be a positive integer")}
  if(vn <= 0){stop("vn must be a positive integer")}
  if(prod(toupper(state.set) %in% datasets::state.abb) != 1){stop("State abbreviation does not match")}
  USAcities <- BSTZINB::USAcities
  if(prod(toupper(cname) %in% toupper(USAcities$county_name)) != 1){stop("County names do not match")}

  qij.mat <- matrix(boot::inv.logit(apply(stfit$Eta1,2,mean)),nrow=ns)
  zinb.summary <- data.frame("County"=cname, "m" = rowMeans(qij.mat))
  zinb.summary <- zinb.summary[order(zinb.summary$m,decreasing = TRUE),]
  zinb.summary$County = factor(zinb.summary$County)
  zinb.summary.sample = zinb.summary[c(1:vn),]

  par(mfrow=c(1,1),mar=c(3,5,1,1))
  p <- ggplot2::ggplot(zinb.summary.sample,aes(x=reorder(zinb.summary$County, zinb.summary$m), y=zinb.summary$m, fill=zinb.summary$County)) + ggplot2::geom_bar(alpha=0.8,stat="identity") +
    ggplot2::xlab("") + ggplot2::ylab("Probability at risk") + ggplot2::ylim(0,1) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "")
  return(p + ggplot2::coord_flip()+ggplot2::theme(axis.title.x = element_text(size=24),
                                axis.text.x = element_text(size=24),
                                axis.text.y = element_text(size=24)))
}

#' @title Time-trend curve over the study time domain for counties in the US
#'
#' @description
#' Produce a time-trend curve over the study time domain for counties in the US
#'
#' @usage TimetrendCurve(stfit,ns,nt,countyname,vn=5,smooth.mode=TRUE)
#'
#' @param stfit fitted object from BSTP, BSTNB or BSTZINB
#' @param ns positive integer, the number of counties
#' @param nt positive integer, the number of timepoints
#' @param countyname character vector of county names to use
#' @param vn positive integer, number of sample counties to use
#' @param smooth.mode logical, should splines be fitted to make it smooth
#'
#' @import datasets
#' @import ggplot2
#' @import dplyr
#' @import splines
#' @import reshape
#' @import utils
#'
#' @return time-trend curves
#' @export
TimetrendCurve = function(stfit,ns,nt,countyname,vn=5,smooth.mode=TRUE){

  if(!is.numeric(ns)){stop("ns must be numeric")}
  if(!is.numeric(nt)){stop("nt must be numeric")}
  if(!is.character(countyname)){stop("countyname must be character vector")}
  if(is.null(stfit$Eta1)){stop("stfit must have an entry named Eta1")}
  if(!is.numeric(vn)){stop("vn must be a positive integer")}
  if(vn <= 0){stop("vn must be a positive integer")}
  if(!is.logical(smooth.mode)){stop("smooth.mode must be TRUE/FALSE")}
  USAcities <- BSTZINB::USAcities
  if(prod(toupper(countyname) %in% toupper(USAcities$county_name)) != 1){stop("County names do not match")}


  time = c(1:nt)
  df   = data.frame(matrix(apply(stfit$Eta1,2,mean),nrow=ns)); rownames(df) <- factor(countyname)
  df2 = data.frame(time,t(df[seq(1,ns,length.out=vn),]))
  if(smooth.mode){
    time = spline(df2[,1],df2[,2])$x
    df2 <- data.frame(time,apply(df2[,-1],2,function(w) {spline(df2[,1],w)$y}))
  }
  dd = reshape::melt(df2,c("time"))

  spd.plot <- ggplot2::ggplot(dd, aes(x=time,y=dd$value)) +
    ggplot2::geom_line(aes(colour = dd$variable, group = dd$variable),linewidth=1.2)+
    ggplot2::geom_line(aes(time,0),color="red",linewidth=1.2,lty=2) + ggplot2::ylab("") + ggplot2::xlab("") +
    ggplot2::theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.x = element_text(size=20))

  return(spd.plot)

}


