USDmapCount = function(state.sel,dat,scol,tcol,tsel=NULL,cname,uplim=NULL){

  if(!is.character(state.sel)){stop("state.sel must be character vector")}
  if(!is.numeric(scol)){stop("scol must be numeric")}
  if(!is.numeric(tcol)){stop("tcol must be numeric")}
  if(!is.null(tsel)){if(!is.numeric(tsel)){stop("tsel must be NULL or numeric")}}
  if(!is.character(cname)){stop("cname must be character vector")}
  if(!is.null(uplim)){if(!is.numeric(uplim)){stop("uplim must be NULL or numeric")}}
  if(!is.data.frame(dat)){stop("dat must be a data frame")}
  if(is.null(dat$y)){stop("dat must have an entry named y")}
  if(is.null(dat$sid)){stop("dat must have an entry named sid")}
  if(is.null(dat$tid)){stop("dat must have an entry named tid")}
  if(prod(toupper(state.sel) %in% state.abb) != 1){stop("State abbreviation does not match")}
  get(USAcities)
  if(prod(toupper(cname) %in% toupper(USAcities$county_name)) != 1){stop("County names do not match")}

  if(is.null(uplim)){
    ulim = max(dat$y)
  }else{
    ulim = uplim
  }
  # library(usmap)
  # all.st <- usmap::us_map(regions="counties")
  # state.map <- all.st[all.st$abbr %in% state.sel,]
  stateabb <- toupper(state.sel)
  statefull <- tolower(state.name[which(state.abb == stateabb)])
  state.map <- map_data("county",region=statefull)
  sid = names(dat) %>% .[scol]
  tid = names(dat) %>% .[tcol]
  if(is.null(tsel)){
    zinb.summary     = dat %>% group_by(sid) %>% summarise(cnt=mean(y))
  }else{
    zinb.summary     = dat %>% group_by(sid) %>% filter(tid==tsel)
    zinb.summary$cnt = zinb.summary$y
    zinb.summary = zinb.summary %>% dplyr::select(sid,cnt)
  }

  zinb.summary$cid = tolower(cname[zinb.summary$sid])
  map.df  <- left_join(state.map, zinb.summary, by = c("subregion"="cid"))
  p <- ggplot(data = map.df, aes(x = long, y = lat, group=group)) +
    geom_polygon(aes(fill = cnt), color="black") +
    scale_fill_gradientn(colours = viridis::viridis(25))+
    theme_bw() +
    labs(fill= "",
         title = "", x="", y="")
  return(p)

}

qRankPar = function(state.set,ns,nt,cname,stfit,vn=12){

  if(!is.character(state.set)){stop("state.set must be character vector")}
  if(!is.numeric(ns)){stop("ns must be numeric")}
  if(!is.numeric(nt)){stop("nt must be numeric")}
  if(!is.character(cname)){stop("cname must be character vector")}
  if(is.null(stfit$Eta1)){stop("stfit must have an entry named Eta1")}
  if(!is.numeric(vn)){stop("vn must be a positive integer")}
  if(vn <= 0){stop("vn must be a positive integer")}
  if(prod(toupper(state.set) %in% state.abb) != 1){stop("State abbreviation does not match")}
  get(USAcities)
  if(prod(toupper(cname) %in% toupper(USAcities$county_name)) != 1){stop("County names do not match")}



  qij.mat <- matrix(inv.logit(apply(stfit$Eta1,2,mean)),nrow=ns)
  zinb.summary <- data.frame(County=cname, m = rowMeans(qij.mat))
  zinb.summary <- zinb.summary[order(zinb.summary$m),]
  zinb.summary$County = factor(zinb.summary$County)
  zinb.summary.sample = zinb.summary[floor(seq(1,ns,length.out=vn)),]

  par(mfrow=c(1,1),mar=c(3,5,1,1))
  p <- ggplot(zinb.summary.sample,aes(x=reorder(County, m), y=m, fill=County)) + geom_bar(alpha=0.8,stat="identity") +
    xlab("") + ylab("Probability at risk") + ylim(0,1) + theme_bw() + theme(legend.position = "")
  return(p + coord_flip()+theme(axis.title.x = element_text(size=24),
                                axis.text.x = element_text(size=24),
                                axis.text.y = element_text(size=24)))
}

qRankParTop = function(state.set,ns,nt,cname,stfit,vn=12){

  if(!is.character(state.set)){stop("state.set must be character vector")}
  if(!is.numeric(ns)){stop("ns must be numeric")}
  if(!is.numeric(nt)){stop("nt must be numeric")}
  if(!is.character(cname)){stop("cname must be character vector")}
  if(is.null(stfit$Eta1)){stop("stfit must have an entry named Eta1")}
  if(!is.numeric(vn)){stop("vn must be a positive integer")}
  if(vn <= 0){stop("vn must be a positive integer")}
  if(prod(toupper(state.set) %in% state.abb) != 1){stop("State abbreviation does not match")}
  get(USAcities)
  if(prod(toupper(cname) %in% toupper(USAcities$county_name)) != 1){stop("County names do not match")}

  qij.mat <- matrix(inv.logit(apply(stfit$Eta1,2,mean)),nrow=ns)
  zinb.summary <- data.frame(County=cname, m = rowMeans(qij.mat))
  zinb.summary <- zinb.summary[order(zinb.summary$m,decreasing = TRUE),]
  zinb.summary$County = factor(zinb.summary$County)
  zinb.summary.sample = zinb.summary[c(1:vn),]

  par(mfrow=c(1,1),mar=c(3,5,1,1))
  p <- ggplot(zinb.summary.sample,aes(x=reorder(County, m), y=m, fill=County)) + geom_bar(alpha=0.8,stat="identity") +
    xlab("") + ylab("Probability at risk") + ylim(0,1) + theme_bw() + theme(legend.position = "")
  return(p + coord_flip()+theme(axis.title.x = element_text(size=24),
                                axis.text.x = element_text(size=24),
                                axis.text.y = element_text(size=24)))
}

TimetrendCurve = function(stfit,ns,nt,countyname,vn=5,smooth.mode=TRUE){

  if(!is.numeric(ns)){stop("ns must be numeric")}
  if(!is.numeric(nt)){stop("nt must be numeric")}
  if(!is.character(countyname)){stop("countyname must be character vector")}
  if(is.null(stfit$Eta1)){stop("stfit must have an entry named Eta1")}
  if(!is.numeric(vn)){stop("vn must be a positive integer")}
  if(vn <= 0){stop("vn must be a positive integer")}
  if(!is.logical(smooth.mode)){stop("smooth.mode must be TRUE/FALSE")}
  get(USAcities)
  if(prod(toupper(countyname) %in% toupper(USAcities$county_name)) != 1){stop("County names do not match")}


  time = c(1:nt)
  df   = data.frame(matrix(apply(stfit$Eta1,2,mean),nrow=ns)); rownames(df) <- factor(countyname)
  df2 = data.frame(time,t(df[seq(1,ns,length.out=vn),]))
  if(smooth.mode){
    library(splines)
    time = spline(df2[,1],df2[,2])$x
    df2 <- data.frame(time,apply(df2[,-1],2,function(w) {spline(df2[,1],w)$y}))
  }
  library(reshape2)
  dd = melt(df2,c("time"))

  spd.plot <- ggplot(dd, aes(x=time,y=value)) +
    geom_line(aes(colour = variable, group = variable),linewidth=1.2)+
    geom_line(aes(time,0),color="red",linewidth=1.2,lty=2) + ylab("") + xlab("") +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.x = element_text(size=20))

  return(spd.plot)

}


