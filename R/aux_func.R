get_adj_mat = function(county.adjacency,Countyvec,Statevec){
  if(!is.data.frame(county.adjacency)){stop("county.adjacency must be a data frame")}
  if(!is.character(Countyvec)){stop("Countyvec must be a character vector")}
  if(!is.character(Statevec)){stop("Statevec must be a character vector")}
  get(USAcities)
  if(prod(Countyvec %in% USAcities$county_name) != 1){"County names do not match"}
  if(prod(Statevec %in% USAcities$state_id) != 1){"State abbreviations do not match"}

  county.names1 <- trimws(unlist(lapply(county.adjacency$Countyname, function(x) unlist(strsplit(as.character(x), ","))[1])))
  county.names2 <- trimws(unlist(lapply(county.adjacency$neighborname, function(x) unlist(strsplit(as.character(x), ","))[1])))
  state.abbs1 <- trimws(unlist(lapply(county.adjacency$Countyname,
                                      function(x) tail(unlist(strsplit(as.character(x), ",")),n=1))))
  state.abbs2 <- trimws(unlist(lapply(county.adjacency$neighborname,
                                      function(x) tail(unlist(strsplit(as.character(x), ",")),n=1))))

  county.names1.full <- unlist(lapply(1:length(county.names1), function(x) paste(county.names1[x], state.abbs1[x], sep="-")))
  county.names2.full <- unlist(lapply(1:length(county.names2), function(x) paste(county.names2[x], state.abbs2[x], sep="-")))

  county.list <- paste(Countyvec,Statevec,sep="-")
  n.county <- length(county.list)

  county.adj.mat <- matrix(0, nrow=n.county, ncol=n.county)
  for (i in 1:dim(county.adjacency)[1]) {
    inx.county <- which(county.list == county.names1.full[i])
    inx.neighbor.county <- which(county.list == county.names2.full[i])
    county.adj.mat[inx.county, inx.neighbor.county] = 1
    county.adj.mat[inx.neighbor.county, inx.county] = 1
  }
  return(county.adj.mat)
}



ResultTableSummary = function(bstfit){
  if(is.null(bstfit$Alpha)){stop("bstfit must have a named component Alpha")}
  if(is.null(bstfit$Beta)){stop("bstfit must have a named component Beta")}

  alphamat = apply(bstfit$Alpha,c(1,2),mean)
  betamat  = apply(bstfit$Beta,c(1,2),mean)
  temp <- data.frame(cbind(alphamat,betamat))
  # colnames(temp) <- c(paste0("a",(c(1:ncol(alphamat))-1)),paste0("b",(c(1:ncol(betamat))-1)))
  colnames(temp) <- c(paste("a",colnames(alphamat),sep="."),paste("b",colnames(betamat),sep="."))

  temp %>% tbl_summary(statistic = list(all_continuous() ~ "{mean} ({p5}, {p95})", all_categorical() ~ "{n} ({p}%)"),
                       digits = list(all_continuous() ~ c(3,3))) %>%
    modify_header(label = "**Coefficients**",
                  stat_0 = '**Estimates**') %>%
    modify_caption("**BSTZINB Model Coefficients**") %>%
    bold_labels() %>%
    modify_footnote(
      all_stat_cols() ~ "Point estimates (90% credible intervals)")
}


ResultTableSummary2 = function(y, X, A, nt, nchain=3, nsim=10, nburn=0, nthin=1){

  if(!is.vector(y)){stop("y must be a vector")}
  if(!is.matrix(X)){stop("X must be a matrix")}
  if(!is.matrix(A)){stop("A must be a matrix")}
  if(nt <= 0){stop("nt must be positive integer")}
  if(nchain < 1){stop("nchain must be a positive integer")}
  if(nsim < 1){stop("nsim must be a positive integer")}
  if(nburn < 0){stop("nburn must be a non-negative integer")}
  if(nthin < 1){stop("nthin must be a positive integer")}
  if(sum(is.na(y))>0){naind <- which(is.na(y)); if(length(unique(y[-naind]))!=2)stop("y must have two categories")} else{if(length(unique(y))!=2)stop("y must have two categories")}
  y <- as.numeric(y)
  if(!is.numeric(X)){stop("X must be numeric")}
  if(!is.numeric(A)){stop("A must be numeric")}

  stfit4 = BSTZINB(y, X, A, nt, LinearT=FALSE, nchain, nsim, nburn, nthin)
  alphamat = apply(stfit4$Alpha,c(1,2),mean)
  betamat  = apply(stfit4$Beta,c(1,2),mean)
  temp4 <- data.frame(matrix(NA,nrow(alphamat),ncol(alphamat)+ncol(betamat)+2))
  # colnames(temp4) <- c(paste0("a",(c(1:ncol(alphamat))-1)),paste0("b",(c(1:ncol(betamat))-1)))
  colnames(temp4) <- c("a.t",paste("a",colnames(alphamat),sep="."),"b.t",paste("b",colnames(betamat),sep="."))
  temp4[,paste("a",colnames(alphamat),sep=".")] <- alphamat
  temp4[,paste("b",colnames(betamat),sep=".")]  <- betamat
  DIC4 = compute.ZINB.DIC(stfit4,(nsim-nburn)/nthin,nchain)
  ind4 = rep(NA,ncol(temp4)); names(ind4) <-  colnames(temp4)
  ind4[paste("a",colnames(alphamat),sep=".")] <- conv.test(stfit4$Alpha)
  ind4[paste("b",colnames(betamat),sep=".")] <- conv.test(stfit4$Beta)

  stfit3 = BSTZINB(y, X, A, nt, LinearT=TRUE, nchain, nsim, nburn, nthin)
  alphamat = apply(stfit3$Alpha,c(1,2),mean)
  betamat  = apply(stfit3$Beta,c(1,2),mean)
  temp3 <- data.frame(matrix(NA,nrow(temp4),ncol(temp4)))
  colnames(temp3) <- colnames(temp4)
  temp3[,paste("a",colnames(alphamat),sep=".")] <- alphamat
  temp3[,paste("b",colnames(betamat),sep=".")]  <- betamat
  DIC3 = compute.ZINB.DIC(stfit3,(nsim-nburn)/nthin,nchain)
  ind3 = rep(NA,length(ind4)); names(ind3) <- names(ind4)
  ind3[paste("a",colnames(alphamat),sep=".")] <- conv.test(stfit3$Alpha)
  ind3[paste("b",colnames(betamat),sep=".")] <- conv.test(stfit3$Beta)

  stfit2 = BSTNB(y, X, A, nt, nchain, nsim, nburn, nthin)
  betamat  = apply(stfit2$Beta,c(1,2),mean)
  temp2 <- data.frame(matrix(NA,nrow(temp4),ncol(temp4)))
  colnames(temp2) <- colnames(temp4)
  temp2[,paste("b",colnames(betamat),sep=".")] <- betamat
  DIC2 = compute.NB.DIC(stfit2,(nsim-nburn)/nthin,nchain)
  ind2 = rep(NA,length(ind4)); names(ind2) <- names(ind4)
  ind2[paste("b",colnames(betamat),sep=".")] <- conv.test(stfit2$Beta)

  stfit1 = BSTP(y, X, A, nt, nchain, nsim, nburn, nthin)
  betamat  = apply(stfit1$Beta,c(1,2),mean)
  temp1 <- data.frame(matrix(NA,nrow(temp4),ncol(temp4)))
  colnames(temp1) <- colnames(temp4)
  temp1[,paste("b",colnames(betamat),sep=".")] <- betamat
  DIC1 = compute.NB.DIC(stfit1,(nsim-nburn)/nthin,nchain)
  ind1 = rep(NA,length(ind4)); names(ind1) <- names(ind4)
  ind1[paste("b",colnames(betamat),sep=".")] <- conv.test(stfit1$Beta)

  tabout <- data.frame(rbind(temp4,temp3,temp2,temp1))
  tabout$var <- c(rep("4BSTZINB(NLT)",nrow(temp4)),
                  rep("3BSTZINB(LT)",nrow(temp3)),
                  rep("2BSTNB",nrow(temp2)),
                  rep("1BSTP",nrow(temp1)))

  table <- tabout%>% tbl_summary(by = var,missing = "no",
                                 digits = list(all_continuous() ~ c(3,3)))%>%
    modify_header(label = "**Coefficients**",
                  stat_1 = '**BSTP**',
                  stat_2 = '**BSTNB**',
                  stat_3 = '**BSTZINB(LT)**',
                  stat_4 = '**BSTZINB(NLT)**') %>%
    modify_footnote(
      all_stat_cols() ~ paste("Point estimates (90% credible intervals)",
                              "DIC1 = ",round(DIC1,1),";",
                              "DIC2 = ",round(DIC2,1),";",
                              "DIC3 = ",round(DIC3,1),";",
                              "DIC4 = ",round(DIC4,1)))

  table$table_body$stat_1[table$table_body$stat_1 == "NA (NA, NA)"] <- NA
  table$table_body$stat_2[table$table_body$stat_2 == "NA (NA, NA)"] <- NA
  table$table_body$stat_3[table$table_body$stat_3 == "NA (NA, NA)"] <- NA
  table$table_body$stat_4[table$table_body$stat_4 == "NA (NA, NA)"] <- NA



  table.fin <- table %>% as_gt %>%
    tab_style(
      style = list(cell_text(color = "darkred")),
      locations = cells_body(columns=stat_4,rows=(ind4==FALSE))
    ) %>%
    tab_style(
      style = list(cell_text(color = "darkred")),
      locations = cells_body(columns=stat_3,rows=(ind3==FALSE))
    )%>%
    tab_style(
      style = list(cell_text(color = "darkred")),
      locations = cells_body(columns=stat_2,rows=(ind2==FALSE))
    )%>%
    tab_style(
      style = list(cell_text(color = "darkred")),
      locations = cells_body(columns=stat_1,rows=(ind1==FALSE))
    )

  return(table.fin)

}


compute.ZINB.DIC = function(fit,lastit,nchain){

  if(is.null(fit$Eta1)){stop("fit must have a named component Eta1")}
  if(is.null(fit$Eta2)){stop("fit must have a named component Eta2")}
  if(is.null(fit$R)){stop("fit must have a named component R")}
  if(!is.numeric(lastit) | lastit <= 0){stop("lastit must be a positive integer")}
  if(!is.numeric(nchain) | nchain <= 0){stop("nchain must be a positive integer")}

  computeD.avg = function(fit){
    eta1.mean <- apply(fit$Eta1,2,mean)
    eta2.mean <- apply(fit$Eta2,2,mean)
    r.mean    <- mean(fit$R)
    pi = gtools::inv.logit(eta1.mean)
    q <- pmax(0.01,pmin(0.99,1/(1+exp(eta2.mean))))
    I = fit$I[dim(fit$I)[1],,nchain]
    dNB = dnbinom(y[I==1],r.mean,q[I==1],log=T)
    comp1 = log(1-pi[I==0])
    comp2 = (log(pi[I==1])+dNB)
    # (-2)*(sum(comp1)+sum(comp2))
    (-2)*sum(comp2)
  }
  computeD.indiv = function(fit,iter,chain){
    eta1 <- fit$Eta1[iter,,chain]
    eta2 <- fit$Eta2[iter,,chain]
    r    <- fit$R[iter,chain]
    pi = gtools::inv.logit(eta1)
    q <- pmax(0.01,pmin(0.99,1/(1+exp(eta2))))
    I = fit$I[iter,,chain]
    dNB = dnbinom(y[I==1],r,q[I==1],log=T)
    comp1 = log(1-pi[I==0])
    comp2 = (log(pi[I==1])+dNB)
    # (-2)*(sum(comp1)+sum(comp2))
    (-2)*sum(comp2)
  }
  Dmat = matrix(0,lastit,nchain)
  for(iter in 1:lastit){
    for(chain in 1:nchain){
      Dmat[iter,chain] = computeD.indiv(fit,iter,chain)
    }
  }

  comp1 = computeD.avg(fit)
  comp2 = mean(Dmat)
  DIC = comp1 + 2*( comp2 - comp1)
  return(DIC)
}


compute.NB.DIC = function(fit,lastit,nchain){

  if(is.null(fit$Eta1)){stop("fit must have a named component Eta1")}
  if(is.null(fit$R)){stop("fit must have a named component R")}
  if(!is.numeric(lastit) | lastit <= 0){stop("lastit must be a positive integer")}
  if(!is.numeric(nchain) | nchain <= 0){stop("nchain must be a positive integer")}

  computeD.avg = function(fit){
    eta.mean <- apply(fit$Eta1,2,mean)
    r.mean   <- mean(fit$R)
    q <- pmax(0.01,pmin(0.99,1/(1+exp(eta.mean))))
    dNB = dnbinom(y,r.mean,q,log=T)
    (-2)*sum(dNB)
  }
  computeD.indiv = function(fit,iter,chain){
    eta  <- fit$Eta1[iter,,chain]
    r    <- fit$R[iter,chain]
    q <- pmax(0.01,pmin(0.99,1/(1+exp(eta))))
    dNB = dnbinom(y,r,q,log=T)
    (-2)*sum(dNB)
  }
  Dmat = matrix(0,lastit,nchain)
  for(iter in 1:lastit){
    for(chain in 1:nchain){
      Dmat[iter,chain] = computeD.indiv(fit,iter,chain)
    }
  }

  comp1 = computeD.avg(fit)
  comp2 = mean(Dmat)
  DIC = comp1 + 2*( comp2 - comp1)
  return(DIC)

}


conv.test = function(params,nchain=3,thshold=1.96){

  if(!is.numeric(nchain) | nchain <= 0){stop("nchain must be a positive integer")}
  if(!is.numeric(thshold) | thshold <= 0){stop("thshold must be a positive number")}

  if(length(dim(params))==3){
    p = dim(params)[2]
    out = rep(NA,p)
    for(colid in 1:p){
      mcmc.temp = mcmc(params[,colid,])
      testout = geweke.diag(mcmc.temp) %>% unlist
      out[colid] <- as.logical(prod(abs(testout[1:nchain])<thshold))
    }
    return(out)
  }else{
    mcmc.temp = mcmc(params)
    testout = geweke.diag(mcmc.temp) %>% unlist
    as.logical(prod(abs(testout[1:nchain])<thshold))
  }

}
