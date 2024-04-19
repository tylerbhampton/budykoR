#' @title Simulated Data with Fu Budyko omega=2
#' @return 20 years of simulated data from a Budyko Fu formulation
#' and omega parameter value equal to 2 and average vertical fit error equal
#' to 0.05
#' @docType data
#' @md

"testdata"

#' @title Blank Budyko Curve Graph
#' @return Generates a blank Budyko Curve with energy and water limit
#' lines
#' @docType data
#' @examples
#' library(ggplot2)
#' blankBC+coord_cartesian(xlim=c(0,2))+geom_vline(xintercept=1,lty=2)+
#'     geom_text(data=data.frame(
#'         labels=c("Energy Limited","Water Limited"),
#'         x=c(0.7,1.3),y=c(0.25,0.25)
#'         ),aes(x=x,y=y,label=labels))
#' fit=as.data.frame(c("param"=2,"err"=0.5,"hshift"=0))
#' blankBC+
#'     geom_line(data=budyko_sim(fit=fit,method="Fu"))
#' @md

"blankBC"

#' @title Simulate Budyko Curve
#' @return Given a Budyko Curve type and fit parameter,
#' simulates a smooth curve in Budyko-space
#' @param p [numeric] parameter value for custom Budyko curve
#' @param fit [list] dataframe object with names() method and entries
#' "param" for parameter p, "err" for vertical fit error, and "hshift" for
#' optional hshift value (default to zero)
#' @param method A string. One of the following: Fu, Turc-Pike, Wang-Tank
#' @param res [numeric] resolution of fit, defaults to 0.01 units
#' @param hshift A boolean.  T/F (default F), whether to test for horizontal shift
#' @param hs [numeric] optional horizontal shift value
#' @param silent A boolean.  (default F) display error message?
#' @examples
#' library(ggplot2)
#' fit=as.data.frame(c("param"=2,"err"=0.5,"hshift"=0))
#' blankBC+
#'     geom_line(data=budyko_sim(fit=fit,method="Fu"))+
#'     coord_cartesian(xlim=c(0,5))+
#'     theme_bw()
#' @md
#' @export

budyko_sim=function(p=2,fit=NULL,method="budyko",res=0.01,hshift=FALSE,hs=NULL,silent=FALSE){
  if(!is.null(fit)){p=fit[1,]}
  if(!is.null(fit) & is.null(hs) & hshift==TRUE){hs=fit[4,]}
  if(!is.null(fit) & is.null(method)){
    if(tolower(names(fit))%in%c("fu","turc-pile","wang-tang")){
      method=names(fit)
    }
  }
  PET.P=seq(0,20,res)
  if(hshift==FALSE){
    if(tolower(method)=="budyko"){AET.P=((1-exp(-PET.P))*PET.P*tanh(1/PET.P))^(1/p)}
    if(tolower(method)=="fu"){AET.P=1+PET.P-(1+(PET.P)^p)^(1/p)}
    if(tolower(method)=="turc-pike"){AET.P=(1+(PET.P)^(-p))^(-1/p)}
    if(tolower(method)=="wang-tang"){AET.P=(1+PET.P-((1+PET.P)^2-4*p*(2-p)*PET.P)^(1/2))/(2*p*(2-p))}
    if(tolower(method)=="zhang"){AET.P=(1+p*PET.P)/(1+p*PET.P+(PET.P^(-1)))}
  }else{
    if(tolower(method)=="fu"){AET.P=1+(PET.P-hs)-(1+(PET.P-hs)^p)^(1/p)}
    if(tolower(method)=="turc-pike"){AET.P=(1+(PET.P-hs)^(-p))^(-1/p)}
    if(tolower(method)=="wang-tang"){AET.P=(1+PET.P-((1+PET.P)^2-4*p*(2-p)*PET.P)^(1/2))/(2*p*(2-p))}
    if(tolower(method)=="zhang"){AET.P=(1+p*PET.P)/(1+p*PET.P+(PET.P^(-1)))}
    AET.P[PET.P<hs]=NA
  }
  return(data.frame(PET.P,AET.P=round(AET.P,5)))
}


#' @title Fit Budyko Curve to Data
#' @return Given a dataset of Aridity Index and Evaporative Index,
#' find the best-fit Budyko Curve. Returns a "fit" object dataframe.
#' @param data [list] dataframe requiring the following two columns: "PET.P" for
#' potential evapotranspiration divided by precipitation (units normalized,
#' dimensionless), and "AET.P" for actual evapotranspiration divided by precipitation
#' or one minus the runoff coefficient (Q/P)
#' @param method A string. One of the following: Fu, Turc-Pike, Wang-Tang
#' @param dif A string. One of the following: "nls" (default) "rsq" "mae". Determines whether to fit automatically with stats::nls, or manually, with best fit determined by either Pearson's R squared or Minimum Mean Absolute Error
#' @param res [numeric] resolution of fit, defaults to 0.01 units
#' @param hshift A boolean.  T/F (default F), whether to test for horizontal shift
#' @param hs [numeric] optional horizontal shift value
#' @param silent A boolean.  (default F) display messages?
#' @examples
#' library(ggplot2)
#' fit=budyko_fit(data=testdata,method="Fu")
#' blankBC+
#'   geom_line(data=budyko_sim(fit=fit))+
#'   geom_point(data=testdata)+
#'   coord_cartesian(xlim=c(0,5))
#' @md
#' @export

budyko_fit=function(data,method,dif="nls",res=NULL,hshift=FALSE,hs=NULL,silent=FALSE){
  if(!(tolower(method) %in% c("wang-tang","fu","turc-pike","zhang")) & silent==FALSE){
    print("Error: unrecognized method")
  }else{
    if(silent==FALSE){print(paste0("Budyko Fit: method = ",method))}
    if(is.null(res)){res=0.01}
    data=data[,c("AET.P","PET.P")]
    data = subset(data,!is.na(AET.P))
    data = subset(data,PET.P<=20)
    data = subset(data,!(AET.P>PET.P))
    data=data[order(data$PET.P),]

    if(dif=="nls"){
      formula = list(
        "fu"="AET.P~1+PET.P-(1+(PET.P)^p)^(1/p)",
        "turc-pike"="AET.P~(1+(PET.P)^(-p))^(-1/p)",
        "wang-tang"="AET.P~(1+PET.P-((1+PET.P)^2-4*p*(2-p)*PET.P)^(1/2))/(2*p*(2-p))",
        "zhang"="AET.P~(1+p*PET.P)/(1+p*PET.P+(PET.P^(-1)))"
      )[[tolower(method)]]
      startval = list(
        "fu"=list(p=2.7),
        "turc-pike"=list(p=1.9),
        "wang-tang"=list(p=0.59),
        "zhang"=list(p=1.4)
      )[[tolower(method)]]
      budykoNLS = stats::nls(formula = formula,
                             start = startval,
                             data = data)
      paramval = summary(budykoNLS)$coefficients[[1]]
      fitmae = mean(abs(stats::residuals(budykoNLS)))
      fitrsq = summary(stats::lm(data$AET.P~fitted(budykoNLS)))$r.squared
      fitdat=as.data.frame(c(
        "param"=paramval,
        "mae"=round(fitmae,4),
        "rsq"=round(fitrsq,4),
        "hs"=0))
    }
    if(dif %in% c("rsq","mae")){
      if(tolower(method)=="fu"){testval=c(seq(1,3,res),seq(3,10,min(0.1,res*10)))}
      if(tolower(method)=="turc-pike"){testval=c(seq(0,2,res),seq(2,10,min(0.1,res*10)))}
      if(tolower(method)=="wang-tang"){testval=c(seq(-10,(-res),min(0.1,res*10)),seq(res,1,res))}
      if(tolower(method)=="zhang"){testval=c(seq(-0.05,0.99,0.01),seq(1,5,0.1))}
      if(hshift==FALSE){
        fiterr=do.call("rbind",lapply(testval,function(p){
          fit=as.data.frame(c("param"=p,"mae"=0,"rsq"=0,hshift=0))
          test=budyko_sim(fit=fit,method=method,res=res)
          AETPest=sapply((data$PET.P),function(z){test$AET.P[round(test$PET.P,-log10(res))==round(z,-log10(res))]})
          return(data.frame(
            rsq=1-summary(stats::lm(AETPest~data$AET.P))$r.squared,
            mae=mean(abs(data$AET.P-AETPest))
          ))
        }))
        wch=which.min(fiterr[,dif])[1]
        fitdat=as.data.frame(c(
          "param"=testval[wch],
          "mae"=round(fiterr$mae[wch],4),
          "rsq"=1-round(fiterr$rsq[wch],4),
          "hs"=0))
      }else{
        if(is.null(hs)){
          hs=c(seq(0,1,0.1),seq(1,ifelse(mean(data$PET.P,na.rm=TRUE)<=1,1,ceiling(mean(data$PET.P,na.rm=TRUE))),0.1))
        }else{hs=hs}
        
        difM=matrix(data=NA,nrow=length(hs),ncol=length(testval))
        for(r in 1:nrow(difM)){
          difM[r,]=sapply(testval,function(p){
            fit=as.data.frame(c("param"=p,"err"=0,hshift=hs[r]))
            test=budyko_sim(fit=fit,method=method,res=res,hshift=TRUE)
            difma=mean(abs(data$AET.P-sapply((data$PET.P),function(z){test$AET.P[round(test$PET.P,-log10(res))==round(z,-log10(res))]})))
            return(difma)
          })
        }
        fitdat=as.data.frame(c("param"=testval[which(difM==min(difM,na.rm=TRUE),arr.ind=TRUE)[1,][2]],"err"=round(min(difM,na.rm=TRUE)[1],4),"hshift"=hs[which(difM==min(difM,na.rm=TRUE),arr.ind=TRUE)[1,][1]]))
      }
    }
    names(fitdat)=method
    return(fitdat)
  }
}




#' @title Generate errorbounds from a Budyko fit
#' @return Given a "fit" dataframe object from budyko::budyko_fit,
#' generate a dataframe of error bounds that can be fed to plot or
#' ggplot2::ggplot. name (as in names()) of fit must be the method
#' @param fit [list] dataframe object with names() method and entries
#' "param" for parameter p, "err" for vertical fit error, and "hshift" for
#' optional hshift value (default to zero)
#' @param res [numeric] resolution of fit, defaults to 0.01 units
#' @param hshift A boolean.  T/F (default F), whether to test for horizontal shift.
#' hs value must be defined in fit
#' @param dif A string. One of the following: "nls" (default) "rsq" "mae". Determines whether 
#' to fit automatically with stats::nls, or manually, with best fit determined by either 
#' Pearson's R squared or Minimum Mean Absolute Error
#' @param data [list] dataframe requiring the following two columns: "PET.P" for
#' potential evapotranspiration divided by precipitation (units normalized,
#' dimensionless), and "AET.P" for actual evapotranspiration divided by precipitation
#' or one minus the runoff coefficient (Q/P)
#' @param alpha [numeric] If using dif="nls", confidence value for bounds. Default to 0.05 (95% confidence)
#' @param method A string. One of the following: Fu, Turc-Pike, Wang-Tang
#' @examples
#' library(ggplot2)
#' fitdata = budyko_errbounds(data=testdata,method="Fu")
#' blankBC+
#'   geom_polygon(data=subset(fitdata,key!="fit"),aes(lty=key),
#'                col=1,fill="transparent")+
#'   geom_line(data=subset(fitdata,key=="fit"),aes(lty="fit"))+
#'   scale_linetype_manual(values=c("fit"=1,"confidence"=2,"prediction"=3))+
#'   geom_point(data=testdata)+
#'   coord_cartesian(xlim=c(0,5),ylim=c(0,1))
#' @md
#' @export

budyko_errbounds=function(fit=NULL,res=0.1,hshift=FALSE,
                          dif="nls",data=NULL,alpha=0.05,method=NULL){
  if(dif=="nls"){
    data=data[,c("AET.P","PET.P")]
    data = subset(data,!is.na(AET.P))
    data = subset(data,PET.P<=20)
    data = subset(data,!(AET.P>PET.P))
    data=data[order(data$PET.P),]
    formula = list(
      "fu"="AET.P~1+PET.P-(1+(PET.P)^p)^(1/p)",
      "turc-pike"="AET.P~(1+(PET.P)^(-p))^(-1/p)",
      "wang-tang"="AET.P~(1+PET.P-((1+PET.P)^2-4*p*(2-p)*PET.P)^(1/2))/(2*p*(2-p))",
      "zhang"="AET.P~(1+p*PET.P)/(1+p*PET.P+(PET.P^(-1)))"
    )[[tolower(method)]]
    startval = list(
      "fu"=list(p=2.7),
      "turc-pike"=list(p=1.9),
      "wang-tang"=list(p=0.59),
      "zhang"=list(p=1.4)
    )[[tolower(method)]]
    budykoNLS = stats::nls(formula = formula,
                           start = startval,
                           data = data)
    fiterrorbounds = do.call("rbind",lapply(1:2,function(i){
      int = c("confidence","prediction")[i]
      dt = data.frame(PET.P=seq(0,10,0.01))
      cd = cbind(dt,investr::predFit(budykoNLS,dt,interval=int,level=(1-alpha)))
      if(i==1){jval=1:3}
      if(i==2){jval=2:3}
      do.call("rbind",lapply(jval,function(j){
        lim = c("fit","lwr","upr")[j]
        cd2 = cd[,c("PET.P",lim)]
        #cd = data.table::setorder(cd,PET.P)
        if(j==3){cd2 = as.data.frame(purrr::map_df(cd2,rev))}
        names(cd2)[2]="AET.P"
        cd2$key = ifelse(j==1,"fit",int)
        return(cd2)
      }))
    }))
  }
  if(dif %in% c("rsq","mae")){
    err=fit[2,]
    if(hshift==TRUE){hs=fit[4,]}else{hs=0}
    fitdata=budyko_sim(fit=fit,method = tolower(names(fit)),hshift = hshift,res = res)
    fiterrorbounds=data.frame(
      PET.P=c(fitdata$PET.P,rev(fitdata$PET.P)),
      AET.P=c((fitdata$AET.P+err),rev(fitdata$AET.P-err))
    )
    #fiterrorbounds$AET.P[fiterrorbounds$AET.P>=1]=1
    #fiterrorbounds$AET.P[fiterrorbounds$AET.P<=0]=0
    for(p in seq(0,1,0.01)){fiterrorbounds$AET.P[fiterrorbounds$PET.P==p & fiterrorbounds$AET.P>=p]=p}
  }
  return(fiterrorbounds)
}

