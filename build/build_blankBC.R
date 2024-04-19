blankBC=ggplot2::ggplot(data=data.frame(PET.P=seq(0,20,1),AET.P=c(0,rep(1,20))),ggplot2::aes(x=PET.P,y=AET.P))+
  ggplot2::geom_line()+
  ggplot2::xlab("Aridity Index")+
  ggplot2::ylab("Evaporative Index")
usethis::use_data(blankBC,overwrite = TRUE)