library(ggrain)
library(patchwork)
library(RColorBrewer)

plotRain <-function(effects4couple,start,end,intvar,outvar, lim){
  intlevels<-nrow(effects4couple[[1]][,,1])
  outlevels<- length(effects4couple[[1]][1,1,])
  bootnumb<- length(effects4couple)
  effsarray <-  array(as.numeric(unlist(effects4couple)), dim=c(intlevels, intlevels, outlevels, bootnumb))
  df<-data.frame(effect=c(),lev_out=c()) 
  meanlev<-c()
  for (l in 1:outlevels){
    effect<-effsarray[start,end,l,]
    meanlev<-c(meanlev,mean(effect))
    lev_out<-factor(rep(l,length(effect)))
    df<-rbind(df,data.frame(effect,lev_out))
  }
  col<-brewer.pal(n = outlevels, name = "Dark2")
  p<-ggplot(df, aes(1,effect, fill = lev_out, color =lev_out)) +
    geom_rain(alpha = .5,rain.side = 'l',  boxplot.args = list(
      color="black",show.legend=TRUE, outlier.shape=NA), 
      boxplot.args.pos = list(position = ggpp::position_dodgenudge(x = .1), width = 0.1))+ 
    theme(panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),axis.title.x = element_blank(), 
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      panel.background = element_rect("white","white")) + 
    labs(y="OCE", fill = "Level outcome variable")+
    scale_fill_brewer(palette = 'Dark2') +
    scale_color_brewer(palette = 'Dark2')+ 
    guides(color="none", fill = "none")+ggtitle(paste0("Int ",intvar), subtitle=paste0(" Out ", outvar))+
    coord_cartesian(clip = "on", ylim= c(-lim,lim))+scale_y_continuous(breaks=seq(-lim,lim,0.2))
  for (i in 1:outlevels){
    p<-p+geom_segment(y=meanlev[i], yend=meanlev[i],x=0.85,xend=0.6, linewidth=0.8, color=col[i])
  }
  return(p)
}
plotRain3 <-function(effects4couple,start,end,intvar,outvar, lim){
  intlevels<-nrow(effects4couple[[1]][,,1])
  outlevels<- length(effects4couple[[1]][1,1,])
  bootnumb<- length(effects4couple)
  effsarray <-  array(as.numeric(unlist(effects4couple)), dim=c(intlevels, intlevels, outlevels, bootnumb))
  df<-data.frame(effect=c(),lev_out=c()) 
  meanlev<-c()
  for (l in 1:outlevels){
    effect<-effsarray[start,end,l,]
    meanlev<-c(meanlev,mean(effect))
    lev_out<-factor(rep(l,length(effect)))
    df<-rbind(df,data.frame(effect,lev_out))
  }
  col<-brewer.pal(n = outlevels, name = "Dark2")
  aes(x=lev_out,y=effect, fill = lev_out, color =lev_out)
  p<-ggplot(df, aes(x=1,y=effect, fill = lev_out, color =lev_out)) +
    geom_rain(alpha = .5,rain.side = 'r',  boxplot.args = list(
      color="black",show.legend=TRUE, outlier.shape=NA), 
      boxplot.args.pos = list(position = ggpp::position_dodgenudge(x = .1), width = 0.1))+ 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x=element_text(),
          panel.background = element_rect("white","white")) + 
    labs(y="OCE",x="",fill = "Level outcome variable")+
    scale_fill_brewer(palette = 'Dark2') +
    scale_color_brewer(palette = 'Dark2')+ 
    guides(color="none", fill = "none")+ggtitle(paste0("Int ",intvar, " from ", start, " to ", end), subtitle=paste0(" Out ", outvar))+
    coord_flip(clip = "on", ylim= c(-lim,lim))
  for (i in 1:outlevels){
    p<-p+geom_segment(x=1.15, xend=1.4,y=meanlev[i],yend=meanlev[i], linewidth=0.8, color=col[i])
  }
  return(p)
}

plotRain2 <-function(effects4couple,start,end,intvar,outvar, lim){
  intlevels<-nrow(effects4couple[[1]][,,1])
  outlevels<- length(effects4couple[[1]][1,1,])
  bootnumb<- length(effects4couple)
  effsarray <-  array(as.numeric(unlist(effects4couple)), dim=c(intlevels, intlevels, outlevels, bootnumb))
  df<-data.frame(effect=c(),lev_out=c()) 
  meanlev<-c()
  for (l in 1:outlevels){
    effect<-effsarray[start,end,l,]
    meanlev<-c(meanlev,mean(effect))
    lev_out<-factor(rep(l,length(effect)))
    df<-rbind(df,data.frame(effect,lev_out))
  }
  col<-brewer.pal(n = outlevels, name = "Dark2")
  p<-ggplot(df, aes(x=lev_out,y=effect, fill = lev_out, color =lev_out)) +
    geom_rain(alpha = .5,rain.side = 'r',  boxplot.args = list(
      color="black",show.legend=TRUE, outlier.shape=NA), 
      boxplot.args.pos = list(position = ggpp::position_dodgenudge(x = .2), width = 0.1),
      point.args.pos = rlang::list2(position = position_jitter(width = 0.04, height = 0)),violin.args.pos = list(
        side = "r",
        width = 0.7, position = position_nudge(x = 0.3)))+ 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #axis.title.x = element_blank(), 
          # axis.text.x = element_blank(), 
          # axis.ticks.x = element_blank(), 
          panel.background = element_rect("white","white")) + 
    labs(y="OCE", x="Outcome level", fill = "Level outcome variable")+
    scale_fill_brewer(palette = 'Dark2') +
    scale_color_brewer(palette = 'Dark2')+ 
    guides(color="none", fill = "none")+ggtitle(paste0("Int ",intvar), subtitle=paste0(" Out ", outvar))+
    coord_flip(clip = "on", ylim= c(-lim,lim))#+scale_y_continuous(breaks=seq(-lim,lim,0.1))
  for (i in 1:outlevels){
    p<-p+geom_segment(y=meanlev[i], yend=meanlev[i],x=i+0.1,xend=i+0.6, linewidth=0.5, color=col[i])
  }
  return(p)
}

plotRain4 <-function(effects4couple,start,end,intvar,outvar, lim){
  intlevels<-nrow(effects4couple[[1]][,,1])
  outlevels<- length(effects4couple[[1]][1,1,])
  bootnumb<- length(effects4couple)
  effsarray <-  array(as.numeric(unlist(effects4couple)), dim=c(intlevels, intlevels, outlevels, bootnumb))
  df<-data.frame(effect=c(),lev_out=c()) 
  meanlev<-c()
  l <-3
    effect<-effsarray[start,end,l,]
    meanlev<-c(meanlev,mean(effect))
    lev_out<-factor(rep(l,length(effect)))
    df<-rbind(df,data.frame(effect,lev_out))
  col<-"#D95F02"
  p<-ggplot(df, aes(x=effect))+
    geom_density(alpha=.5, fill="#D95F02")+
    #ggplot(df, aes(x=-0.3,y=effect, fill = "#D95F02", color ="#D95F02")) +
    # geom_rain(alpha = .5,rain.side = 'r',  boxplot.args = list(
    #   color="black",show.legend=TRUE, outlier.shape=NA), 
    #   boxplot.args.pos = list(position = ggpp::position_dodgenudge(x = .2), width = 0.1),
    #   point.args.pos = rlang::list2(position = position_jitter(width = 0.04, height = 0)),violin.args.pos = list(
    #     side = "r",
    #     width = 0.7, position = position_nudge(x = 0.3)))+ 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          panel.background = element_rect("white","white"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) + 
    labs(x="OCE")+
    #scale_fill_brewer(palette = 'Dark2') +
    #scale_color_brewer(palette = 'Dark2')+ 
    guides(color="none", fill = "none")+ggtitle(paste0("Int ",intvar, " from ", start, " to ", end), subtitle=paste0(" Out ", outvar, " Level 3 "))
    #coord_cartesian(clip = "on", xlim= c(-lim,lim))#+scale_y_continuous(breaks=seq(-lim,lim,0.1))
    p<-p+geom_segment(x=meanlev, xend=meanlev,y=-0.0005,yend=60, linewidth=0.5, color="#D95F02")
  return(p)
}

plotRain5 <-function(effects4couple,start,end,intvar,outvar, lim){
  intlevels<-nrow(effects4couple[[1]][,,1])
  outlevels<- length(effects4couple[[1]][1,1,])
  bootnumb<- length(effects4couple)
  effsarray <-  array(as.numeric(unlist(effects4couple)), dim=c(intlevels, intlevels, outlevels, bootnumb))
  df<-data.frame(effect=c(),lev_out=c()) 
  meanlev<-c()
  l <- round(outlevels/2)
  effect<-effsarray[start,end,l,]
  meanlev<-c(meanlev,mean(effect))
  lev_out<-factor(rep(l,length(effect)))
  df<-rbind(df,data.frame(effect,lev_out))
  col<-brewer.pal(n = outlevels, name = "Dark2")
  p<-ggplot(df, aes(x=effect))+
    geom_density(alpha=.5, fill=col[l])+
    #ggplot(df, aes(x=-0.3,y=effect, fill = "#D95F02", color ="#D95F02")) +
    # geom_rain(alpha = .5,rain.side = 'r',  boxplot.args = list(
    #   color="black",show.legend=TRUE, outlier.shape=NA), 
    #   boxplot.args.pos = list(position = ggpp::position_dodgenudge(x = .2), width = 0.1),
    #   point.args.pos = rlang::list2(position = position_jitter(width = 0.04, height = 0)),violin.args.pos = list(
    #     side = "r",
    #     width = 0.7, position = position_nudge(x = 0.3)))+ 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          panel.background = element_rect("white","white"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) + 
    labs(x="OCE")+
    #scale_fill_brewer(palette = 'Dark2') +
    #scale_color_brewer(palette = 'Dark2')+ 
    guides(color="none", fill = "none")+ggtitle(paste0("Int ",intvar, " from ", start, " to ", end), subtitle=paste0(" Out ", outvar, " Level ", l))+
  coord_cartesian(xlim= lim)#+scale_y_continuous(breaks=seq(-lim,lim,0.1))
  p<-p+geom_segment(x=meanlev, xend=meanlev,y=-0.0005,yend=80, linewidth=0.5, color=col[l])
  return(p)
}


plotRaindouble <-function(effects4couple1,effects4couple2,start,end,intvar,outvar, lim){
  intlevels<-nrow(effects4couple1[[1]][,,1])
  outlevels<- length(effects4couple1[[1]][1,1,])
  bootnumb<- length(effects4couple1)
  effsarray1 <-  array(as.numeric(unlist(effects4couple1)), dim=c(intlevels, intlevels, outlevels, bootnumb))
  effsarray2 <-  array(as.numeric(unlist(effects4couple2)), dim=c(intlevels, intlevels, outlevels, bootnumb))
  df1<-data.frame(effect=c(),lev_out=c(), method=c()) 
  df2<-data.frame(effect=c(),lev_out=c(), method=c()) 
  meanlev1<-c()
  meanlev2<-c()
  for (l in 1:outlevels){
    effect1<-effsarray1[start,end,l,]
    effect2<-effsarray2[start,end,l,]
    meanlev1<-c(meanlev1,mean(effect1))
    meanlev2<-c(meanlev2,mean(effect2))
    lev_out1<-factor(rep(l,length(effect1)))
    lev_out2<-factor(rep(l,length(effect2)))
    method1 = factor(rep("OSEM",length(effect1)))
    method2 = factor(rep("Param",length(effect2)))
    df1<-rbind(df1,data.frame(effect=effect1,lev_out=lev_out1, method=method1))
    df2<-rbind(df2,data.frame(effect=effect2,lev_out=lev_out2, method=method2))
  }
  df<-rbind(df1,df2)
  #col<-c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02" ,"#A6761D","#666666")
  col <- rev(c(
    rgb(31, 120, 180, maxColorValue = 255, alpha = 150),  # Blue
    rgb(51, 160, 44, maxColorValue = 255, alpha = 150),   # Green
    rgb(227, 26, 28, maxColorValue = 255, alpha = 150),   # Red
    rgb(255, 127, 0, maxColorValue = 255, alpha = 150),   # Orange
    rgb(106, 61, 154, maxColorValue = 255, alpha = 150),  # Purple
    rgb(166, 206, 227, maxColorValue = 255, alpha = 150), # Light Blue
    rgb(178, 223, 138, maxColorValue = 255, alpha = 150), # Light Green
    rgb(253, 191, 111, maxColorValue = 255, alpha = 150)  # Light Orange
  ))
  p<-ggplot(df, aes(x=method,y=effect, fill = lev_out, color =lev_out)) +
    geom_rain(alpha = .5,rain.side = 'r'#,  boxplot.args = list(
      # color="black",show.legend=TRUE, outlier.shape=NA), 
      # boxplot.args.pos = list(position = ggpp::position_dodgenudge(x = .1), width = 0.1)
      )+ 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect("white","white")) + 
    labs(y="OCE", fill = "Level outcome variable")+
    scale_colour_manual(values = col[1:outlevels])+
    scale_fill_manual(values = col[1:outlevels])+
    guides(color="none", fill ="none")+ggtitle(paste0("Int ",intvar, " from ", start, " to ", end), subtitle=paste0(" Out ", outvar))+
    coord_flip(clip = "on", ylim= c(-lim,lim))#+scale_y_continuous(breaks=seq(-lim,lim,0.2))
  for (i in 1:outlevels){
    p<-p+geom_segment(y=meanlev1[i], yend=meanlev1[i],x=1.15,xend=1.55, linewidth=0.8, color=col[i])+geom_segment(y=meanlev2[i], yend=meanlev2[i],x=2.15,xend=2.55, linewidth=0.8, color=col[i])
  }
  return(p)
}


# Function to extract legend
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend) }

# To Plot separately OSEM and Param simulations 
PlotTot<- function(int,lim,n){
  plot<-list()
  for (j in 1:n){
    if (j == int){
      blank<-data.frame()
      plot[[j]]<-ggplot(blank) + geom_blank()+
        theme(panel.background = element_rect(fill="white","white"),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+ggtitle(paste0("Int ",int), subtitle=paste0(" Out ", int))
    } else{
      start = 1
      end = nrow(effects4couple[[int]][[j]][[1]][,,1])
      plot[[j]]<-plotRain2(effects4couple=effects4couple[[int]][[j]], start, end, intvar=int, outvar=j, lim = lim)
    }
  }
  return(plot)
}

# Generate plots for plot matrix 24x24 of effects on outcome variable 
# of shifting intervention from its lowest level (start = 1) to its highest (end = nrow(effects4couple[[int]][[j]][[1]][,,1]))
PlotTotpres<- function(int,lim){
  plot<-list()
  for (j in 1:n){
    if (j == int){
      blank<-data.frame()
      plot[[j]]<-ggplot(blank) + geom_blank()+
        theme(panel.background = element_rect(fill="white","white"),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+ggtitle(paste0("Int ",int), subtitle=paste0(" Out ", int))
    } else{
      start = 1
      end = nrow(effects4couple[[int]][[j]][[1]][,,1])
      plot[[j]]<-plotRain5(effects4couple=effects4couple[[int]][[j]], start, end, intvar=int, outvar=j, lim = lim)
    }
  }
  return(plot)
}

PlotTotdouble<- function(int,lim){
  plot<-list()
  for (j in 1:5){ #should be n
    if (j == int){
      blank<-data.frame()
      plot[[j]]<-ggplot(blank) + geom_blank()+
        theme(panel.background = element_rect(fill="white","white"),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+ggtitle(paste0("Int ",int), subtitle=paste0(" Out ", int))
    } else{
      start = 1
      end = nrow(effects4couple_OSEM[[int]][[j]][[1]][,,1])
      plot[[j]]<-plotRaindouble(effects4couple1=effects4couple_OSEM[[int]][[j]],effects4couple2=effects4couple_param[[int]][[j]], start, end, intvar=int, outvar=j, lim = lim)
    }
  }
  return(plot)
}

# Plot the effects of any level-switch of the intervention variable on all levels of the outcome variable
PlotLocal<- function(effects4couple, intvar, outvar, lim){
  count<-1
  plot<-list()
  intlevels<-nrow(effects4couple[[1]][,,1])
  for (ii in c(1:(intlevels-1))){
    for (jj in c((ii+1):intlevels)){
      start = ii
      end = jj
      plot[[count]]<-plotRain2(effects4couple=effects4couple, start, end, intvar, outvar, lim = lim)+ggtitle(paste0("Int ", intvar, " from ",ii, " to ", jj))
      count<-count+1
      }
  }
  return(plot)
}

PlotTotdoublemean<- function(int,lim){
  plot<-list()
  for (j in 1:5){
    col<-c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02" ,"#A6761D","#666666")
    if (j == int){
      blank<-data.frame()
      plot[[j]]<-ggplot(blank) + geom_blank()+
        theme(panel.background = element_rect(fill="white","white"),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+ggtitle(paste0("Int ",int), subtitle=paste0(" Out ", int))
    } else{
      start = 1
      end = nrow(effects4couple_OSEM[[int]][[j]][[1]][,,1])
      levout<-length(effects4couple_OSEM[[int]][[j]][[1]][1,1,])
      plot[[j]]<-plotRaindouble(effects4couple1=effects4couple_OSEM[[int]][[j]],effects4couple2=effects4couple_param[[int]][[j]], start, end, intvar=int, outvar=j, lim = lim)
      for (lev in 1:levout){
        plot[[j]]<-plot[[j]]+
          geom_segment(y=Effects_Owen[[int]][[j]][1,2,lev], yend=Effects_Owen[[int]][[j]][1,2,lev],x=1.15,xend=2.55, linewidth=0.8, color=col[lev],linetype="dashed")
      }
    }
  }
  return(plot)
}
