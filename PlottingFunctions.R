# Venom Transcriptome Plotting Functions & Colors
# Functions by: Andrew J. Mason & Rhett M. Rautsaw

library(RColorBrewer)
library(ggpubr)

toxin_colors<-c("#d30b94","3FTx",
                "#201923","BPP",
                "#772b9d","CRISP",
                "#0ec434","CTL",
                "#a4a4a4","Ficolin",
                "#ffffff","FusedToxin",
                "#29bdab","HYAL",
                "#ffdab9","KUN",#e68f66
                "#ffc413","LAAO",
                "#fcff5d","Lacta",
                "#f47a22","MYO",
                "#fffac8","NGF",#ffcba5
                "#5d4c86","NUC",
                "#7cd676","PDE",
                "#8B5A2B","PDI",
                "#f52727","PLA2",
                "#991919","PLA2_neuro",
                "#c3a5b4","PLB",
                "#3998f5","SVMMP",
                "#277da7","SVMPI",
                "#3750db","SVMPII",
                "#2f2aa0","SVMPIII",
                "#228c68","SVSP",
                "#aaffc3","VEGF",#946aa2
                "#aa6e28","Vespryn",#c56133
                "#cccccc","VF",
                "#f07cab","Waprin")
#                "#235b54","zOpenColor1",
#                "#37294f","zOpenColor2",
#                "#96341c","zOpenColor3",
#                "#632819","zOpenColor4",
#                "#b732cc","zOpenColor5",
#                "#8ad8e8","zOpenColor6")

toxin_colors_df<-matrix(toxin_colors,nrow=length(toxin_colors)/2,ncol=2,byrow=T)
toxin_colors_df<-as.data.frame(cbind(toxin_colors_df,rep(1,nrow(toxin_colors_df))))
toxin_colors_df$V2 <- factor(toxin_colors_df$V2, levels = toxin_colors_df$V2)
toxin_colors_df$V3 <- as.numeric(toxin_colors_df$V3)
#ggbarplot(toxin_colors_df, "V2","V3", fill=toxin_colors_df$V1, width = 1, xlab="",ylab="", main="toxin colors") + rotate_x_text(angle = 45)
toxin_colors<-palette(as.vector(toxin_colors_df$V1))
toxin_colors<-palette(as.vector(toxin_colors_df$V1))
names(toxin_colors)<-toxin_colors_df$V2

ToxinNon_colors <- colorRampPalette(c('black','gray80'))

## Full transcriptome barplot
FullTranscriptomePlot<-function(df=TPM_df2,id="Average",class="class",print=TRUE){
  df[[id]]<-log(df[[id]])
  transcript<-names(df)[1]
  
  df[[class]]<-factor(df[[class]],levels=c("Toxin","Nontoxin",levels(df[[class]])[grep("Toxin|Nontoxin",levels(df[[class]]),invert=T)]))
  
  A<-ggbarplot(df, transcript, id, sort.val = "desc", sort.by.groups = F,
                size=0, width=1, fill=class, palette = ToxinNon_colors(nlevels(df[[class]])),legend="right",
                ylab="ln(TPM)", title=paste(id,"Venom Gland Transcriptome")) + 
    rremove("xlab")+rremove("x.axis")+rremove("x.text")+rremove("x.ticks")

  if(print==TRUE){
    print(A)
  }
  if(print==FALSE){
    return(A)
  }
}

## Toxin barplot
ToxinBarplot<-function(df=TPM_df2,id="Average",class="class",toxin_class="toxin_class",colors=toxin_colors,print=TRUE){
  df[[id]]<-log(df[[id]])
  df<-df[df[class]=="Toxin",]
  #df<-subset(df,df[[class]]=="Toxin")
  transcript<-names(df)[1]
  labels <- df[[transcript]]
  labels<-sub("^.*[^*]","",labels)
  labels<-sub("...","*",labels)
  
  A<-ggbarplot(df, transcript, id, sort.val = "desc", sort.by.groups = F,
               size=0.5, width=1, fill=toxin_class, palette = colors,legend="right",
               ylab="ln(TPM)", title=paste(id,"Venom Gland Transcriptome"), label=labels) + 
               rremove("xlab")+rremove("x.axis")+rremove("x.text")+rremove("x.ticks")
  if(print==TRUE){
    print(A)
  }
  if(print==FALSE){
    return(A)
  }
}

## Expression piecharts
ExpressionPie<-function(df=TPM_df2,id="Average",class="class",toxin_class="toxin_class",colors=toxin_colors,print=TRUE){
  tmp_sum<-aggregate(df[[id]],by=list(class=df[[class]]),FUN=sum)
  tmp_sum$label <- paste0(tmp_sum$class," ",round((tmp_sum$x/sum(tmp_sum$x))*100,2),"%")
  
  tmp_sum[[class]]<-factor(tmp_sum[[class]],levels=c("Toxin","Nontoxin",levels(tmp_sum[[class]])[grep("Toxin|Nontoxin",levels(tmp_sum[[class]]),invert=T)]))
  
  A<-ggpie(tmp_sum,"x",fill = "class", palette=ToxinNon_colors(nlevels(as.factor(df[[class]]))), label="label",lab.pos = "in", lab.font = "white")+rremove("legend")

  df2<-df[df[class]=="Toxin",]
  #df2<-subset(df,df[[class]]=="Toxin")
  tmp_sum<-aggregate(df2[[id]],by=list(class=df2[[toxin_class]]),FUN=sum)
  tmp_sum$label <- paste0(tmp_sum$class," ",round((tmp_sum$x/sum(tmp_sum$x))*100,2),"%")
  tmp_sum<-tmp_sum[order(tmp_sum$x),]
  tmp_sum$class <- factor(tmp_sum$class, levels = tmp_sum$class)
  B<-ggpie(tmp_sum,"x",fill = "class", palette=colors, label="label",lab.pos = "in", lab.font = "white")+rremove("legend")
  
  C<-A + B
  if(print==TRUE){
    print(C)
  }
  if(print==FALSE){
    return(C)
  }
}

## FancyFigure
FancyFigure <- function(df=TPM_df2,id="Average",class="class",toxin_class="toxin_class",colors=toxin_colors){
  A <- FullTranscriptomePlot(df,id,class,print=FALSE)
  B <- ToxinBarplot(df,id,class,toxin_class,colors,print=FALSE)
  C <- ExpressionPie(df,id,class,toxin_class,colors,print=FALSE)
  
  D<-(B + C) / A
  print(D)
}

## Pairwise scatterplot of clr transformed expression data
TransCompPlot<-function(dat,x,y){
  clr_dat <- data.frame(cbind(dat[,1],dat[,x],dat[,y]))
  for (i in 2:3){
    clr_dat[,i] <- as.numeric(as.character(clr_dat[,i]))}
  rownames(clr_dat) = clr_dat[,1 ] # the first row will be the header
  clr_dat = clr_dat[,-1]
  colnames(clr_dat) = c(colnames(dat[x]),colnames(dat[y]))
  clr_dat <- t(clr_dat)
  clr_dat <- clr(clr_dat)
  clr_dat <- as.data.frame(t(clr_dat))
  clr_dat <-cbind(dat[,2:3],clr_dat)
  
  Nontoxin<-subset(clr_dat,clr_dat[[1]]=="Nontoxin")
  Toxins<-subset(clr_dat,clr_dat[[1]]=="Toxin")
  
  varcovar<-matrix(c(var(Nontoxin[[3]]),cov(Nontoxin[[3]],Nontoxin[[4]]),
                     cov(Nontoxin[[4]],Nontoxin[[3]]),var(Nontoxin[[4]])),2,2)
  e1<-eigen(varcovar)
  s<-(2.58*(sqrt(e1$values[2])))
  lm1<-lm(Nontoxin[[4]]~Nontoxin[[3]])
  lm2<-lm(Toxins[[4]]~Toxins[[3]])
  
  R2_Non<-summary(lm1)$r.squared
  R_non<-cor(Nontoxin[[3]],Nontoxin[[4]],method="pearson")
  p_non<-cor(Nontoxin[[3]],Nontoxin[[4]],method="spearman")
  
  R2_Tox<-summary(lm2)$r.squared
  R_Tox<-cor(Toxins[[3]],Toxins[[4]],method="pearson")
  p_Tox<-cor(Toxins[[3]],Toxins[[4]],method="spearman")
  
  print("Nontoxin")
  print(cat("n=",(length(Nontoxin[[3]]))))
  print(cat("p=",p_non))
  print(cat("R=",R_non))
  print(cat("R2=",R2_Non))
  print(" ")
  print("Toxins")
  print(cat("n=",(length(Toxins[[3]]))))
  print(cat("p=",p_Tox))
  print(cat("R=",R_Tox))
  print(cat("R2=",R2_Tox))
  
  MinX=min(c(min(Nontoxin[[3]]),min(Toxins[[3]])))
  MinY=min(c(min(Nontoxin[[4]]),min(Toxins[[4]])))
  
  
  ggplot(data=clr_dat, aes(x=clr_dat[[3]], y=clr_dat[[4]]))+
    geom_point(data=Nontoxin,aes(x=Nontoxin[[3]],y=Nontoxin[[4]]),color="grey63",size=2)+
    #geom_abline(slope=1,intercept = 0,color="black",size=1)+
    geom_abline(slope=(lm1$coefficients[2]),intercept=lm1$coefficients[1],color="red",size=1)+
    geom_segment(x=1.1*MinX,y=1.1*(lm1$coefficient[2]*MinX)+((lm1$coefficients[1])+s),
                 yend=(1.1*max(clr_dat[[4]])), 
                 xend=((1.1*max(clr_dat[[4]]))-(lm1$coefficient[1]+s))/(lm1$coefficients[2]), 
                 size=1,linetype="dashed",alpha=0.1)+
    geom_segment(aes(x=(1.15*(MinY-(lm1$coefficients[1]-s))/(lm1$coefficients[2])),
                     y=1.15*MinY,xend=(1.1*max(clr_dat[[3]])),
                     yend=(lm1$coefficients[2]*max(1.1*clr_dat[[3]]))+(lm1$coefficients[1]-s)),
                 color="black",size=1,linetype="dashed",alpha=0.1)+
    # geom_point(data=Toxins,aes(x=Toxins[[3]],y=Toxins[[4]],color=factor(Toxins[[2]])),
    #            size=3)+
    geom_point(data=Toxins,aes(x=Toxins[[3]],y=Toxins[[4]],fill=factor(Toxins[[2]])),
               size=3,pch=21,colour="black")+
    scale_fill_manual(values=toxin_colors, name="Toxin Class") +
    #geom_text(data=Toxins,aes(x=Toxins[[3]],y=Toxins[[4]],label=rownames(Toxins)))+
    labs(x=colnames(Toxins[3]),y=colnames(Toxins[4]))+
    theme(legend.background=element_rect(colour="black"),legend.key=element_rect(colour="black",fill="white"),
          axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
          axis.title=element_text(size=14),
          panel.background = element_rect(fill = "white"),line=element_line(colour="black",size=1),
          plot.title=element_text(hjust=0.5,size=20))+
    coord_fixed()
}


