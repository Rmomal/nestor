#' ggplot heatmap
#'
#' @param data Input matrix.
#' @param no.names Boolean controlling the display of variable names and ticks.
#' @param order Defines a specific display order of the heatmap's rows.
#' @return A heatmap of the input data
#' @export
#' @importFrom reshape2 melt
#' @import ggplot2
#' @examples Sigma=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=FALSE)$Sigma
#' ggimage(Sigma)
ggimage<-function(data, no.names=FALSE, order=NULL){
  num=FALSE
  melted_data <- reshape2::melt(data)
  if(is.null(order)){p=ncol(data) ; order = 1:p}
 if(is.null(colnames(data))){
   level=1:p
   num=TRUE
   }else{level=colnames(data)}
  melted_data$Var1 <- factor( melted_data$Var1, levels = level[order])
  melted_data$Var2 <- factor( melted_data$Var2, levels = level[order])
  if(num){
    melted_data[,1:2]=apply(melted_data[,1:2],2, function(x) as.numeric(as.character(x)))
  }
  g=ggplot2::ggplot(melted_data, ggplot2::aes(x=.data$Var1, y=.data$Var2, fill=.data$value)) + ggplot2::theme_light()+
    ggplot2::labs(x="",y="")+ ggplot2::geom_tile() +ggplot2::guides(fill=FALSE)+
    ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.5))+
    ggplot2::coord_fixed()
  if(no.names) g=g+theme(axis.text=element_blank(),axis.ticks =element_blank())
  g
}

#' wraper of auc from ROCR
#' @param pred Predictions.
#' @param label Labels.
#' @return The Area Under the Curve value.
#' @export
#' @importFrom ROCR prediction performance
#' @examples  data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=FALSE)
#' PLNfit<-norm_PLN(data$Y)
#' MO<-PLNfit$MO
#' SO<-PLNfit$SO
#' sigma_O=PLNfit$sigma_O
#' initClique=data$TC
#' initList=initVEM(data$Y,cliqueList=initClique,sigma_O, MO,r=1 )
#' nestorFit=nestorFit(data$Y, MO,SO, initList=initList, maxIter=3,verbatim=1)
#' res=auc(nestorFit$Pg,data$G)

auc<-function(pred,label){
  prediction<-ROCR::prediction(as.vector(pred[upper.tri(pred)]),as.vector(label[upper.tri(label)]))
  ROC_auc <- round(ROCR::performance(prediction,"auc")@y.values[[1]],digits=2)
  return(ROC_auc)
}
#' Computes precision and recall statistics separating observed from hidden variables, and (FPR,FNR) for hidden variables .
#'
#' @param probs Matrix of estimated edges probabilities.
#' @param G Original graph.
#' @param r Number of missing actors.
#' @param thresh Required threshold for criteria computations, default to 0.5.
#'
#' @return \itemize{
#' \item{PPV:}{ precision of the whole data.}
#' \item{PPVH:}{ precision regarding observed data.}
#' \item{PPVO:}{ precision regarding hidden data.}
#' \item{TPR:}{ recall of the whole data.}
#' \item{TPRH:}{ recall regarding observed data.}
#' \item{TPRO:}{ recall regarding hidden data.}
#' \item{FPRH:}{ false Positive Rate of hidden data.}
#' \item{FNRH:}{ false Negative Rate of hidden data.}
#' }
#' @export
#'
#' @examples data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=FALSE)
#' PLNfit<-norm_PLN(data$Y)
#' MO<-PLNfit$MO
#' SO<-PLNfit$SO
#' sigma_O=PLNfit$sigma_O
#' #-- use true clique for example
#' initClique=data$TC
#' #-- initialize the VEM
#' initList=initVEM(data$Y,cliqueList=initClique,sigma_O, MO,r=1 )
#' nestorFit=nestorFit(data$Y, MO,SO, initList=initList, maxIter=3,verbatim=1 )
#' #-- obtain criteria
#' ppvtpr(nestorFit$Pg,r=1, data$G)
ppvtpr<-function(probs,G,r, thresh=0.5){
  q=ncol(probs) ;
  PPV=round(sum((G!=0)*(probs>thresh))/(sum((G!=0)*(probs>thresh))+ sum((G==0)*(probs>thresh))),2)#TP/(TP+FP)
  TPR=round(sum((G!=0)*(probs>thresh))/sum(G!=0), 2)
  if(r>0){
    h=(q-r):q
    PPVH=round(sum((G[h,]!=0)*(probs[h,]>thresh))/(sum((G[h,]!=0)*(probs[h,]>thresh))+ sum((G[h,]==0)*(probs[h,]>thresh))),2)
    PPVO=round(sum((G[-h,-h]!=0)*(probs[-h,-h]>thresh))/(sum((G[-h,-h]!=0)*(probs[-h,-h]>thresh))+  sum((G[-h,-h]==0)*(probs[-h,-h]>thresh))),2)
    TPRH=round(sum((G[h,]!=0)*(probs[h,]>thresh))/sum(G[h,]!=0), 2)
    TPRO=round(sum((G[-h,-h]!=0)*(probs[-h,-h]>thresh))/sum(G[-h,-h]!=0), 2)
    P=(G[h,]!=0) ; N=(G[h,]==0)
    FP=sum((probs[h,]>thresh)*(G[h,]==0))
    TN=sum((probs[h,]<thresh)*(G[h,]==0))
    FPRH=FP/(FP+TN)
    FNRH=1-TPRH
  }else{
    PPVO=PPV ; TPRO=TPR
    PPVH=TPRH=FPRH=FNRH=NULL
  }
  return(list(PPV=PPV,PPVH=PPVH,PPVO=PPVO,TPR=TPR,TPRH=TPRH,TPRO=TPRO,
              FPRH=FPRH,FNRH=FNRH))
}
#' Comparaitve plot of estimated edge probabilities and original graph.
#'
#' @param P Edges probabilities matrix.
#' @param G Adjacency matrix of the original graph.
#' @param r Number of missing actors.
#' @param thresh Required threshold for criteria computations, default to 0.5.
#' @param no.names Boolean controlling the display of variable names and ticks.
#' @return Two heatmaps with performance criteria computed for the specified threshold as title.
#' @export
#' @examples data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=FALSE)
#' PLNfit<-norm_PLN(data$Y)
#' MO<-PLNfit$MO
#' SO<-PLNfit$SO
#' sigma_O=PLNfit$sigma_O
#' #-- use true clique for example
#' initClique=data$TC
#' #-- initialize the VEM
#' initList=initVEM(data$Y,cliqueList=initClique,sigma_O, MO,r=1 )
#' nestorFit=nestorFit(data$Y, MO,SO, initList=initList, maxIter=3,verbatim=1 )
#' #-- obtain criteria
#' plotPerf(nestorFit$Pg, data$G,r=1)
plotPerf<-function(P,G,r,thresh=0.5, no.names=FALSE){
  # plots heatmaps for the chosen threshold and print verdicts as title
  q=ncol(G)
  criteria=ppvtpr(P,G,r,thresh)
  PPV=criteria$PPV ;PPVH=criteria$PPVH ; PPVO=criteria$PPVO
  TPR=criteria$TPR ;TPRH=criteria$TPRH ; TPRO=criteria$TPRO
  p1<-ggimage(P,no.names=no.names)+labs(title=paste0("G hat"))
  p2<-ggimage(G,no.names=no.names)+labs(title="True G")
  auc<-round(auc(pred = P, label = G),3)
  if(r>0){
    grid.arrange(p1,p2,ncol=2, top=paste0("Recall=",TPR," (Obs=",TPRO," , Hid=",TPRH,
                                          ")\n Precision=",PPV," (Obs=",PPVO," , Hid=",PPVH,")",
                                          "\n AUC=",auc))
  }else{
    grid.arrange(p1,p2,ncol=2, top=paste0("Recall=",TPR,", Precision=",PPV,",\nAUC=",auc))
  }

}

#' Plot function for nestorFit convergence
#'
#' @param nestorFit  Fit from the nestorFit() function.
#' @return Visualiaztion of nestorFit convergence.
#' @export
#' @import ggplot2
#' @importFrom tidyr gather
#' @importFrom dplyr rename
#' @importFrom tibble rowid_to_column
#' @importFrom gridExtra grid.arrange
#' @examples data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=FALSE)
#' PLNfit<-norm_PLN(data$Y)
#' MO<-PLNfit$MO
#' SO<-PLNfit$SO
#' sigma_O=PLNfit$sigma_O
#' #-- use true clique for example
#' initClique=data$TC
#' #-- initialize the VEM
#' initList=initVEM(data$Y,cliqueList=initClique,sigma_O, MO,r=1 )
#' nestorFit=nestorFit(data$Y, MO,SO, initList=initList, maxIter=3 , verbatim=1)
#' #-- obtain criteria
#' plotConv(nestorFit)
plotConv<-function(nestorFit){
  trackJ=(nrow(nestorFit$lowbound)>(2*nestorFit$finalIter))
  mytheme.dark <-function(legend){
    list= list(theme_light(), scale_color_brewer(legend,palette="Dark2"),
               scale_fill_brewer(legend,palette="Dark2"),
               theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                     plot.title = element_text(hjust = 0.5)))
    return(list)}

  mytheme <- list(theme_light(), scale_fill_brewer("",palette="Dark2"),#scale_colour_hp_d(option = "RavenClaw", name = ""), #scale_color_manual("",values=mypal),
                  theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                        plot.title = element_text(hjust = 0.5)))

    g1<-nestorFit$features  %>% dplyr::select(.data$diffPg, .data$diffOmega) %>%
      dplyr::rename(Prob=.data$diffPg, Omega=.data$diffOmega) %>%
      tibble::rowid_to_column() %>%
      tidyr::gather(key, values, -.data$rowid) %>%
      ggplot2::ggplot(aes(.data$rowid,.data$values, color=.data$key))+ geom_point()+geom_line() + facet_wrap(~.data$key, scales="free")+
      labs(x="",y="", title="Parameters")+ mytheme.dark("")+guides(color=FALSE)
    if(trackJ){
      g2<- nestorFit$lowbound %>% tibble::rowid_to_column() %>%  tidyr::gather(key,value,-.data$rowid,-.data$parameter) %>%
        ggplot2::ggplot(aes(.data$rowid,.data$value, group=.data$key))+geom_line()+geom_point(aes(color=as.factor(.data$parameter)), size=2, alpha=0.8)+
        facet_wrap(~key, scales="free")+  labs(x="sub-iteration",y="", title="Lower bound and components")+mytheme.dark("")
    }else{ g2<- nestorFit$lowbound %>% rowid_to_column() %>%
      ggplot2::ggplot(aes(.data$rowid,.data$J ))+geom_line()+geom_point(aes(color=as.factor(.data$parameter)),size=2, alpha=0.8)+
      labs(x="iteration",y="", title="Lower bound")+mytheme+ scale_color_manual("",values="#2976d6")+
      guides(color=FALSE)}
    gridExtra::grid.arrange(g1,g2, ncol=1)

}
