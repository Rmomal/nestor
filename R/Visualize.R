#' ggplot heatmap
#'
#' @param data input matrix
#'
#' @return a heatmap
#' @export
#' @importFrom reshape2 melt
#' @import ggplot2
#' @examples Sigma=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=FALSE)$Sigma
#' ggimage(Sigma)
ggimage<-function(data){
  melted_data <- reshape2::melt(data)
  ggplot2::ggplot(melted_data, aes(x=.data$Var1, y=.data$Var2, fill=.data$value)) + theme_light()+labs(x="",y="")+
    geom_tile() +guides(fill=FALSE)+ theme(plot.title = element_text(size=10, hjust=0.5))+ coord_fixed()
}
#' wraper of auc from ROCR
#'
#' @param pred predictions
#' @param label labels
#'
#' @return The Area Under the Curve value
#' @export
#' @importFrom ROCR prediction performance
#' @examples data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=FALSE)
#' PLNfit<-norm_PLN(data$Y)
#' MO<-PLNfit$MO
#' SO<-PLNfit$SO
#' sigma_O=PLNfit$sigma_O
#' #-- use true clique for example
#' initClique=data$TC
#' #-- initialize the VEM
#' initList=initVEM(data$Y,cliqueList=initClique,sigma_O, MO,r=1 )
#' nestorFit=nestor(data$Y, MO,SO, initList=initList, maxIter=3,verbatim=1)
#' #-- obtain criteria
#' auc(nestorFit$Pg,data$G)
auc<-function(pred,label){
  prediction<-ROCR::prediction(as.vector(pred[upper.tri(pred)]),as.vector(label[upper.tri(label)]))
  ROC_auc <- round(ROCR::performance(prediction,"auc")@y.values[[1]],digits=2)
  return(ROC_auc)
}
#' Computes precision and recall statistics separating observed from hidden variables, and (FPR,FNR) for hidden variables .
#'
#' @param probs matrix of estimated edges probabilities
#' @param G original graph
#' @param r number of missing actors
#' @param thresh required threshold for criteria computations, default to 0.5
#'
#' @return \itemize{
#' \item{PPV}{ Precision of the whole data}
#' \item{PPVH}{ Precision regarding observed data}
#' \item{PPVO}{ Precision regarding hidden data}
#' \item{TPR}{ Recall of the whole data}
#' \item{TPRH}{ Recall regarding observed data}
#' \item{TPRO}{ Recall regarding hidden data}
#' \item{FPRH}{ False Positive Rate of hidden data}
#' \item{FNRH}{ False Negative Rate of hidden data}
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
#' nestorFit=nestor(data$Y, MO,SO, initList=initList, maxIter=3,verbatim=1 )
#' #-- obtain criteria
#' ppvtpr(nestorFit$Pg,r=1, data$G)
ppvtpr<-function(probs,G,r, thresh=0.5){
  q=ncol(probs) ; h=(q-r):q
  PPV=round(sum((G!=0)*(probs>thresh))/(sum((G!=0)*(probs>thresh))+ sum((G==0)*(probs>thresh))),2)#TP/(TP+FP)
  PPVH=round(sum((G[h,]!=0)*(probs[h,]>thresh))/(sum((G[h,]!=0)*(probs[h,]>thresh))+ sum((G[h,]==0)*(probs[h,]>thresh))),2)
  PPVO=round(sum((G[-h,-h]!=0)*(probs[-h,-h]>thresh))/(sum((G[-h,-h]!=0)*(probs[-h,-h]>thresh))+  sum((G[-h,-h]==0)*(probs[-h,-h]>thresh))),2)
  TPR=round(sum((G!=0)*(probs>thresh))/sum(G!=0), 2)
  TPRH=round(sum((G[h,]!=0)*(probs[h,]>thresh))/sum(G[h,]!=0), 2)
  TPRO=round(sum((G[-h,-h]!=0)*(probs[-h,-h]>thresh))/sum(G[-h,-h]!=0), 2)
  P=(G[h,]!=0) ; N=(G[h,]==0)
  FP=sum((probs[h,]>thresh)*(G[h,]==0))
  TN=sum((probs[h,]<thresh)*(G[h,]==0))
  FPRH=FP/(FP+TN)
  FNRH=1-TPRH
  return(list(PPV=PPV,PPVH=PPVH,PPVO=PPVO,TPR=TPR,TPRH=TPRH,TPRO=TPRO,
              FPRH=FPRH,FNRH=FNRH))
}
#' Comparaitve plot of estimated edge probabilities and original graph.
#'
#' @param P edges probabilities matrix
#' @param G adjacency matrix of the original graph
#' @param r number of missing actors
#' @param thresh required threshold for criteria computations, default to 0.5
#'
#' @return two heatmaps with performance criteria computed for the specified threshold as title.
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
#' nestorFit=nestor(data$Y, MO,SO, initList=initList, maxIter=3,verbatim=1 )
#' #-- obtain criteria
#' plotPerf(nestorFit$Pg, data$G,r=1)
plotPerf<-function(P,G,r,thresh=0.5){
  # plots heatmaps for the chosen threshold and print verdicts as title
  q=ncol(G)
  criteria=ppvtpr(P,G,r,thresh)
  PPV=criteria$PPV ;PPVH=criteria$PPVH ; PPVO=criteria$PPVO
  TPR=criteria$TPR ;TPRH=criteria$TPRH ; TPRO=criteria$TPRO
  p1<-ggimage(P)+labs(title=paste0("G hat"))
  p2<-ggimage(G)+labs(title="True G")
  auc<-round(auc(pred = P, label = G),3)
  grid.arrange(p1,p2,ncol=2, top=paste0("Recall=",TPR," (Obs=",TPRO," , Hid=",TPRH,
                                        ")\n Precision=",PPV," (Obs=",PPVO," , Hid=",PPVH,")",
                                        "\n AUC=",auc))
}

#' Plot function for nestor convergence
#'
#' @param nestorFit a fit from the nestor() function
#' @return visualiaztion of nestor convergence
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
#' nestorFit=nestor(data$Y, MO,SO, initList=initList, maxIter=3 , verbatim=1)
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
