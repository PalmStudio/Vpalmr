#' Plot the architectural parameters estimations
#'
#' @description Plot the estimations of the palm architectural parameters from field measurements
#' for control
#'
#' @param progeny Character. The progeny name
#' @param pdfdir  The path to the folder where to write the plots
#'
#' @return A series of plots
#'
#' @importFrom dplyr ungroup group_by summarise "%>%"
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' archi_param= estimate_archi(progeny= "DY",map= 47)
#' plot_params(archi_param)
#'}
#' @export
#'
#'
plot_params= function(progeny,pdfdir){

  #----     create the pdf file ----
  pdf(paste(pdfdir,'/',progeny,'_',map,'MAP',Sys.Date(),'.pdf',sep=''))

  # Parameters:

  #rank for simulation (visualisation of the functions)
  rankSimu=seq(0,60,1)

  #relative position on rachis (visualisation of the functions)
  PositionRelativeRachis=round(seq(0,1,0.001),3)

  #number of emmitted leaves simulated (visualisation of the functions)
  LeafNumberSimu=seq(0,300,1)

  #rachis length simulated (visualisation of the functions)
  rachisLengthSimu=seq(50,700,1)

  #relative value
  relSimu=seq(0,1,0.01)

  #####     Phyllotaxis   #####
  par(mar=c(5,5,1,1))
  plot(data=DataAll, Phylo~Progeny,
       xlab=expression(Progeny),ylab=expression(italic(phi) (degree)),
       cex.lab=1.5)
  points(data=DataAll,Phylo~Progeny,col=Progeny,pch=as.numeric(Progeny)+14)

  # Stem diameter
  plot(data=DataAll[DataAll$Progeny==progeny,], StemDiameter ~ RachisLength,
       col=Progeny,pch=20,main=paste(progeny),
       xlab='Rachis Length (cm)',ylab='Stem Basis Diameter (cm)',cex.lab=1.3)
  lines(f.sigmo(X= rachisLengthSimu, max=finalStemDiam,
                 slope=StemDiamSlope,infl=StemDiamInfl)~rachisLengthSimu,
         col=Progeny)

  # Stem height
  plot(data=StemH, StemHeight17 ~TotalEmitted,col=Progeny,pch=20,
       ylab='Stem height (cm)',xlab='Number of emitted leaves from planting',
       xlim=c(0,200),ylim=c(0,150),cex.lab=1.3)
  points(f.expo(X = LeafNumberSimu,y0=H0,coef = coefStemHeight)~LeafNumberSimu,
         col=StemH $Progeny,type='l')

  # Rachis length
  LeafNumberCalib=seq(nbLeafEmitted-nbFronds_M,nbLeafEmitted,1)

  plot(data= DataAll[DataAll$Progeny==progeny,],
       RachisLength~ LeafNumber,
       ylim=c(0,max(DataAll$RachisLength,na.rm=T)+100),
       xlim=c(0,max(DataAll$LeafNumber,na.rm=T)+40),
       col=Progeny,pch=20,ylab='Rachis Length (cm)',
       xlab='Number of emmitted leaves from planting',main=paste(progeny),cex.lab=1.3)
  points(data=Rachis,rachisLength_intercept + LeafNumberCalib* rachisLength_slope~ LeafNumberCalib,
         type='l',col=Progeny)
  legend('topleft',c(paste('Palm age=',nbLeafEmitted,'leaves emitted since planting')),bty='n')

  # petiol/rachis ratio
  plot(data=Pet, RatioPetiole~ RachisLength,ylim=c(0,0.5),xlim=c(0,max(DataAll$RachisLength,na.rm=T)),col=Progeny,pch=20,ylab='Ratio Petiol/Rachis length',xlab='Rachis Length (cm)',main=paste(progeny))
  if (nrow(Pet)>0){
    abline(h=mean(unique(Pet$RatioPetiole),na.rm=T),col=Pet$Progeny)}
  if (nrow(Pet)==0){
    abline(h=mean(DataAll[DataAll$Progeny== progeny & DataAll$FrondRank>17,]$RatioPetiole,na.rm=T),col= DataAll [DataAll$Progeny==progeny,]$Progeny)
    legend('top','No data available for the nbLeafEmitted requested',bty='n')}

  # B point position
  plot(data= Bpos, PosB ~ RachisLength,ylim=c(0,1),xlim=c(0,max(DataAll$RachisLength,na.rm=T)),col=Progeny,pch=20,ylab='Relative position of B point',xlab='Rachis Length (cm)',main=paste(progeny))
  if (nrow(Bpos)>0){
    abline(h=mean(unique(Bpos$PosB),na.rm=T),col= Bpos $Progeny)}
  if (nrow(Bpos)==0){
    abline(h=mean(DataAll[DataAll$Progeny==progeny,]$PosB,na.rm=T),
           col= DataAll[DataAll$Progeny==progeny,]$Progeny)
    legend('top','No data available for the nbLeafEmitted requested',bty='n')}

  # Number of leaflets
  plot(data= Nb, Nb_leaflets ~ RachisLength,
       xlim=c(0,max(DataAll$RachisLength,na.rm=T)),
       ylim=c(0,max(DataAll$Nb_leaflets,na.rm=T)),
       col=Progeny,pch=20,xlab='Rachis Length (cm)',
       ylab='Number of leaflets',main=paste(progeny))
  points(f.sigmo(X= rachisLengthSimu, max = coef.nbLeaflets_mean['nbMax'],
                 slope = coef.nbLeaflets_mean['nbSlope'],
                 infl = coef.nbLeaflets_mean['nbInfl'])~rachisLengthSimu,
         col=Nb$Progeny,type='l')

  plot(data=Dec, Decli_C~Rank,main=paste(progeny),col=Progeny,xlab='Rank',ylab='Declination C (degree)',ylim=c(0,120),xlim=c(0,60),pch=20)
  points(f.linear(X=rankSimu, intercept = decliC_intercept ,
                  slope = decliC_slope)~ rankSimu,
         col=Dec$Progeny,type='l')


  #-----close pdf----
  dev.off()

}
