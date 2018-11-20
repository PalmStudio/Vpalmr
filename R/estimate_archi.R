#' Estimate architectural parameters
#'
#' @description Estimate palm architectural parameters from field measurements, and export
#'              them as VPalm inputs.
#'
#' @param progeny Character. The progeny name
#' @param map     Integer. The desired age of the Palm for simulation, in Month After Planting.
#' @param Dec     Declination angle data.frame for all progenies
#' @param output  The path to the folder where to write the outputs
#'
#' @details This function is highly adapted to a given input format. Please make sure to make
#'          the input formats as required.
#'
#' @return VPalm inputs
#'
#' @importFrom dplyr ungroup group_by summarise "%>%"
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' estimate_archi(progeny= "DY",map= 47)
#'}
#' @export
#'
estimate_archi= function(progeny,map,Dec){

  # selection of plant age in number of leaves emitted from MAP (Month After Planting)
  nbLeafEmitted= Parameter[Parameter$Progeny==progeny & Parameter$MAP==map,]$nbLeaves

  message(paste(nbLeafEmitted,'leaves are emitted at',map,'month after planting'))

  #number of fronds in 3D mock ups
  nbFronds_M=45

  #___________________CONTROL_PARAMETERS___________________#####

  #epsilon value for matrix inversion
  epsilon=10^-6

  #number of iterations for parameters estimation
  param_controle= nlme::lmeControl(maxIter=500000, niterEM=25000)


  ####____________PARAMETRES_ESTIMATION__________________________#####


  ####################STEM_SCALE#################


  #####     Stem diameter   #####
  StemD=DataAll[DataAll$Progeny==progeny,]

  StemDiam.nls=nls(data=StemD, StemDiameter~f.sigmo(X= RachisLength,max= finalStemDiam,slope= StemDiamSlope,infl= StemDiamInfl),start=c(finalStemDiam=300, StemDiamSlope=0.01, StemDiamInfl=200))
  residStemDiam=summary(StemDiam.nls)$sigma

  #####     Stem height   #####
  StemH=DataAll[DataAll$Progeny==progeny,]
  H0=5

  model.stemHeight=nls(data=StemH[!is.na(StemH$StemHeight17),], StemHeight17 ~ f.expo(X = TotalEmitted,y0=5,coef =  coefStem),start=c(coefStem =0.1))

  coefStemHeight=coef(model.stemHeight)
  residStemHeight=summary(model.stemHeight)$sigma

  ##################################LEAF_SCALE##########################################################

  #####     Rachis length   #####
  Rachis= DataAll[DataAll$TotalEmitted<= nbLeafEmitted+30 &
                    DataAll$TotalEmitted>= nbLeafEmitted-30 &
                    DataAll$Progeny==progeny & !is.na(DataAll$RachisLength),]

  rachisLength.lme=lme(RachisLength ~ LeafNumber,data= Rachis[!is.na(Rachis$RachisLength),],random=~1+ LeafNumber |TreeNumber,method='ML',control=param_controle)

  ###estimation of parameters
  rachisLength_intercept=summary(rachisLength.lme)$coefficients$fixed[1]
  rachisLength_slope=summary(rachisLength.lme)$coefficients$fixed[2]
  coef.rachisLength_mean=c(rachisLength_intercept, rachisLength_slope)
  SigmaR_rachisLength=summary(rachisLength.lme)$sigma

  label.rachisLength.lme =colnames(summary(rachisLength.lme)$coefficients$random$TreeNumber)

  SdG1_rachisLength <-as.numeric(VarCorr(rachisLength.lme)['(Intercept)','StdDev'])
  SdG2_rachisLength <-as.numeric(VarCorr(rachisLength.lme)['LeafNumber','StdDev'])
  corG_rachisLength <-as.numeric(VarCorr(rachisLength.lme)['LeafNumber','Corr'])

  ###variance-covariance matrix
  MatG_rachisLength=
    matrix(data=c(SdG1_rachisLength^2, SdG1_rachisLength* SdG2_rachisLength* corG_rachisLength,
                  SdG1_rachisLength* SdG2_rachisLength* corG_rachisLength, SdG2_rachisLength ^2),
           nrow=length(coef.rachisLength_mean),ncol=length(coef.rachisLength_mean),
           dimnames=list(label.rachisLength.lme, label.rachisLength.lme))


  ##### Ratio petiol/rachis  #####
  Pet= DataAll[DataAll$TotalEmitted<= nbLeafEmitted+60 &
                 DataAll$TotalEmitted>= nbLeafEmitted-60 &
                 DataAll$FrondRank>17 & DataAll$Progeny==progeny,]

  ###check for error inputs
  if (nrow(na.omit(Pet[Pet$RatioPetiole>0.5,]))>0){
    Pet$ID=paste('tree ',Pet$TreeNumber,' leaf ', Pet$LeafIndex,sep='')
    print('!!!!!PROBABLE ERROR IN PETIOLE LENGTH OR RACHIS LENGTH FOR:')
    print(paste(Pet[Pet$RatioPetiole>0.5,]$ID))}


  #####     B point position  #####
  Bpos= DataAll[DataAll$TotalEmitted<= nbLeafEmitted+60 & DataAll$TotalEmitted>= nbLeafEmitted-60 &  DataAll$FrondRank>17 & DataAll$Progeny==progeny & !is.na(DataAll $PosB),]

  #####     Number of leaflets   #####
  Nb=DataAll[DataAll$Progeny==progeny & !is.na(DataAll$Nb_leaflets),]

  nbLeaflets.nls=nls(data=Nb, Nb_leaflets ~ f.sigmo(X= RachisLength,max= nbMax,slope= nbSlope,infl= nbInfl),start=list(nbMax =100,nbSlope=0.01,nbInfl= 140),upper=c(nbMax=180,nbSlope=10,nbInfl=300),lower= list(nbMax=0,nbSlope=-10,nbInfl=0),algorithm='port')

  nbMax =coef(nbLeaflets.nls)['nbMax']
  nbSlope =coef(nbLeaflets.nls)['nbSlope']
  nbInfl =coef(nbLeaflets.nls)['nbInfl']
  coef.nbLeaflets_mean=c(nbMax, nbSlope, nbInfl)

  SigmaR_nbLeaflets=summary(nbLeaflets.nls)$sigma

  #####    Declination at C point   #####

  Dec=Dec[Dec$Progeny==progeny,]

  decliC.lme= lme(data= Dec[!is.na(Dec$Decli_C),],
                  Decli_C ~ Rank,random=~1+ Rank |TreeNumber,
                  method='ML',control=param_controle)

  decliC_intercept= summary(decliC.lme)$coefficients$fixed[1]
  decliC_slope= summary(decliC.lme)$coefficients$fixed[2]
  coef.decliC_mean= c(decliC_intercept, decliC_slope)
  SigmaR_decliC= summary(decliC.lme)$sigma

  label.decliC.lme= colnames(summary(decliC.lme)$coefficients$random$TreeNumber)

  SdG1_decli <-as.numeric(VarCorr(decliC.lme)['(Intercept)','StdDev'])
  SdG2_decli <-as.numeric(VarCorr(decliC.lme)['Rank','StdDev'])
  corG_decli<-as.numeric(VarCorr(decliC.lme)['Rank','Corr'])

  ###variance-covariance matrix
  MatG_decliC=matrix(data=c(SdG1_decli^2,SdG1_decli*SdG2_decli* corG_decli,SdG1_decli*SdG2_decli* corG_decli,SdG2_decli^2),nrow=length(coef.decliC_mean),ncol=length(coef.decliC_mean),dimnames=list(label.decliC.lme, label.decliC.lme))

  #####     Leaf curavture   #####
  Curve<-fread("1-Data/Archi/LeafCurvature_SMSE14.csv",header=TRUE,sep=";",dec=",", data.table= F)
  Curve$Obs_Date=as.Date(Curve$Obs_Date,'%d/%m/%y')

  Curve$Y_distance_soil_cm= Curve$Y_distance_cm+Curve$Height_O

  # Changing the reference point to C point, and removing O point
  Curve%<>%
    group_by(TreeNumber,Rank)%>%
    mutate(Y_distance_cm= Y_distance_cm-Y_distance_cm[Point=='C'],
           X_distance_cm= X_distance_cm-X_distance_cm[Point=='C'],
           Z_distance_cm= Z_distance_cm-Z_distance_cm[Point=='C'])%>%
    filter(Point!='O' & Complet!='Non' & Manipe!='Yes')

  #rachis length estimation
  Curve%<>%
    group_by(TreeNumber,Rank)%>%
    mutate(PositionRachisEstim=
             object_lenght(X_distance_cm,Y_distance_cm,
                           method= "smooth.spline")$Distance%>%cumsum,
           RachisLength= max(PositionRachisEstim),
           RelativePositionRachisEstim= PositionRachisEstim/RachisLength)
  # NB: The rachis lenght can also be found by finding the value of
  # PositionRachisEstim at point A

  # Leaflet angles at point A and C
  dataCurve=
    Curve[Curve$Progeny==pro,]%>%
    group_by(Progeny, TreeNumber, Rank)%>%
    summarise(angC= find_angle(X_distance_cm,Y_distance_cm)[1],
              angA= find_angle(X_distance_cm,Y_distance_cm)[2],
              rachisLength= mean(RachisLength))
  # Could use Curve to compute it for all progenies


  paramTree=NULL

  for (t in unique(CurvePro$TreeNumber)){

    data=CurvePro[CurvePro$TreeNumber==t,]

    optimDist<-NULL
    parameters<-list(coefCurv=0.5,decAInfl=30)
    try({
      optimDist<-optim(par=list(coefCurv=0.5,decAInfl=30),fn=SumDist,method="L-BFGS-B")
      #print(optimDist)
      parameters<-list(decAInfl= optimDist$par[["decAInfl"]],
                       coefCurv= optimDist$par[["coefCurv"]],
                       value= optimDist$value/nrow(data),
                       conv=optimDist$convergence
      )
    })
    result=data.frame(Progeny=progeny,TreeNumber=t,parameters)
    paramTree=rbind(paramTree,result)
  }

  plot(NA,ylim=c(-100,500),xlim=c(0,500),main=paste(progeny),pch=20,ylab='Y',xlab='X')

  for (t in unique(CurvePro$TreeNumber)){

    data=CurvePro[CurvePro$TreeNumber==t,]

    for (r in unique(data$Rank)){
      data_sub=data[data$Rank==r,]
      points(data=data_sub,Y_distance_cm ~ X_distance_cm,type='o',pch=20)
      points(f.leafCurvature(positionRelativeRachis= relSimu, angC= dataCurve[dataCurve$Progeny==progeny & dataCurve$TreeNumber==paste(t) & dataCurve$Rank==r,]$angC, angA= f.sigmo(dataCurve[dataCurve$Progeny==progeny & dataCurve$TreeNumber==paste(t) & dataCurve$Rank==r,]$angC,max=160,slope=0.02,infl=paramTree[paramTree$Tree==paste(t),]$decAInfl), coefCurv=paramTree[paramTree$Tree==paste(t),]$coefCurv, rachisLength=dataCurve[dataCurve$Progeny==progeny & dataCurve$TreeNumber==paste(t) & dataCurve$Rank==r,]$rachisLength),col=2,type='l')
      points(f.leafCurvature(positionRelativeRachis= relSimu, angC= dataCurve[dataCurve$Progeny==progeny & dataCurve$TreeNumber==paste(t) & dataCurve$Rank==r,]$angC, angA= f.sigmo(dataCurve[dataCurve$Progeny==progeny & dataCurve$TreeNumber==paste(t) & dataCurve$Rank==r,]$angC,max=160,slope=0.02,infl=paramTree[paramTree$Tree==paste(t),]$decAInfl), coefCurv=paramTree[paramTree$Tree==paste(t),]$coefCurv, rachisLength=dataCurve[dataCurve$Progeny==progeny & dataCurve$TreeNumber==paste(t) & dataCurve$Rank==r,]$rachisLength)[round(data_sub$RelativePositionRachisEstim,2)*100,],col=2,type='p',pch=20)

    }

  }
  legend('topright',c('obs','simu'),pch=20,col=c(1,2),lty=1,bty='n')


  #####     Declination at A point   #####
  decliA= dataCurve[dataCurve$Progeny==progeny,]


  decMaxA =180
  decSlopeA =0.01

  decliA.nls <- nls(data= decliA, angA ~ f.sigmo(X=angC,max= decMaxA, slope = decSlopeA, infl= decInflA),start=list(decInflA =mean(paramTree$decAInfl)),upper=c(decInflA =120),lower= list(infl=0),algorithm='port')

  decInflA =coef(decliA.nls)['decInflA']

  SigmaR_decliA=summary(decliA.nls)$sigma

  plot(data= decliA, angA~ angC,col=CurvePro$Progeny,pch=20,ylim=c(0,180),xlim=c(0,120),ylab='Declination at A point (degree)',xlab='Declination at C point (degree)',main=paste(progeny))

  points(f.sigmo(X=seq(0,120,1),max= decMaxA, slope = decSlopeA, infl= decInflA)~ seq(0,120,1),col= CurvePro $Progeny,type='l')

  #####     Rachis twist  #####
  Tor<- read.csv('1-Data/Archi/AnglesC&A_SMSE_Nov14.csv',header=T,sep=";",dec=',')
  Tor$Observation_Date=as.Date(Tor $Observation_Date,format='%d/%m/%Y')

  Tor=Tor[Tor$Progeny==progeny,]

  ###absolute value of the twist
  Tor$TwistA=abs(Tor$TwistA)


  #####     Rachis deviation  #####
  Dev=Curve[Curve$Point=='A' & !is.na(Curve$Z_distance_cm) & Curve$Progeny==progeny,]

  #degree of deviation
  Dev$DevA_deg=rep(0,nrow(Dev))
  for (t in unique(Dev$TreeNumber)){
    for (r in unique(Dev[Dev$TreeNumber==t,]$Rank)){
      Dev[Dev$TreeNumber==t & Dev$Rank==r,]$DevA_deg=f.devA(devCmObs= Dev[Dev$TreeNumber==t & Dev$Rank==r,]$Z_distance_cm, rachisLength= Dev[Dev$TreeNumber==t & Dev$Rank==r,]$RachisLength)
    }
  }


  #############################LEAFLET_SCALE###################################

  #####     Leaflets position on rachis    #####
  ####calibration from the MAP corresponding to nbLeafEmitted
  MAPmin=min(DataAll[DataAll$TotalEmitted<= nbLeafEmitted+20 & DataAll$TotalEmitted>= nbLeafEmitted-20 & DataAll$Progeny==progeny & DataAll$areaData=='Yes',]$MonthAfterPlanting)
  MAPmax=max(DataAll[DataAll$TotalEmitted<= nbLeafEmitted+20 & DataAll$TotalEmitted>= nbLeafEmitted-20 & DataAll$Progeny==progeny & DataAll$areaData=='Yes',]$MonthAfterPlanting)

  ##leaf area data available doesn't match because nbLeafEmitted too small--> take the first data
  if (MAPmin<min(unique(Area$MAP))){
    Area=Area[Area$MAP==min(unique(Area$MAP)),]
  }

  ##leaf area data available doesn't match because nbLeafEmitted too much--> take the last data
  if (MAPmax>max(unique(Area$MAP))){
    Area=Area[Area$MAP==max(unique(Area$MAP)),]
  }

  #leaf area data available matches with nbLeafEmitted
  if (nrow(Area[Area$MAP<= MAPmax & Area $MAP>= MAPmin,])>0){
    Area=Area[Area$MAP<= MAPmax & Area $MAP>= MAPmin,]
  }


  ###Number of leaflets per leaf
  Dispo=Area


  Dispo$TotalLeaflets=rep(NA,nrow(Dispo))

  for (t in unique(Dispo$TreeNumber)){
    for (l in unique(Dispo[Dispo$TreeNumber==t,]$LeafIndex)){
      Dispo[Dispo$TreeNumber==t & Dispo$LeafIndex ==l,]$TotalLeaflets=sum(Dispo[Dispo$TreeNumber==t & Dispo$LeafIndex ==l & Dispo$PositionOnLeaflet==0,]$NbLeaflets)
    }
  }

  ###Leaflet rank on rachis
  Dispo$LeafletRank=rep(NA,nrow(Dispo))
  for (t in unique(Dispo$TreeNumber)){
    for (l in unique(Dispo[Dispo$TreeNumber==t,]$LeafIndex)){

      Dispo[Dispo$TreeNumber==t & Dispo$LeafIndex ==l & Dispo$Section==1,]$LeafletRank=Dispo[Dispo$TreeNumber==t & Dispo$LeafIndex ==l & Dispo$Section==1,]$LeafletRankOnSection

      for (s in 2:10){
        Dispo[Dispo$TreeNumber==t & Dispo$LeafIndex ==l & Dispo$Section==s,]$LeafletRank=sum(Dispo[Dispo$TreeNumber==t & Dispo$LeafIndex ==l & Dispo$PositionOnLeaflet==0 & Dispo$Section<s,]$NbLeaflets)+Dispo[Dispo$TreeNumber==t & Dispo$LeafIndex ==l & Dispo$PositionOnLeaflet==0 & Dispo$Section==s,]$LeafletRankOnSection
      }
    }
  }

  ###Relative leafletRank
  Dispo$RelativeLeafletRank=Dispo$LeafletRank/Dispo$TotalLeaflets

  ###Relative position of leaflet on rachis
  Dispo$RelativePositionRachis=Dispo$PositionOnRachis/Dispo$LeafLength

  ###get one line per leaflet (1 line per section)
  Dispo=Dispo[Dispo$PositionOnLeaflet==0,]

  Dispo$LeafNumber=rep(NA,nrow(Dispo))

  ###affect LeafNumber from LeafIndex
  for (t in unique(Dispo$TreeNumber)){
    for (l in  unique(Dispo[Dispo$TreeNumber==t,]$LeafIndex)){
      Dispo[Dispo$TreeNumber==t & Dispo$LeafIndex==l,]$LeafNumber=min(DataAll[DataAll$TreeNumber==t & DataAll$LeafIndex==l,]$LeafNumber,na.rm=T)
    }
  }

  Dispo=Dispo[Dispo$Progeny==progeny,]

  ###
  plot(data=Dispo, RelativePositionRachis ~ RelativeLeafletRank,col=Progeny,xlim=c(0,1),ylim=c(0,1),pch=1,xlab='Leaflet relative rank',ylab='Leafelt relative position',main=paste(progeny))

  dispo.nls=nls(data=Dispo[Dispo$Progeny==progeny,], RelativePositionRachis ~f.leaflet_dispo(relRank= RelativeLeafletRank, coefDispo= coefDispo),start=list(coefDispo=2),control=param_controle)
  coefDispo=coef(dispo.nls)
  coefDispo_SD=summary(dispo.nls)$parameters[,'Std. Error']
  residDispo=summary(dispo.nls)$sigma
  points(f.leaflet_dispo(relRank=relSimu, coefDispo= coefDispo)~ relSimu,col=Dispo$Progeny,lwd=2,type='l')

  plot(data=Dispo, PositionOnRachis ~ RelativeLeafletRank,col=Progeny,ylim=c(0,max(Dispo$PositionOnRachis,na.rm=T)),xlim=c(0,1),pch=20,xlab='Leaflet relative rank',ylab='Leafelt position on rachis (cm)',main=paste(progeny))
  for (t in unique(Dispo$TreeNumber)){
    for (l in unique(Dispo[Dispo$TreeNumber==t,]$LeafIndex)){
      rachisLength=unique(Dispo[Dispo$TreeNumber==t & Dispo$LeafIndex==l,]$LeafLength)
      points(f.leaflet_dispo(relRank=relSimu, coefDispo= coefDispo)* rachisLength ~ relSimu,col=Dispo$Progeny,lwd=1,type='l')
    }
  }
  legend('topleft',c(paste('Palm age=',nbLeafEmitted,'leaves emitted since planting')),bty='n')


  #####     Leaflet length  Bpoint#####
  LfbAll=DataAll[!is.na(DataAll$LeafletBLength),]

  LfB= DataAll[DataAll$Progeny==progeny ,]

  plot(data=LfB, LeafletBLength~ RachisLength,pch=20,col=LfB$Progeny,ylab='Lealet length at B point (cm)',xlab='Rachis length (cm)',ylim=c(0,100),xlim=c(0,500),main=paste(progeny))
  leaflet_length_B.lme=lme(LeafletBLength ~ RachisLength,data= LfB[!is.na(LfB$LeafletBLength),],random=~1+ RachisLength |TreeNumber,method='ML',control=param_controle)

  #extraction OF parametres
  Length_B_intercept=summary(leaflet_length_B.lme)$coefficients$fixed[1]
  Length_B_slope=summary(leaflet_length_B.lme)$coefficients$fixed[2]
  Length_B_cov=vcov(leaflet_length_B.lme)
  coef.length_B_mean=c(Length_B_intercept,Length_B_slope)
  SigmaR_length_B=summary(leaflet_length_B.lme)$sigma

  label.leaflet_length_B.lme=colnames(summary(leaflet_length_B.lme)$coefficients$random$TreeNumber)

  SdG1_Length_B<-as.numeric(VarCorr(leaflet_length_B.lme)['(Intercept)','StdDev'])
  SdG2_Length_B<-as.numeric(VarCorr(leaflet_length_B.lme)['RachisLength','StdDev'])
  corG_Length_B<-as.numeric(VarCorr(leaflet_length_B.lme)['RachisLength','Corr'])

  ###variance-covariance matrix
  MatG_Length_B=matrix(data=c(SdG1_Length_B^2,SdG1_Length_B*SdG2_Length_B*corG_Length_B,SdG1_Length_B*SdG2_Length_B*corG_Length_B,SdG2_Length_B^2),nrow=length(coef.length_B_mean),ncol=length(coef.length_B_mean),dimnames=list(label.leaflet_length_B.lme, label.leaflet_length_B.lme))

  points(Length_B_intercept+Length_B_slope* rachisLengthSimu~ rachisLengthSimu,col= LfB$Progeny,type='l')

  #chol(MatG_Length_B)

  #####     Leaflet width Bpoint #####
  par(mar=c(5,5,1,1))
  plot(data=LfB, LeafletBWidth~ RachisLength,pch=20,col=LfB$Progeny,ylab='Lealet width at B point (cm)',xlab='Rachis length (cm)',ylim=c(0,10),xlim=c(0,500),main=paste(progeny))

  leaflet_width_B.lme=lme(LeafletBWidth ~ RachisLength,data=LfB[!is.na(LfB$LeafletBLength),],random=~1+ RachisLength |TreeNumber,method='ML',control=param_controle)

  #extraction des parametres
  width_B_intercept=summary(leaflet_width_B.lme)$coefficients$fixed[1]
  width_B_slope=summary(leaflet_width_B.lme)$coefficients$fixed[2]
  width_B_cov=vcov(leaflet_width_B.lme)
  coef.width_B_mean=c(width_B_intercept,width_B_slope)
  SigmaR_width_B=summary(leaflet_width_B.lme)$sigma

  label.leaflet_width_B.lme=colnames(summary(leaflet_width_B.lme)$coefficients$random$TreeNumber)


  SdG1_width_B <-as.numeric(VarCorr(leaflet_width_B.lme)['(Intercept)','StdDev'])
  SdG2_width_B <-as.numeric(VarCorr(leaflet_width_B.lme)['RachisLength','StdDev'])
  corG_width_B <-as.numeric(VarCorr(leaflet_width_B.lme)['RachisLength','Corr'])

  ###variance-covariance matrix
  MatG_width_B=matrix(data=c(SdG1_width_B^2,SdG1_width_B*SdG2_width_B*corG_width_B,SdG1_width_B*SdG2_width_B*corG_width_B,SdG2_width_B^2),nrow=length(coef.width_B_mean),ncol=length(coef.width_B_mean),dimnames=list(label.leaflet_width_B.lme, label.leaflet_width_B.lme))

  points(width_B_intercept+width_B_slope* rachisLengthSimu~ rachisLengthSimu,col= LfB$Progeny,type='l')


  #####   Leaflet relative width #####
  LfAll=Area

  ##subset of Area with just leaflet length (PositionOnLeaflet=leaflet tip)
  LfAll=LfAll[LfAll$Width==0 & LfAll$PositionOnLeaflet!=0,]
  LfAll$Position_rachis_rel= LfAll$PositionOnRachis/LfAll$LeafLength
  LfAll$Length= LfAll$PositionOnLeaflet

  #get the longest leaflet per leaf
  LfAll$Max_length=rep(NA,nrow(LfAll))
  for (i in LfAll $TreeNumber){
    for (l in unique(LfAll[LfAll $TreeNumber==i,]$LeafIndex)){
      LfAll[LfAll $TreeNumber==i & LfAll $LeafIndex==l,]$Max_length=max(LfAll[LfAll$TreeNumber==i & LfAll $LeafIndex==l,]$Length)
    }
  }

  #get the relative leaflet length per leaf
  LfAll $Relative_length= LfAll $Length/LfAll $Max_length

  ###check for data input errors
  LfAll $leafID=paste(LfAll $TreeNumber,'Leaf', LfAll $LeafIndex)
  if (nrow(LfAll[LfAll $Relative_length ==1 & LfAll $Position_rachis_rel<0.2 | LfAll $Relative_length ==1 & LfAll $Position_rachis_rel>0.8,])>0){

    print('!!!!!!!PROBABLE ERRORS IN MAX LENGTH FOR:')
    print(paste('treeName',LfAll[LfAll $Relative_length ==1 & LfAll $Position_rachis_rel<0.2 | LfAll $Relative_length ==1 & LfAll $Position_rachis_rel>0.8,]$leafID, 'section', LfAll[LfAll $Relative_length ==1 & LfAll $Position_rachis_rel<0.2 | LfAll $Relative_length ==1 & LfAll $Position_rachis_rel>0.8,]$Section),sep=' ')
  }

  ####
  Lf=LfAll[LfAll$Progeny==progeny,]

  Length_obs=Lf$Relative_length
  Position_rachis_obs=Lf$Position_rachis_rel


  leafleftLength.nlme=nlme(data=Lf,Relative_length~f.leaflet_length(X=Position_rachis_rel,Ymax=1,Y0=L0,Yfin=Lfin,X_Ymax=Pos_Lmax),start=(list(fixed=c(L0=mean(Length_obs[Position_rachis_obs==min(Position_rachis_obs,na.rm=T)]),Lfin=mean(Length_obs[Position_rachis_obs==max(Position_rachis_obs,na.rm=T)]),Pos_Lmax=mean(Position_rachis_obs[Length_obs==max(Length_obs,na.rm=T)])))),fixed=L0+Lfin+Pos_Lmax~1,random=L0+Lfin+Pos_Lmax~1|TreeNumber,control=param_controle)

  L0=summary(leafleftLength.nlme)$coefficients$fixed['L0']
  Lfin=summary(leafleftLength.nlme)$coefficients$fixed['Lfin']
  Pos_Lmax=summary(leafleftLength.nlme)$coefficients$fixed['Pos_Lmax']
  Leaflet_relative_length_cov=vcov(leafleftLength.nlme)
  coef.relative_leaflet_length_mean=c(L0,Lfin,Pos_Lmax)
  SigmaR_relative_leaflet_length=summary(leafleftLength.nlme)$sigma

  label.leafleftLength.nlme=colnames(summary(leafleftLength.nlme)$coefficients$random$TreeNumber)

  SdG1_leafletLength <-as.numeric(VarCorr(leafleftLength.nlme)['L0','StdDev'])
  SdG2_leafletLength <-as.numeric(VarCorr(leafleftLength.nlme)['Lfin','StdDev'])
  SdG3_leafletLength <-as.numeric(VarCorr(leafleftLength.nlme)['Pos_Lmax','StdDev'])
  corG12_leafletLength<-as.numeric(VarCorr(leafleftLength.nlme)['Lfin','Corr'])
  corG13_leafletLength<-as.numeric(VarCorr(leafleftLength.nlme)['Pos_Lmax','Corr'])
  corG23_leafletLength<-as.numeric(VarCorr(leafleftLength.nlme)['Pos_Lmax',length(summary(leafleftLength.nlme)$coefficients$fixed)+1])

  ###variance-covariance matrix
  MatG_leafletLength=matrix(data=c(SdG1_leafletLength^2,SdG1_leafletLength*SdG2_leafletLength*corG12_leafletLength,SdG1_leafletLength*SdG3_leafletLength*corG13_leafletLength,SdG1_leafletLength*SdG2_leafletLength*corG12_leafletLength,SdG2_leafletLength^2,SdG2_leafletLength*SdG3_leafletLength*corG23_leafletLength,SdG1_leafletLength*SdG3_leafletLength*corG13_leafletLength,SdG2_leafletLength*SdG3_leafletLength*corG23_leafletLength,SdG3_leafletLength^2),nrow=length(coef.relative_leaflet_length_mean),ncol=length(coef.relative_leaflet_length_mean),dimnames=list(label.leafleftLength.nlme,label.leafleftLength.nlme))

  plot(data= Lf, Relative_length~ Position_rachis_rel,xlim=c(0,1),ylim=c(0,1),pch=20,col=Lf$Progeny,main=paste(progeny),xlab='Relative position on rachis',ylab='Relative leaflet length')

  points(f.leaflet_length(X= relSimu,Ymax=1,Y0=L0,Yfin=Lfin,X_Ymax=Pos_Lmax)~ relSimu,col=Lf$Progeny,type='l')
  legend('bottom',c(paste('Palm age=',nbLeafEmitted,'leaves emitted since planting')),bty='n')


  # chol(MatG_leafletLength)
  # eigen(MatG_leafletLength)

  #####     Leaflet relative widths  #####
  WidthAll<- Area

  #
  WidthAll=WidthAll[WidthAll$PositionOnLeaflet!=0,]
  WidthAll$Leaflet_length=rep(NA,nrow(WidthAll))
  WidthAll$Leaflet_max_width=rep(NA,nrow(WidthAll))

  for (i in levels(WidthAll $TreeNumber)){
    for (l in unique(WidthAll[WidthAll$TreeNumber==i,]$LeafIndex)){
      for (s in unique(WidthAll[WidthAll$TreeNumber==i & WidthAll$LeafIndex==l,]$Section)){
        WidthAll[WidthAll$TreeNumber==i & WidthAll$LeafIndex==l & WidthAll$Section==s,]$Leaflet_length=max(WidthAll[WidthAll$TreeNumber==i & WidthAll$LeafIndex==l & WidthAll$Section==s,]$PositionOnLeaflet)
        WidthAll[WidthAll$TreeNumber==i & WidthAll$LeafIndex==l & WidthAll$Section==s,]$Leaflet_max_width=max(WidthAll[WidthAll $TreeNumber==i & WidthAll $LeafIndex==l & WidthAll $Section==s,]$Width)
      }
    }
  }

  WidthAll$Position_rachis_rel=WidthAll$PositionOnRachis/WidthAll$LeafLength

  WidthAll$Width_rel=WidthAll$Width/WidthAll$Leaflet_max_width


  WidthAll=WidthAll[WidthAll$Width==0,]

  WidthAll$Max_max_width=rep(NA,nrow(WidthAll))
  for (i in unique(WidthAll$TreeNumber)){
    for (l in unique(WidthAll[WidthAll$TreeNumber==i,]$LeafIndex)){
      WidthAll[WidthAll$TreeNumber==i & WidthAll$LeafIndex==l,]$Max_max_width=max(WidthAll[WidthAll$TreeNumber==i & WidthAll$LeafIndex==l,]$Leaflet_max_width)
    }
  }

  WidthAll$Relative_max_width=WidthAll$Leaflet_max_width/WidthAll$Max_max_width

  ###check for data input errors
  WidthAll$leafID=paste(WidthAll$TreeNumber,'Leaf',WidthAll$LeafIndex)
  if (nrow(WidthAll[WidthAll$Relative_max_width ==1 & WidthAll$Position_rachis_rel<0.3 | WidthAll$Relative_max_width ==1 & WidthAll$Position_rachis_rel>0.8,])>0){

    print('!!!!!!!PROBABLE ERRORS IN WIDTH MAX LENGTH FOR:')
    print(paste(WidthAll[WidthAll$Relative_max_width ==1 & WidthAll$Position_rachis_rel<0.3 | WidthAll$Relative_max_width ==1 & WidthAll$Position_rachis_rel>0.8,]$leafID,WidthAll[WidthAll$Relative_max_width ==1 & WidthAll$Position_rachis_rel<0.3 |WidthAll$Relative_max_width ==1 & WidthAll$Position_rachis_rel>0.8,]$Section))
  }


  Width=WidthAll[WidthAll$Progeny==progeny ,]

  ####ajustment
  Width_obs=Width$Relative_max_width
  Position_rachis_obs=Width$Position_rachis_rel

  leafleftWidth.nlme=nlme(data=Width,Width_obs~f.leaflet_max_width(X=Position_rachis_obs,Ymax=1,Y0=W0,Yfin=Wfin,X_Ymax=Pos_Wmax),start=(list(fixed=c(W0=0.2,Wfin=mean(Width_obs[Position_rachis_obs==max(Position_rachis_obs,na.rm=T)]),Pos_Wmax=mean(Position_rachis_obs[Width_obs==max(Width_obs,na.rm=T)])))),fixed=W0+Wfin+Pos_Wmax~1,random=W0+Wfin+Pos_Wmax~1|TreeNumber,control=param_controle)

  W0=summary(leafleftWidth.nlme)$coefficients$fixed['W0']
  Wfin=summary(leafleftWidth.nlme)$coefficients$fixed['Wfin']
  Pos_Wmax=summary(leafleftWidth.nlme)$coefficients$fixed['Pos_Wmax']
  Leaflet_relative_Width_cov=vcov(leafleftWidth.nlme)
  coef.relative_leaflet_Width_mean=c(W0,Wfin,Pos_Wmax)
  SigmaR_relative_leaflet_Width=summary(leafleftWidth.nlme)$sigma

  label.leafleftWidth.nlme=colnames(summary(leafleftWidth.nlme)$coefficients$random$TreeNumber)

  SdG1_leafletWidth <-as.numeric(VarCorr(leafleftWidth.nlme)['W0','StdDev'])
  SdG2_leafletWidth <-as.numeric(VarCorr(leafleftWidth.nlme)['Wfin','StdDev'])
  SdG3_leafletWidth <-as.numeric(VarCorr(leafleftWidth.nlme)['Pos_Wmax','StdDev'])
  corG12_leafletWidth<-as.numeric(VarCorr(leafleftWidth.nlme)['Wfin','Corr'])
  corG13_leafletWidth<-as.numeric(VarCorr(leafleftWidth.nlme)['Pos_Wmax','Corr'])
  corG23_leafletWidth<-as.numeric(VarCorr(leafleftWidth.nlme)['Pos_Wmax',length(summary(leafleftWidth.nlme)$coefficients$fixed)+1])

  ###variance-covariance matrix
  MatG_leafletWidth=matrix(data=c(SdG1_leafletWidth^2,SdG1_leafletWidth*SdG2_leafletWidth*corG12_leafletWidth,SdG1_leafletWidth*SdG3_leafletWidth*corG13_leafletWidth,SdG1_leafletWidth*SdG2_leafletWidth*corG12_leafletWidth,SdG2_leafletWidth^2,SdG2_leafletWidth*SdG3_leafletWidth*corG23_leafletWidth,SdG1_leafletWidth*SdG3_leafletWidth*corG13_leafletWidth,SdG2_leafletWidth*SdG3_leafletWidth*corG23_leafletWidth,SdG3_leafletWidth^2),nrow=length(coef.relative_leaflet_Width_mean),ncol=length(coef.relative_leaflet_Width_mean),dimnames=list(label.leafleftWidth.nlme,label.leafleftWidth.nlme))

  plot(data= Width, Relative_max_width ~ Position_rachis_rel,xlim=c(0,1),ylim=c(0,1),pch=20,col=Progeny,main=paste(progeny),xlab='Relative position on rachis',ylab='Relative leaflet maximum width')

  points(f.leaflet_max_width(X= relSimu,Ymax=1,Y0=W0,Yfin=Wfin,X_Ymax=Pos_Wmax)~ relSimu,col=Width$Progeny,type='l')
  legend('bottom',c(paste('Palm age=',nbLeafEmitted,'leaves emitted since planting')),bty='n')

  #####     Leaflet axial angle  #####
  LftAngle<- read.csv("LeafDispositionComp_SMSE14.csv",header=T,sep=";",dec=',')
  LftAngle$Axial=AlKashi(LftAngle$Axial_cm,32.5)

  LftAngle$Rachis_length=rep(NA,nrow(LftAngle))
  for (i in 1:length(unique(LftAngle$TreeNumber))){
    ind=unique(LftAngle$TreeNumber)[i]
    LftAngle[LftAngle$TreeNumber==ind,]$Rachis_length=max(LftAngle[LftAngle$TreeNumber==ind,]$Position,na.rm=T)
  }

  ###Relative position of leaflet on rachis
  LftAngle$Position_rel=LftAngle$Position/LftAngle$Rachis_length

  LftAngle$Nb_leaflets=rep(NA,nrow(LftAngle))
  for (i in 1:length(unique(LftAngle$TreeNumber))){
    ind=unique(LftAngle$TreeNumber)[i]
    LftAngle[LftAngle$TreeNumber==ind,]$Nb_leaflets=max(LftAngle[LftAngle$TreeNumber==ind,]$Leaflet_rank,na.rm=T)
  }


  LftAngle=LftAngle[LftAngle$Progeny==progeny,]

  axialAngle.nlme=nlme(data= LftAngle[!is.na(LftAngle $Axial),],Axial~f.axialAngle(X=Position_rel,angleC=angleC,slopeC=slopeC,angleA=angleA),start=list(fixed=c(angleC=101,slopeC=-2,angleA=10)),fixed=angleC+slopeC+angleA~1,random=angleC+slopeC+angleA~1|TreeNumber,control=param_controle)

  angleC=summary(axialAngle.nlme)$coefficients$fixed['angleC']
  slopeC=summary(axialAngle.nlme)$coefficients$fixed['slopeC']
  angleA=summary(axialAngle.nlme)$coefficients$fixed['angleA']
  coef.axialAngle_mean=c(angleC,slopeC, angleA)
  axialAngle_cov=vcov(axialAngle.nlme)
  SigmaR_axialAngle=summary(axialAngle.nlme)$sigma

  label.axialAngle.nlme=colnames(summary(axialAngle.nlme)$coefficients$random$TreeNumber)

  SdG1_axialAngle <-as.numeric(VarCorr(axialAngle.nlme)['angleC','StdDev'])
  SdG2_axialAngle<-as.numeric(VarCorr(axialAngle.nlme)['slopeC','StdDev'])
  SdG3_axialAngle <-as.numeric(VarCorr(axialAngle.nlme)['angleA','StdDev'])
  corG12_axialAngle<-as.numeric(VarCorr(axialAngle.nlme)['slopeC','Corr'])
  corG13_axialAngle<-as.numeric(VarCorr(axialAngle.nlme)['angleA','Corr'])
  corG23_axialAngle<-as.numeric(VarCorr(axialAngle.nlme)['angleA',length(summary(axialAngle.nlme)$coefficients$fixed)+1])

  ###variance-covariance matrix
  MatG_axialAngle=matrix(data=c(SdG1_axialAngle^2,SdG1_axialAngle*SdG2_axialAngle*corG12_axialAngle,SdG1_axialAngle*SdG3_axialAngle*corG13_axialAngle,SdG1_axialAngle*SdG2_axialAngle*corG12_axialAngle,SdG2_axialAngle^2,SdG2_axialAngle*SdG3_axialAngle*corG23_axialAngle,SdG1_axialAngle*SdG3_axialAngle*corG13_axialAngle,SdG2_axialAngle*SdG3_axialAngle*corG23_axialAngle,SdG3_axialAngle^2),nrow=length(coef.axialAngle_mean),ncol=length(coef.axialAngle_mean),dimnames=list(label.axialAngle.nlme, label.axialAngle.nlme))

  plot(data= LftAngle, Axial~ Position_rel,col=Progeny,pch=20,xlim=c(0,1),ylim=c(0,120),ylab='Leaflet axial angle (deg)',xlab='Relative position on rachis',main=paste(progeny))
  points(f.axialAngle(X=relSimu,angleC=angleC,slopeC=slopeC,angleA=angleA)~relSimu,col=LftAngle$Progeny,type='l')

  #####     Leaflet radial angle  #####
  LftAngle$Alpha=90-AlKashi(LftAngle$Radial_cm,32.5)

  #estimation of the projected radia angle
  LftAngle$Radial_deg=atan(sin(LftAngle$Axial*pi/180)*tan(LftAngle$Alpha*pi/180))*180/pi

  ###mean radial angle vs relative position on rachis
  positionRelative=round(seq(0.1,1,0.1),1)

  dataAngle=NULL
  for (t in unique(LftAngle$Type))
    for (i in positionRelative){
      test=LftAngle[LftAngle$Type==t & LftAngle$Section/10==i,]
      dataAngle_sub=data.frame(type=t,position=i,mean=mean(test$Radial_deg,na.rm=T),sd=sd(test$Radial_deg,na.rm=T))
      dataAngle=rbind(dataAngle,dataAngle_sub)
    }

  dataAngle$up=dataAngle$mean+2*dataAngle$sd
  dataAngle$dwn=dataAngle$mean-2*dataAngle$sd

  ####

  radialHighSup.nls=nls(data=na.omit(dataAngle[dataAngle$type==1,]),up~f.radialAngle(x=position,A0=A0,Amax=Amax,Xm=0.5),start=list(A0=10,Amax=max(na.omit(dataAngle[dataAngle$type==1,])$up)),algorithm='port',control=list(maxiter=5000000,minFactor=0.000000000001,warnOnly=T),trace=F)

  radialHighInf.nls=nls(data=na.omit(dataAngle[dataAngle$type==0,]),up~f.radialAngle(x=position,A0=A0,Amax=Amax,Xm=0.5),start=list(A0=10,Amax=max(na.omit(dataAngle[dataAngle$type==0,])$up)),algorithm='port',control=list(maxiter=5000000,minFactor=0.000000000001,warnOnly=T),trace=F)

  radialLowSup.nls=nls(data=na.omit(dataAngle[dataAngle$type==0,]),dwn~f.radialAngle(x=position,A0=A0,Amax=Amax,Xm=0.5),start=list(A0=-10,Amax=max(na.omit(dataAngle[dataAngle$type==0,])$up)),algorithm='port',control=list(maxiter=5000000,minFactor=0.000000000001,warnOnly=T),trace=F)

  radialLowInf.nls=nls(data=na.omit(dataAngle[dataAngle$type==-1,]),dwn~f.radialAngle(x=position,A0=A0,Amax=Amax,Xm=0.5),start=list(A0=-10,Amax=max(na.omit(dataAngle[dataAngle$type==-1,])$up)),algorithm='port',control=list(maxiter=5000000,minFactor=0.000000000001,warnOnly=T),trace=F)

  plot(data= LftAngle, Radial_deg ~ Position_rel,col=Type+2,pch=20,xlim=c(0,1),ylim=c(-90,90),ylab='Leaflet axial angle (deg)',xlab='Relative position on rachis',main=paste(progeny))

  polygon(c(relSimu,rev(relSimu)),c(f.radialAngle(x= relSimu,A0=coef(radialHighSup.nls)['A0'],Amax=coef(radialHighSup.nls)['Amax'],Xm=0.5),rev(f.radialAngle(x= relSimu,A0=coef(radialHighInf.nls)['A0'],Amax=coef(radialHighInf.nls)['Amax'],Xm=0.5))),col='lightgreen',border=NA)

  polygon(c(relSimu,rev(relSimu)),c(f.radialAngle(x= relSimu,A0=coef(radialHighInf.nls)['A0'],Amax=coef(radialHighInf.nls)['Amax'],Xm=0.5),rev(f.radialAngle(x= relSimu,A0=coef(radialLowSup.nls)['A0'],Amax=coef(radialLowSup.nls)['Amax'],Xm=0.5))),col='coral',border=NA)

  polygon(c(relSimu,rev(relSimu)),c(f.radialAngle(x= relSimu,A0=coef(radialLowSup.nls)['A0'],Amax=coef(radialLowSup.nls)['Amax'],Xm=0.5),rev(f.radialAngle(x= relSimu,A0=coef(radialLowInf.nls)['A0'],Amax=coef(radialLowInf.nls)['Amax'],Xm=0.5))),col='grey',border=NA)

  points(data= LftAngle, Radial_deg ~ Position_rel,col=Type+2,pch=20)

  legend('bottomleft',title='Observation',c('up','med','dwn'),pch=20,col=c(3,2,1),bty='n')
  legend('bottom',title='Simulation range',c('up','med','dwn'),pch=NA,lwd=10,col=c('lightgreen','coral','grey'),bty='n')

  #####     Frequency of leaflets type  #####

  ###RDistibution of leaflets depending on their radial angle
  Rep=as.data.frame(table(LftAngle$Type,LftAngle$Section))

  colnames(Rep)=c('Type','Section','Freq')
  RepSection=as.data.frame(table(LftAngle$Section))
  colnames(RepSection)=c('Section','FreqTotal')
  Rep=cbind(Rep,FreqSection=rep(RepSection$FreqTotal,each=3))
  Rep$Position_rel=(as.numeric(Rep$Section))/10
  Rep$Prop=Rep$Freq/Rep$FreqSection

  plot(data=Rep[Rep$Type==1,] ,Prop~Position_rel,type='o',xlim=c(0,1),col=3,pch=20,ylim=c(0,1),main=paste(progeny),ylab='Leaflets relative frequency',xlab='Relative position on rachis')
  points(data=Rep[Rep$Type==0,] ,Prop~Position_rel,type='o',col=2,pch=20)
  points(data=Rep[Rep$Type==-1,] ,Prop~Position_rel,type='o',col=1,pch=20)
  legend('topleft',title='Observation',c('up','med','dwn'),pch=20,col=c(3,2,1),bty='n')


  #####     Leaflet shape  #####

  Shape<-Area
  Shape$Leaflet_length=rep(NA,nrow(Shape))
  Shape$Leaflet_max_width=rep(NA,nrow(Shape))

  ##leaflet length and max width
  for (i in levels(Shape$TreeNumber)){
    for (s in unique(Shape[Shape$TreeNumber==i,]$Section)){
      Shape[Shape$TreeNumber==i & Shape$Section==s,]$Leaflet_length=max(Shape[Shape$TreeNumber==i & Shape$Section==s,]$PositionOnLeaflet)
      Shape[Shape$TreeNumber==i & Shape$Section==s,]$Leaflet_max_width=max(Shape[Shape$TreeNumber==i & Shape$Section==s,]$Width)
    }
  }

  #normalized data
  Shape$Position_rachis_rel=Shape$PositionOnRachis/Shape$LeafLength
  Shape$Position_leaflet_rel=Shape$PositionOnLeaflet/Shape$Leaflet_length
  Shape$Width_rel=Shape$Width/Shape$Leaflet_max_width


  #adjustment for each section of the leaf (1 leaflet per section) --> retur the value of xm et ym for each leaflet
  dataXmYm=NULL

  for (t in unique(Shape$TreeNumber)){
    for (s in unique(Shape[Shape$TreeNumber==t,]$Section)){

      Leaflet_width_obs=Shape[Shape$TreeNumber==t  & Shape$Section==s,]$Width_rel
      Leaflet_position_obs=Shape[Shape$TreeNumber==t & Shape$Section==s,]$Position_leaflet_rel

      ajust=try(nls(Leaflet_width_obs~f.beta(x=Leaflet_position_obs,xm=xm_estim,ym=ym_estim),start=list(xm_estim=0.5,ym_estim=0.5),control=list(maxiter=5000000,minFactor=0.000000000001,warnOnly=T),trace=F),silent=T)
      if (inherits(ajust,'try-error')){
        #print(paste('no fitting for treeNumber',t,'section',s))
        next
      }
      dataXmYm_sub=data.frame(Progeny=unique(Shape[Shape$TreeNumber==t,]$Progeny),TreeNumber=t, PositionRelative=s/10,xm=coef(ajust)['xm_estim'],ym=coef(ajust)['ym_estim'])
      dataXmYm=rbind(dataXmYm,dataXmYm_sub)
    }
  }


  par(mfcol=c(1,1))
  plot(data= dataXmYm[dataXmYm$Progeny==progeny,],xm~ PositionRelative,ylim=c(0,1),pch=20,col=1,main=paste(progeny),ylab='Parameter value',xlab='Relative position on rachis')
  points(data= dataXmYm[dataXmYm$Progeny==progeny,],ym~ PositionRelative,pch=20,col=2)
  xm.lm=lm(data=dataXmYm[dataXmYm$Progeny==progeny,],xm~PositionRelative)
  ym.lm=lm(data=dataXmYm[dataXmYm$Progeny==progeny,],ym~PositionRelative)
  abline(xm.lm,col=1)
  abline(ym.lm,col=2)

  legend('topright',title='Parameters',c('xm','ym'),col=c(1,2),bty='n',pch=20)

  ######################################################NERVE SHAPE#########################################################
  #####     Nerve width at C point   #####
  PetioleSectionC=read.csv("Petiole_SMSE14.csv",header=T,sep=";",dec=',')
  PetioleSectionC$Obs_Date=as.Date(PetioleSectionC$Obs_Date,'%d/%m/%y')
  PetioleSectionC$Petiole_width_C_cm=PetioleSectionC$Petiole_width_C_mm/10
  PetioleSectionC$Petiole_height_C_cm=PetioleSectionC$Petiole_height_C_mm/10
  PetioleSectionC=PetioleSectionC[!is.na(PetioleSectionC$Petiole_width_C_cm) & !is.na(PetioleSectionC$CtoA),]
  PetioleSectionC=PetioleSectionC[PetioleSectionC$Progeny==progeny,]

  ####Ajustment + graph
  plot(data=PetioleSectionC,Petiole_width_C_cm~CtoA,ylab='Width section at C point (cm)',xlab='Rachis length (cm)',pch=20,main=paste(progeny),cex.lab=1.3,col=Progeny)

  petioleWidthC.lme=lme(data=PetioleSectionC,Petiole_width_C_cm~CtoA,random=~1+CtoA |TreeNumber,method='ML',control=param_controle)
  widthC_intercept=summary(petioleWidthC.lme)$coefficients$fixed[1]
  widthC_slope=summary(petioleWidthC.lme)$coefficients$fixed[2]
  coef.widthC_mean=c(widthC_intercept,widthC_slope)
  widthC_cov=vcov(petioleWidthC.lme)
  SigmaR_widthC=summary(petioleWidthC.lme)$sigma

  label.petioleWidthC.lme=colnames(summary(petioleWidthC.lme)$coefficients$random$TreeNumber)

  SdG1_widthC <-as.numeric(VarCorr(petioleWidthC.lme)['(Intercept)','StdDev'])
  SdG2_widthC<-as.numeric(VarCorr(petioleWidthC.lme)['CtoA','StdDev'])
  corG_widthC<-as.numeric(VarCorr(petioleWidthC.lme)['CtoA','Corr'])

  ###variance-covariance matrix
  MatG_widthC=matrix(data=c(SdG1_widthC^2,SdG1_widthC*SdG2_widthC*corG_widthC,SdG1_widthC*SdG2_widthC*corG_widthC,SdG2_widthC^2),nrow=length(coef.widthC_mean),ncol=length(coef.widthC_mean),dimnames=(list(label.petioleWidthC.lme, label.petioleWidthC.lme)))

  points(widthC_intercept+ widthC_slope*seq(0,700,1)~seq(0,700,1),type='l',col= PetioleSectionC $Progeny)

  #####     Nerve height   #####
  RachisHeight=read.csv("Torsion_SMSE14.csv",header=T,sep=";",dec=',')
  RachisHeight$Obs_Date=as.Date(RachisHeight$Obs_Date,'%d/%m/%y')
  RachisHeight$Position_rachis_rel=RachisHeight$Position_rachis/RachisHeight$Rachis_length
  RachisHeight$Petiole_width_cm=RachisHeight$Petiole_width/10
  RachisHeight$Petiole_height_cm=RachisHeight$Petiole_height/10
  RachisHeight=RachisHeight[!is.na(RachisHeight$Petiole_width_cm),]
  RachisHeight$Petiole_max_width=rep(NA,nrow(RachisHeight))

  for (p in unique(RachisHeight$TreeNumber)){
    RachisHeight[RachisHeight$TreeNumber==p,]$Petiole_max_width=max(RachisHeight[RachisHeight$TreeNumber==p,]$Petiole_width,na.rm=T)
  }
  RachisHeight$Petiole_max_height=rep(NA,nrow(RachisHeight))
  for (p in unique(RachisHeight$TreeNumber)){
    RachisHeight[RachisHeight$TreeNumber==p,]$Petiole_max_height=max(RachisHeight[RachisHeight$TreeNumber==p,]$Petiole_height,na.rm=T)
  }

  RachisHeight$Petiole_relative_width=RachisHeight$Petiole_width/RachisHeight$Petiole_max_width
  RachisHeight$Petiole_relative_height=RachisHeight$Petiole_height/RachisHeight$Petiole_max_height

  rachisRelativeHeight.nlme <-nlme(data=RachisHeight[!is.na(RachisHeight$Petiole_relative_height),], Petiole_relative_height ~f.rachisRelativeHeight(PositionRelativeRachis=Position_rachis_rel,a),fixed=a~1,random=a~1|TreeNumber,start=list(fixed=c(a=0.1)),control=param_controle)

  coef.rachisRelativeHeight_mean=summary(rachisRelativeHeight.nlme)$coefficients$fixed
  rachisRelativeHeight_cov=vcov(rachisRelativeHeight.nlme)
  SigmaR_rachisRelativeHeight=summary(rachisRelativeHeight.nlme)$sigma

  SdG1_rachisRelativeHeight<-intervals(rachisRelativeHeight.nlme)$reStruct$TreeNumber[1,'est.']

  plot(data=RachisHeight[!is.na(RachisHeight$Petiole_relative_height),], Petiole_relative_height ~ Position_rachis_rel,pch=20,ylab='Petiol relative height',xlab='Relative position on rachis',cex.lab=1.3)
  points(f.rachisRelativeHeight(seq(0,1,0.1), coef.rachisRelativeHeight_mean)~seq(0,1,0.1),type='l')

  #-----close pdf----
  dev.off()

  ####____________PARAMETERS OUTPUT  _________________#####

  ###parameters file with inter-tree variability

  ArchiOutput=function(Progeny,nbLeafEmitted,seed){
    #-----------------------------------FIXED_PARAMETERS-----------------------------------------------

    trunkBending_M=0
    trunkBending_SD=0

    # internodeFinalLength=2

    laminaAngle = 140

    ###leafletStifness=0--> soft (pending) leaflets
    #leafletStifness=10000--> erect leaflets

    leafletStiffness =    1500
    leafletStiffness_SD = 7000

    ##nerve shape section at leaf extremities
    frondBaseWidth = 30
    frondtipWidth = 0.3
    frondBaseHeight = 10
    frondTipHeight = 0.5

    nbInflorescences = 0 ###not yet implemented

    ###twist
    rachisTwistCoef=3
    rachisDevCoef=2

    #-------------------------------FIXED_PARAMETERS_PER_PROGENY---------------------------------------
    ##################
    #stem phylotaxis #
    ##################
    frondPhyllo_M= mean(Phylo$Phylo,na.rm=T)
    frondPhyllo_SD= sd(Phylo$Phylo,na.rm=T)

    ###############
    # stem height #
    ###############
    residStemHeight

    ################
    #stem diameter #
    ################

    stemDiamMax=coef(StemDiam.nls)['finalStemDiam']
    stemDiamSlope= coef(StemDiam.nls)['StemDiamSlope']
    stemDiamInfl= coef(StemDiam.nls)['StemDiamInfl']


    ##################
    #number of frond #
    ##################
    nbFronds_M
    nbFronds_SD=0


    #############
    #frond twist#
    #############
    rachisTwistFinalAngle_M=mean(Tor$TwistA,na.rm=T)
    rachisTwistFinalAngle_SD=sd(Tor$TwistA,na.rm=T)

    #################
    #frond deviation#
    #################
    rachisDevFinalAngle_M=mean(Dev$DevA_deg)
    rachisDevFinalAngle_SD=sd(Dev$DevA_deg)

    #######################################
    #parameters ratio petiole/rachis final#
    #######################################
    petioleRachisRatio_M=ifelse(is.na(mean(unique(Pet$RatioPetiole),na.rm=T)),mean(DataAll[DataAll$Progeny== Progeny & DataAll$FrondRank>17 ,]$RatioPetiole,na.rm=T),mean(unique(Pet$RatioPetiole),na.rm=T))
    petioleRachisRatio_SD=ifelse(is.na(sd(unique(Pet$RatioPetiole),na.rm=T)),sd(DataAll[DataAll$Progeny== Progeny & DataAll$FrondRank>17 ,]$RatioPetiole,na.rm=T),sd(unique(Pet$RatioPetiole),na.rm=T))

    #####################################
    #B point relative position on rachis#
    #####################################
    pointBrelativePosition_M=ifelse(is.na(mean(Bpos$PosB,na.rm=T)),mean(DataAll[DataAll$Progeny== Progeny,]$PosB,na.rm=T),mean(Bpos$PosB,na.rm=T))
    pointBrelativePosition_SD=ifelse(is.na(sd(Bpos$PosB,na.rm=T)),sd(DataAll[DataAll$Progeny== Progeny,]$PosB,na.rm=T),sd(Bpos$PosB,na.rm=T))

    ################
    #leaflet number#
    ################
    nbMax= coef.nbLeaflets_mean['nbMax']
    nbSlope = coef.nbLeaflets_mean['nbSlope']
    nbInfl= coef.nbLeaflets_mean['nbInfl']
    nbLeaflets_SDP=round(SigmaR_nbLeaflets)

    #################################
    #parameters Leaflet radial angle#
    #################################
    leafletRadialHighA0Sup=coef(radialHighSup.nls)['A0']
    leafletRadialHighAmaxSup=coef(radialHighSup.nls)['Amax']

    leafletRadialHighA0Inf=coef(radialHighInf.nls)['A0']
    leafletRadialHighAmaxInf=coef(radialHighInf.nls)['Amax']

    leafletRadialLowA0Sup=coef(radialLowSup.nls)['A0']
    leafletRadialLowAmaxSup=coef(radialLowSup.nls)['Amax']

    leafletRadialLowA0Inf=coef(radialLowInf.nls)['A0']
    leafletRadialLowAmaxInf=coef(radialLowInf.nls)['Amax']

    ############################
    #Frequency of leaflets type#
    ############################
    leafletFrequencyHigh=Rep[Rep$Type==1,]$Prop
    leafletFrequencyLow=Rep[Rep$Type==-1,]$Prop

    ###############
    #leaflet shape#
    ###############
    xm_intercept=coef(xm.lm)['(Intercept)']
    xm_slope=coef(xm.lm)['PositionRelative']
    ym_intercept=coef(ym.lm)['(Intercept)']
    ym_slope=coef(ym.lm)['PositionRelative']

    #################################
    #ratio height/width sectionC	#
    #################################
    heightWidthCRatio=mean(PetioleSectionC$Petiole_height_C_cm/PetioleSectionC$Petiole_width_C_cm)


    #------------------------------------PARAMETERS_PER_TREE---------------------------------------------

    ######
    #Seed#
    ######
    set.seed(seed)


    ###############
    #rachis length#
    ###############
    diag(MatG_rachisLength)=diag(MatG_rachisLength)+epsilon

    coef.rachisLength_sd<-matrix(data=rnorm(n=ncol(MatG_rachisLength)),ncol=ncol(MatG_rachisLength),dimnames=list(Progeny= Progeny,Effet=rownames(MatG_rachisLength)))%*% chol(MatG_rachisLength)

    coef.rachisLength_simu= coef.rachisLength_mean + coef.rachisLength_sd

    rachisLengthIntercept=coef.rachisLength_simu[,'(Intercept)']
    rachisLengthSlope=coef.rachisLength_simu[,'LeafNumber']
    rachisLength_SDP=SigmaR_rachisLength



    #####################
    #C angle declination#
    #####################

    coef.decliC_sd<-
      matrix(data=rnorm(n=ncol(MatG_decliC)),ncol=ncol(MatG_decliC),
             dimnames=list(Progeny= Progeny,Effet=rownames(MatG_decliC)))%*%
      chol(MatG_decliC)

    coef.decliC_simu= coef.decliC_mean + coef.decliC_sd

    decliCintercept =coef.decliC_simu[Progeny,'(Intercept)']
    decliCslope =coef.decliC_simu[Progeny,'Rank']

    cPointAngle_SDP = SigmaR_decliC

    ########################
    #declination frond tip #
    ########################
    decMaxA
    decSlopeA
    decInflA

    APointAngle_SDP = SigmaR_decliA

    #####################
    #rachis curvature	#
    #####################

    coefCurve=rnorm(n=1,mean=mean(paramTree$coefCurv),sd=sd(paramTree$coefCurv))

    ###################
    #leaflets position#
    ###################

    coefDispo
    Dispo_SDP= residDispo

    ########################
    #leaflet length b point#
    ########################
    diag(MatG_Length_B)=diag(MatG_Length_B)+2*epsilon

    coef.length_B_sd<-matrix(data=rnorm(n=ncol(MatG_Length_B)),ncol=ncol(MatG_Length_B),dimnames=list(Progeny=progeny,Effet=rownames(MatG_Length_B)))%*% chol(MatG_Length_B)

    coef.length_B_simu=coef.length_B_mean + coef.length_B_sd

    lenfletLengthAtBIntercept=coef.length_B_simu[,'(Intercept)']
    leafletLengthAtBSlope=coef.length_B_simu[,'RachisLength']

    ###################################
    #parameters  leaflet width b point#
    ###################################
    diag(MatG_width_B)=diag(MatG_width_B)+2*epsilon

    coef.width_B_sd<-matrix(data=rnorm(n=ncol(MatG_width_B)),ncol=ncol(MatG_width_B),dimnames=list(Progeny=Progeny,Effet=rownames(MatG_width_B)))%*% chol(MatG_width_B)

    coef.width_B_simu= coef.width_B_mean + coef.width_B_sd

    bWidthIntercept= coef.width_B_simu[,'(Intercept)']
    bWidthSlope =coef.width_B_simu[,'RachisLength']

    ##########################
    # relative leaflet length#
    ##########################
    ###when probleme of non semi positive matrix, add epsilon on the diagonal to make it closer to indentity matrice
    diag(MatG_leafletLength)=diag(MatG_leafletLength)+5*epsilon
    coef.leafletLength_sd<-matrix(data=rnorm(n=ncol(MatG_leafletLength)),ncol=ncol(MatG_leafletLength),dimnames=list(Progeny= Progeny,Effet=rownames(MatG_leafletLength)))%*% chol(MatG_leafletLength)

    coef.leafletLength_simu=coef.relative_leaflet_length_mean + coef.leafletLength_sd

    lengthFirst=max(0,coef.leafletLength_simu[,'L0'])
    lengthLast=max(0,coef.leafletLength_simu[,'Lfin'])
    posLengthMax=coef.leafletLength_simu[,'Pos_Lmax']

    ########################
    #relative leaflet width#
    ########################
    diag(MatG_leafletWidth)=diag(MatG_leafletWidth)+epsilon
    coef.leafletWidth_sd<-matrix(data=rnorm(n=ncol(MatG_leafletWidth)),ncol=ncol(MatG_leafletWidth),dimnames=list(Progeny= Progeny,Effet=colnames(MatG_leafletWidth)))%*% chol(MatG_leafletWidth)

    coef.leafletWidth_simu= coef.relative_leaflet_Width_mean + coef.leafletWidth_sd

    widthFirst= max(0,coef.leafletWidth_simu[,'W0'])
    widthLast=max(0,coef.leafletWidth_simu[,'Wfin'])
    posWidthMax=coef.leafletWidth_simu[,'Pos_Wmax']

    ################################
    #parameters Leaflet axial angle#
    ################################
    diag(MatG_axialAngle)=diag(MatG_axialAngle)+10^-2
    coef.axialAngle_sd<-matrix(data=rnorm(n=ncol(MatG_axialAngle)),ncol=ncol(MatG_axialAngle),dimnames=list(Progeny= Progeny,Effet=rownames(MatG_axialAngle)))%*% chol(MatG_axialAngle)

    coef.axialAngle_simu= coef.axialAngle_mean + coef.axialAngle_sd

    leafletAxialAngleC= coef.axialAngle_simu[,'angleC']
    leafletAxialAngleA=coef.axialAngle_simu[,'angleA']
    leafletAxialAngleSlope=coef.axialAngle_simu[,'slopeC']
    leafletAxialAngle_SDP=SigmaR_axialAngle


    #############################
    #parameters width sectionC	#
    #############################

    diag(MatG_widthC)=diag(MatG_widthC)+ epsilon
    coef.widthC_sd<-matrix(data=rnorm(n=ncol(MatG_widthC)),ncol=ncol(MatG_widthC),dimnames=list(Progeny= Progeny,Effet=rownames(MatG_widthC)))%*% chol(MatG_widthC)

    coef.widthC_simu= coef.widthC_mean + coef.widthC_sd

    frondCpointWidthIntercept=coef.widthC_simu[,'(Intercept)']
    frondCpointWidthSlope=coef.widthC_simu[,'CtoA']


    ####################################
    #parameters rachis relative height #
    ####################################
    coef.rachisRelativeHeight_sd=rnorm(n=1,mean=0,sd=SdG1_rachisRelativeHeight)
    coef.rachisRelativeHeight_simu= coef.rachisRelativeHeight_mean + coef.rachisRelativeHeight_sd

    rachisHeightTappering=coef.rachisRelativeHeight_simu['a']

    # -----PARAMETERS_FILE-----

    paramNames=c(
      paste('long seed = '),
      paste('int nbLeafEmitted = '),

      paste('double frondPhyllo_M = '),
      paste('double frondPhyllo_SD = '),

      paste('double H0 = '),
      paste('double coefStemHeight = '),
      paste('double residStemHeight= '),

      paste('double trunkBending_M = '),
      paste('double trunkBending_SD = '),

      paste('int nbFronds_M = '),
      paste('int nbFronds_SD = '),

      paste('double stemDiamMax = '),
      paste('double stemDiamSlope = '),
      paste('double stemDiamInfl = '),
      paste('double residStemDiam = '),

      paste('double decliCintercept ='),
      paste('double decliCslope ='),
      paste('double cPointAngle_SDP ='),

      paste('double decMaxA ='),
      paste('double decSlopeA ='),
      paste('double decInflA ='),
      paste('double APointAngle_SDP ='),

      paste('double rachisTwistFinalAngle_M ='),
      paste('double rachisTwistFinalAngle_SD ='),
      paste('double rachisTwistCoef ='),

      paste('double coefCurve ='),
      #paste('double coefCurve_SD ='),

      paste('double rachisDevFinalAngle_M ='),
      paste('double rachisDevFinalAngle_SD ='),
      paste('double rachisDevCoef ='),

      paste('double petioleRachisRatio_M ='),
      paste('double petioleRachisRatio_SD ='),

      paste('double rachisLengthIntercept ='),
      paste('double rachisLengthSlope ='),
      paste('double rachisLength_SDP ='),

      paste('double laminaAngle = '),

      paste('double nbMax ='),
      paste('double nbSlope ='),
      paste('double nbInfl ='),
      paste('double nbLeaflets_SDP ='),

      paste('double coefDispo= '),
      paste('double Dispo_SDP= '),
      paste('double pointBrelativePosition_M= '),
      paste('double pointBrelativePosition_SD= '),

      paste('double lenfletLengthAtBIntercept = '),
      paste('double leafletLengthAtBSlope = '),

      paste('double lengthFirst = '),
      paste('double lengthLast = '),
      paste('double posLengthMax = '),

      paste('double widthFirst = '),
      paste('double widthLast = '),
      paste('double posWidthMax = '),

      paste('double bWidthIntercept = '),
      paste('double bWidthSlope ='),

      paste('double xm_intercept= '),
      paste('double xm_slope= '),
      paste('double ym_intercept= '),
      paste('double ym_slope= '),

      paste('double leafletAxialAngleC = '),
      paste('double leafletAxialAngleA = '),
      paste('double leafletAxialAngleSlope = '),
      paste('double leafletAxialAngle_SDP = '),

      paste('double leafletStiffness = ' ),
      paste('double leafletStiffness_SD ='),

      paste('double[] leafletFrequencyHigh ='),
      paste('double[] leafletFrequencyLow =' ),

      paste('int nbInflorescences ='),

      paste('double frondBaseWidth = '),
      paste('double frondCpointWidthIntercept = '),
      paste('double frondCpointWidthSlope = '),
      paste('double frondtipWidth = '),
      paste('double frondBaseHeight = '),
      paste('double frondTipHeight = '),
      paste('double heightWidthCRatio = '),
      paste('double rachisHeightTappering = '),

      paste('double leafletRadialHighA0Sup= '),
      paste('double leafletRadialHighAmaxSup= '),
      paste('double leafletRadialHighA0Inf= '),
      paste('double leafletRadialHighAmaxInf= '),
      paste('double leafletRadialLowA0Sup= '),
      paste('double leafletRadialLowAmaxSup= '),
      paste('double leafletRadialLowA0Inf= '),
      paste('double leafletRadialLowAmaxInf= ')
    )

    paramValue=c(
      seed,
      nbLeafEmitted,
      frondPhyllo_M,
      frondPhyllo_SD,

      H0,
      coefStemHeight,
      residStemHeight,

      trunkBending_M,
      trunkBending_SD,

      nbFronds_M,
      nbFronds_SD,

      stemDiamMax,
      stemDiamSlope,
      stemDiamInfl,
      residStemDiam,

      decliCintercept,
      decliCslope,
      cPointAngle_SDP,

      decMaxA,
      decSlopeA,
      decInflA,
      APointAngle_SDP,

      rachisTwistFinalAngle_M,
      rachisTwistFinalAngle_SD,
      rachisTwistCoef,

      coefCurve,
      #coefCurve_SD,

      rachisDevFinalAngle_M,
      rachisDevFinalAngle_SD,
      rachisDevCoef,

      petioleRachisRatio_M,
      petioleRachisRatio_SD,

      rachisLengthIntercept,
      rachisLengthSlope,
      rachisLength_SDP,

      laminaAngle,

      nbMax,
      nbSlope,
      nbInfl,
      nbLeaflets_SDP,

      coefDispo,
      Dispo_SDP,
      pointBrelativePosition_M,
      pointBrelativePosition_SD,

      lenfletLengthAtBIntercept,
      leafletLengthAtBSlope,

      lengthFirst,
      lengthLast,
      posLengthMax,

      widthFirst,
      widthLast,
      posWidthMax,

      bWidthIntercept,
      bWidthSlope,

      xm_intercept,
      xm_slope,
      ym_intercept,
      ym_slope,

      leafletAxialAngleC,
      leafletAxialAngleA,
      leafletAxialAngleSlope,
      leafletAxialAngle_SDP,

      leafletStiffness,
      leafletStiffness_SD,

      paste(list(leafletFrequencyHigh)),
      paste(list(leafletFrequencyLow)),

      nbInflorescences,

      frondBaseWidth,
      frondCpointWidthIntercept,
      frondCpointWidthSlope,
      frondtipWidth,
      frondBaseHeight,
      frondTipHeight,
      heightWidthCRatio,
      rachisHeightTappering,

      leafletRadialHighA0Sup,
      leafletRadialHighAmaxSup,
      leafletRadialHighA0Inf,
      leafletRadialHighAmaxInf,
      leafletRadialLowA0Sup,
      leafletRadialLowAmaxSup,
      leafletRadialLowA0Inf,
      leafletRadialLowAmaxInf
    )


    output=as.matrix(paramValue)
    rownames(output)= paramNames

    return(output)
  }

  ###Mean parameter file per progeny
  ArchiOutputMean=function(Progeny, nbLeafEmitted){
    #-----------------------------------FIXED_PARAMETERS-----------------------------------------------

    trunkBending_M=0
    trunkBending_SD=0

    # internodeFinalLength=2

    laminaAngle = 140

    leafletStiffness =    1500
    leafletStiffness_SD = 7000

    frondBaseWidth = 30
    frondtipWidth = 0.3
    frondBaseHeight = 10
    frondTipHeight = 0.5

    nbInflorescences = 0

    rachisTwistCoef=3
    rachisDevCoef=2

    #-------------------------------FIXED_PARAMETERS_PER_PROGENY---------------------------------------
    ##################
    #stem phylotaxis #
    ##################
    frondPhyllo_M= mean(Phylo$Phylo,na.rm=T)
    frondPhyllo_SD= sd(Phylo$Phylo,na.rm=T)

    ###############
    # stem height #
    ###############
    residStemHeight

    ##################
    #number of frond #
    ##################
    nbFronds_M
    nbFronds_SD=0


    #############
    #frond twist#
    #############
    rachisTwistFinalAngle_M=mean(Tor$TwistA,na.rm=T)
    rachisTwistFinalAngle_SD=sd(Tor$TwistA,na.rm=T)

    #################
    #frond deviation#
    #################
    rachisDevFinalAngle_M=mean(Dev$DevA_deg)
    rachisDevFinalAngle_SD=sd(Dev$DevA_deg)

    #######################################
    #parameters ratio petiole/rachis final#
    #######################################
    petioleRachisRatio_M=ifelse(is.na(mean(unique(Pet$RatioPetiole),na.rm=T)),mean(DataAll[DataAll$Progeny== Progeny & DataAll$FrondRank>17 ,]$RatioPetiole,na.rm=T),mean(unique(Pet$RatioPetiole),na.rm=T))
    petioleRachisRatio_SD=ifelse(is.na(sd(unique(Pet$RatioPetiole),na.rm=T)),sd(DataAll[DataAll$Progeny== Progeny & DataAll$FrondRank>17 ,]$RatioPetiole,na.rm=T),sd(unique(Pet$RatioPetiole),na.rm=T))

    #####################################
    #B point relative position on rachis#
    #####################################
    pointBrelativePosition_M=ifelse(is.na(mean(Bpos$PosB,na.rm=T)),mean(DataAll[DataAll$Progeny== Progeny,]$PosB,na.rm=T),mean(Bpos$PosB,na.rm=T))
    pointBrelativePosition_SD=ifelse(is.na(sd(Bpos$PosB,na.rm=T)),sd(DataAll[DataAll$Progeny== Progeny,]$PosB,na.rm=T),sd(Bpos$PosB,na.rm=T))


    ################
    #leaflet number#
    ################
    nbMax= coef.nbLeaflets_mean['nbMax']
    nbSlope = coef.nbLeaflets_mean['nbSlope']
    nbInfl= coef.nbLeaflets_mean['nbInfl']
    nbLeaflets_SDP=round(SigmaR_nbLeaflets)

    #################################
    #parameters Leaflet radial angle#
    #################################
    leafletRadialHighA0Sup=coef(radialHighSup.nls)['A0']
    leafletRadialHighAmaxSup=coef(radialHighSup.nls)['Amax']

    leafletRadialHighA0Inf=coef(radialHighInf.nls)['A0']
    leafletRadialHighAmaxInf=coef(radialHighInf.nls)['Amax']

    leafletRadialLowA0Sup=coef(radialLowSup.nls)['A0']
    leafletRadialLowAmaxSup=coef(radialLowSup.nls)['Amax']

    leafletRadialLowA0Inf=coef(radialLowInf.nls)['A0']
    leafletRadialLowAmaxInf=coef(radialLowInf.nls)['Amax']

    ############################
    #Frequency of leaflets type#
    ############################
    leafletFrequencyHigh=Rep[Rep$Type==1,]$Prop
    leafletFrequencyLow=Rep[Rep$Type==-1,]$Prop

    ###############
    #leaflet shape#
    ###############
    xm_intercept=coef(xm.lm)['(Intercept)']
    xm_slope=coef(xm.lm)['PositionRelative']
    ym_intercept=coef(ym.lm)['(Intercept)']
    ym_slope=coef(ym.lm)['PositionRelative']

    #################################
    #ratio height/width sectionC	#
    #################################
    heightWidthCRatio=mean(PetioleSectionC$Petiole_height_C_cm/PetioleSectionC$Petiole_width_C_cm)


    #------------------------------------PARAMETERS_PER_TREE---------------------------------------------

    ######
    #Seed#
    ######
    seed=0
    set.seed(seed)

    ########
    #LeafNb#
    ########

    nbLeafEmitted

    ###############
    #rachis length#
    ###############

    coef.rachisLength_simu= coef.rachisLength_mean

    rachisLengthIntercept=coef.rachisLength_simu['(Intercept)']
    rachisLengthSlope=coef.rachisLength_simu['LeafNumber']
    rachisLength_SDP=SigmaR_rachisLength

    ################
    #stem diameter #
    ################

    stemDiamMax=coef(StemDiam.nls)['finalStemDiam']
    stemDiamSlope= coef(StemDiam.nls)['StemDiamSlope']
    stemDiamInfl= coef(StemDiam.nls)['StemDiamInfl']


    #####################
    #C angle declination#
    #####################

    coef.decliC_simu= coef.decliC_mean

    decliCintercept =coef.decliC_simu['(Intercept)']
    decliCslope =coef.decliC_simu['Rank']

    cPointAngle_SDP = SigmaR_decliC

    ########################
    #declination frond tip #
    ########################
    decMaxA
    decSlopeA
    decInflA

    APointAngle_SDP = SigmaR_decliA

    #####################
    #rachis curvature	#
    #####################

    coefCurve=rnorm(n=1,mean=mean(paramTree$coefCurv),sd=sd(paramTree$coefCurv))

    ###################
    #leaflets position#
    ###################

    coefDispo
    Dispo_SDP= residDispo

    ########################
    #leaflet length b point#
    ########################
    coef.length_B_simu=coef.length_B_mean

    lenfletLengthAtBIntercept=coef.length_B_simu['(Intercept)']
    leafletLengthAtBSlope=coef.length_B_simu['RachisLength']

    ###################################
    #parameters  leaflet width b point#
    ###################################
    coef.width_B_simu= coef.width_B_mean

    bWidthIntercept= coef.width_B_simu['(Intercept)']
    bWidthSlope =coef.width_B_simu['RachisLength']

    ##########################
    # relative leaflet length#
    ##########################
    coef.leafletLength_simu=coef.relative_leaflet_length_mean

    lengthFirst=max(0,coef.leafletLength_simu['L0'])
    lengthLast=max(0,coef.leafletLength_simu['Lfin'])
    posLengthMax=coef.leafletLength_simu['Pos_Lmax']

    ########################
    #relative leaflet width#
    ########################
    coef.leafletWidth_simu= coef.relative_leaflet_Width_mean

    widthFirst= max(0,coef.leafletWidth_simu['W0'])
    widthLast=max(0,coef.leafletWidth_simu['Wfin'])
    posWidthMax=coef.leafletWidth_simu['Pos_Wmax']

    ################################
    #parameters Leaflet axial angle#
    ################################
    coef.axialAngle_simu= coef.axialAngle_mean

    leafletAxialAngleC= coef.axialAngle_simu['angleC']
    leafletAxialAngleA=coef.axialAngle_simu['angleA']
    leafletAxialAngleSlope=coef.axialAngle_simu['slopeC']
    leafletAxialAngle_SDP=SigmaR_axialAngle


    #############################
    #parameters width sectionC	#
    #############################

    coef.widthC_simu= coef.widthC_mean
    frondCpointWidthIntercept=coef.widthC_simu['(Intercept)']
    frondCpointWidthSlope=coef.widthC_simu['CtoA']


    ####################################
    #parameters rachis relative height #
    ####################################
    coef.rachisRelativeHeight_sd=rnorm(n=1,mean=0,sd=SdG1_rachisRelativeHeight)
    coef.rachisRelativeHeight_simu= coef.rachisRelativeHeight_mean + coef.rachisRelativeHeight_sd

    rachisHeightTappering=coef.rachisRelativeHeight_simu['a']

    #_______________________________________________________________________________________________________________________________________________PARAMETERS_FILE___________________________________________________________________________________________________________________________________________________________#

    paramNames=c(
      paste('long seed = '),
      paste('int nbLeafEmitted = '),

      paste('double frondPhyllo_M = '),
      paste('double frondPhyllo_SD = '),

      paste('double H0 = '),
      paste('double coefStemHeight = '),
      paste('double residStemHeight= '),

      paste('double trunkBending_M = '),
      paste('double trunkBending_SD = '),

      paste('int nbFronds_M = '),
      paste('int nbFronds_SD = '),

      paste('double stemDiamMax = '),
      paste('double stemDiamSlope = '),
      paste('double stemDiamInfl = '),
      paste('double residStemDiam = '),

      paste('double decliCintercept ='),
      paste('double decliCslope ='),
      paste('double cPointAngle_SDP ='),

      paste('double decMaxA ='),
      paste('double decSlopeA ='),
      paste('double decInflA ='),
      paste('double APointAngle_SDP ='),

      paste('double rachisTwistFinalAngle_M ='),
      paste('double rachisTwistFinalAngle_SD ='),
      paste('double rachisTwistCoef ='),

      paste('double coefCurve ='),
      #paste('double coefCurve_SD ='),

      paste('double rachisDevFinalAngle_M ='),
      paste('double rachisDevFinalAngle_SD ='),
      paste('double rachisDevCoef ='),

      paste('double petioleRachisRatio_M ='),
      paste('double petioleRachisRatio_SD ='),

      paste('double rachisLengthIntercept ='),
      paste('double rachisLengthSlope ='),
      paste('double rachisLength_SDP ='),

      paste('double laminaAngle = '),

      paste('double nbMax ='),
      paste('double nbSlope ='),
      paste('double nbInfl ='),
      paste('double nbLeaflets_SDP ='),

      paste('double coefDispo= '),
      paste('double Dispo_SDP= '),
      paste('double pointBrelativePosition_M= '),
      paste('double pointBrelativePosition_SD= '),

      paste('double lenfletLengthAtBIntercept = '),
      paste('double leafletLengthAtBSlope = '),

      paste('double lengthFirst = '),
      paste('double lengthLast = '),
      paste('double posLengthMax = '),

      paste('double widthFirst = '),
      paste('double widthLast = '),
      paste('double posWidthMax = '),

      paste('double bWidthIntercept = '),
      paste('double bWidthSlope ='),

      paste('double xm_intercept= '),
      paste('double xm_slope= '),
      paste('double ym_intercept= '),
      paste('double ym_slope= '),

      paste('double leafletAxialAngleC = '),
      paste('double leafletAxialAngleA = '),
      paste('double leafletAxialAngleSlope = '),
      paste('double leafletAxialAngle_SDP = '),

      paste('double leafletStiffness = ' ),
      paste('double leafletStiffness_SD ='),

      paste('double[] leafletFrequencyHigh ='),
      paste('double[] leafletFrequencyLow =' ),

      paste('int nbInflorescences ='),

      paste('double frondBaseWidth = '),
      paste('double frondCpointWidthIntercept = '),
      paste('double frondCpointWidthSlope = '),
      paste('double frondtipWidth = '),
      paste('double frondBaseHeight = '),
      paste('double frondTipHeight = '),
      paste('double heightWidthCRatio = '),
      paste('double rachisHeightTappering = '),

      paste('double leafletRadialHighA0Sup= '),
      paste('double leafletRadialHighAmaxSup= '),
      paste('double leafletRadialHighA0Inf= '),
      paste('double leafletRadialHighAmaxInf= '),
      paste('double leafletRadialLowA0Sup= '),
      paste('double leafletRadialLowAmaxSup= '),
      paste('double leafletRadialLowA0Inf= '),
      paste('double leafletRadialLowAmaxInf= ')
    )

    paramValue=c(
      seed,
      nbLeafEmitted,
      frondPhyllo_M,
      frondPhyllo_SD,

      H0,
      coefStemHeight,
      residStemHeight,

      trunkBending_M,
      trunkBending_SD,

      nbFronds_M,
      nbFronds_SD,

      stemDiamMax,
      stemDiamSlope,
      stemDiamInfl,
      residStemDiam,

      decliCintercept,
      decliCslope,
      cPointAngle_SDP,

      decMaxA,
      decSlopeA,
      decInflA,
      APointAngle_SDP,

      rachisTwistFinalAngle_M,
      rachisTwistFinalAngle_SD,
      rachisTwistCoef,

      coefCurve,
      #coefCurve_SD,

      rachisDevFinalAngle_M,
      rachisDevFinalAngle_SD,
      rachisDevCoef,

      petioleRachisRatio_M,
      petioleRachisRatio_SD,

      rachisLengthIntercept,
      rachisLengthSlope,
      rachisLength_SDP,

      laminaAngle,

      nbMax,
      nbSlope,
      nbInfl,
      nbLeaflets_SDP,

      coefDispo,
      Dispo_SDP,
      pointBrelativePosition_M,
      pointBrelativePosition_SD,

      lenfletLengthAtBIntercept,
      leafletLengthAtBSlope,

      lengthFirst,
      lengthLast,
      posLengthMax,

      widthFirst,
      widthLast,
      posWidthMax,

      bWidthIntercept,
      bWidthSlope,

      xm_intercept,
      xm_slope,
      ym_intercept,
      ym_slope,

      leafletAxialAngleC,
      leafletAxialAngleA,
      leafletAxialAngleSlope,
      leafletAxialAngle_SDP,

      leafletStiffness,
      leafletStiffness_SD,

      paste(list(leafletFrequencyHigh)),
      paste(list(leafletFrequencyLow)),

      nbInflorescences,

      frondBaseWidth,
      frondCpointWidthIntercept,
      frondCpointWidthSlope,
      frondtipWidth,
      frondBaseHeight,
      frondTipHeight,
      heightWidthCRatio,
      rachisHeightTappering,

      leafletRadialHighA0Sup,
      leafletRadialHighAmaxSup,
      leafletRadialHighA0Inf,
      leafletRadialHighAmaxInf,
      leafletRadialLowA0Sup,
      leafletRadialLowAmaxSup,
      leafletRadialLowA0Inf,
      leafletRadialLowAmaxInf
    )


    output=as.matrix(paramValue)
    rownames(output)= paramNames

    return(output)
  }


  #____________PARAMETERS_FILE_EXPORTATION_______________#####

  #set working directory
  setwd(parametersFileDirectory)

  ###mean progeny
  filenamePro= file.path("Outputs/ParametersFiles",paste(progeny,'_AverageTree_', map,"MAP.txt",sep = ""))
  capture.output(ArchiOutputMean(Progeny=progeny, nbLeafEmitted = nbLeafEmitted),file = filenamePro)

  ###individuals
  for (i in (1:20)){
    filename= file.path("Outputs/ParametersFiles",paste(progeny,'_','Tree',i,'_', map,"MAP.txt",sep = ""))
    capture.output(ArchiOutput(Progeny=progeny, nbLeafEmitted = nbLeafEmitted,seed=i),file = filename)
  }
}
