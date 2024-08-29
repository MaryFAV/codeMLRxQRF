############################################
# MODELE MLR POUR LA DETECTION DES ARGILES #
# FAVIERE MARYLINE - STAGE UMR LISAH #######
# 2024 #####################################
############################################


################
# PREPARATION ##
################

rm(list=ls())
graphics.off()

# LIBRAIRIES #
library(pls) 
library(raster)       
library(rgdal)      
library(caret)
library(e1071)
library(ade4)         
library(gstat)
library(snow)
library(quantregForest)
library(tibble)
library(rasterVis)
library(sf)
library(dplyr)
data(SAGA_pal)


################
#   ENTREES    #
################
# ENREGISTREMENT HEURE SYSTEME AU DEBUT EXECUTION # 
temps_debut <- Sys.time()

# INITIALISATION DE LA VARIABLE QUI EMPILERA TOUTES LES CARTES EN FONCTION DU NOMBRE ITERATION # 
FinalMap <- stack() 

# INITIALISATION DE LA LISTE DES ERREURS # 
metrics_list <- list() 

# MATRICE QUU CONTIENT DES PREDICTIONS DE VALIDATION # 
NbTest = 36
NbT = 23
iteration = 100

SaveMOPredictionVal <- matrix(0, nrow = NbTest, ncol = iteration)
SaveMOPredictionTest <- matrix(0, nrow = NbT, ncol = iteration)

pathS2="C:/Users/maryf/Documents/code/data_versailles/"
#pathS2="//147.99.18.48/shared-data/Maryline_Faviere/Data_Inde/1_Data/"

# TROUVER LES MESURES #
path="C:/Users/maryf/Documents/code/data_versailles/"  
#path="//147.99.18.48/shared-data/Maryline_Faviere/Data_Inde/1_Data/"

# CHEMIN DE SORTIE #
#pathResults="Y:/Maryline_Faviere/Résultats/MLR/Inde_20170424/"
pathResults="C:/Users/maryf/Documents/CARTES_STAGES/Versailles_CO/BTS/MLR/"

#pathResults="//147.99.18.48/shared-data/Maryline_Faviere/res/"

# DISTANCE DE MAHALANOBIS #
seuilMD=3.5 

# DATE IMAGE S2 #
DateImage="BGP3"

Property="Carorg"
prop_unite="g.kg-1"

# VALEUR MAX DES PLOTS #
LimMaxscatterplot=80 

# Lire le fichier CSV en spécifiant le séparateur décimal
PropSoil_NBSS <- read.table(paste(path, "placettes_CO.csv", sep = ""), 
                            row.names = 1, 
                            head = TRUE, 
                            sep = ",", 
                            dec = ".")

# Convertir les valeurs en valeurs numériques
for (i in 1:ncol(PropSoil_NBSS)) {
  PropSoil_NBSS[, i] <- as.numeric(as.character(PropSoil_NBSS[, i]))
}

# Afficher la structure du dataframe
str(PropSoil_NBSS)

# Transposer le dataframe
PropSoil_NBSS_t <- t(PropSoil_NBSS)

# Supprimer les lignes contenant NA
PropSoil_NBSS_t_clean <- PropSoil_NBSS_t[complete.cases(PropSoil_NBSS_t), ]

# Transposer à nouveau pour obtenir le dataframe original sans les colonnes contenant NA
PropSoil_NBSS <- t(PropSoil_NBSS_t_clean)

# Convertir le résultat en dataframe
PropSoil_NBSS <- as.data.frame(PropSoil_NBSS)

dim(PropSoil_NBSS) # AFFICHAGE DES DIMENSIONS (LIGNES ET COLONNES) #
head(PropSoil_NBSS)# AFFICHAGE DES PREMIERES LIGNES DU TABLEAU #
str(PropSoil_NBSS) # AFFICHAGE DE LA STRUCTURE DE L'OBJET #                   

# EXTRACTION DE LA COLONNE PROPERTY # 
Prop=PropSoil_NBSS[Property]                            

# CONVERSION DES DONNEES EN DONNEES SPATIALES #
coordinates(PropSoil_NBSS) <- c("Direct_N","Direct_E")
coord <- coordinates(PropSoil_NBSS)

SavePerfo=matrix(c(0.0),1,14) # CREATION DE LA MATRICE POUR STOCKER LES PERFORMANCES DU MODELE # 
SavePerfo=as.data.frame(SavePerfo) # CONVERSION DE LA MATRICE EN DATAFRAME # 
colnames(SavePerfo)=c("nbrData", "nbrOutlier","R2cal","RMSEcal","R2val","RMSEval","RPDval","RPIQval","Biasval", "R2Test", "RMSETest", "RPDTest", "RPIQTest", "BiasTest") # ATTRIBUTION DE NOMS AUX COLONNES # 

# OUVERTURE DE SE #
setwd(pathS2) # DEFINIR LE REPERTOIRE ACTUEL #
NameImage=c(paste("mosaiq_mask_",DateImage,sep="")) # CREATION VECTEUR POUR CONTENIR DE NOM DE IMAGE CONCATENE "BARESOIL_S2_" #
rownames(SavePerfo)=paste(DateImage) # AFFECTION DU NOM ET DATE SPECIFIE COMME NOM DE LIGNE DANS SAVEPERFO #

namefile=paste(NameImage,".tif",sep="") # CONCATENATION DU NOM DE FICHIER #
S2 = readGDAL(paste(pathS2,namefile, sep = "")) # LECTURE DES IMAGES RASTER AVEC LA FONCTION READGDAL() # 
S2Band.bare=stack(S2) # EMPILEMENT DE BANDES AVEC LA FONCTION STACK() #


################################
#   PREPARATION DES DONNEES    #
################################

# EXTRACTION DES DONNEES DES PILES DE RASTER # 
Spectres=extract(S2Band.bare, PropSoil_NBSS) # EXTRACTION DES VALEURS SPECTRALES DES PIXELS DE IMAGE DES POINTS DE SOL DANS PROPSOIL_NBSS # 
MO_spectres=data.frame(Spectres) # CREATION DATAFRAME A PARTIR DES VALEURS SPECTRALES EXTRAITES STOCKEE DANS LA VARIABLE PROP #

MO_spectres$Prop=Prop # AJOUTE UNE COLONNE PROP AU DATAFRAME CONTENANT DES VALEURS DE PROPRIETES DU SOL DE LA VARIABLE PROP # 
MO_spectres$coord=coord # AJOUTE UNE COLONN COORD AU DATAFRAME CONTENANT LES COORDONNEES XY DES ECHANTILLONS # 

# REORGANISATION DES COLONNES POUR QUE LES PROPRIETES DU SOL SOIENT LA PREMIERE COLONNE ECT # 
MO_spectresTemp=MO_spectres 
MO_spectresTemp[,1]=MO_spectres[,11]
MO_spectresTemp[,2:11]=MO_spectres[,1:10]

colnames(MO_spectresTemp)[1]=Property
colnames(MO_spectresTemp)[2:11]=colnames(MO_spectres)[1:10]
rownames(MO_spectresTemp)=rownames(Prop)

MO_spectres=MO_spectresTemp

# CONVERTIR LES NOMS DES ECHANTILLONS EN VALEUR # 
MO_spectres <- rownames_to_column(MO_spectres, var = "ID")

# AFFICHAGE DES DIMENSIONS DU DATAFRAME ET AFFICHAGE DES PREMIERES LIGNES # 
print(dim(MO_spectres))
head(MO_spectres)

##############################################
#   SUPPRESSION DES DONNEES DE VEGETATION    #
##############################################
spectresFalseTrue<-is.na(Spectres) # CREATION DE LA MATRICE OU CHAQUE ELEMENT EST TRUE SI LA VALEUR DES SPECTRES EST NA, ET FALSE SINON # 
l=1    # INITIALISATION DU COMPTEUR A 1 # 
MO_spectres_SansNA=MO_spectres    # CREATION COPIE DE SAUVEGARDE DU DATAFRAME # 

# BOUCLE QUI PARCOURT CHAQUE LIGNE DE PROP #
# SI LA VALEUR DE SELECT EST FALSE, ALORS LIGNE EST COPIEE DANS BDBERAMBADISANSNAN A POSITION l, ET INCRIMENTATION DE 1 # 
for (i in 1:dim(Prop)[1]) {
  if (spectresFalseTrue[i,1]==FALSE) {
    MO_spectres_SansNA[l,]<-MO_spectres[i,]
    l=l+1  
  }
}

# LES LIGNES NA SONT SUPPRIME EN PRENANT QUE LES LIGNES JUSQUA POSITION l-1. 
MO_spectres_SansNA=MO_spectres_SansNA[2:l-1,]
head(MO_spectres_SansNA)

# Définir les noms des échantillons comme noms de lignes
rownames(MO_spectres_SansNA) <- MO_spectres_SansNA$ID


########################################
#   PREPARATION DES DONNEES DE TEST    #
########################################

nb_total_samples <- dim(MO_spectres_SansNA)[1] # DETERMINATION DU NOMBRE TOTAL ECHANTILLONS SOL #
SavePerfo[1] <- nb_total_samples # SAUVEGARDE DU NOMBRE ECHANTILLONS # 
matrice_MO_sorted <- MO_spectres_SansNA[order(MO_spectres_SansNA[, "Carorg"]),, drop = FALSE] # TRI DES VALEURS MOS EN ORDRE CROISSANT # 
selected_values <- matrix(nrow = 0, ncol = ncol(matrice_MO_sorted)) # INITIALISATION DE LA MATRICE POUR STOCJER LES VALEURS SELECTIONNEES # 

# SELECTTION DES VALEURS TOUTES LES 7 VALEURS PAR ORDRE CROISSANT # 
for (i in seq(7, nb_total_samples, by = 7)) {
  selected_values <- rbind(selected_values, matrice_MO_sorted[i, ])
}

# Définir les noms des échantillons comme noms de lignes
rownames(selected_values) <- selected_values$ID

selected_values=t(selected_values) # INVERSEMENT DES LIGNES ET COLONNES #
selected_values=t(selected_values) # INVERSEMENT DES LIGNES ET COLONNES #

ExtractNameTest <-as.matrix(selected_values[,1]) # EXTRACTION NOMS DE LIGNE 1 # 
selected_values <- selected_values[, -1]

# Appliquer as.numeric à chaque élément de la matrice
selected_values <- apply(selected_values, c(1, 2), as.numeric)

# Reconvertir la matrice en caractères
selected_values <- as.matrix(selected_values)

selected_indices <- seq(7, nb_total_samples, by = 7) # INDICES DE LIGNES NON SELECTIONNEES POUR DATA TEST # 
non_selected_indices <- setdiff(seq(nb_total_samples), selected_indices) # DIFFERENCE ENTRE LE TOTAL ET INDICES NON SELECTIONNES #
non_selected_values <- matrice_MO_sorted[non_selected_indices, ] # RECUPERATION DES VALEURS DE LIGNE ENTIERE # 
rownames(non_selected_values) <- non_selected_values$ID
non_selected_values=t(non_selected_values) # INVERSEMENT DES LIGNES ET COLONNES #
non_selected_values=t(non_selected_values) # INVERSEMENT DES LIGNES ET COLONNES #

non_selected_values <- non_selected_values[, -1]
# Appliquer as.numeric à chaque élément de la matrice
non_selected_values <- apply(non_selected_values, c(1, 2), as.numeric)

# Reconvertir la matrice en caractères
non_selected_values <- as.matrix(non_selected_values)

ExtractProp_Test<-as.matrix(selected_values[,1]) # ETRACTION DES PROPRIETE DU SOL DE LA 1ERE LIGNE DE BDBERAMBADISANSNAN ET CONVERSION DANS MATRICE # 
ExtractProp_Test=t(ExtractProp_Test)
ExtractSPectres_Test<-as.matrix(selected_values[,2:11]) # EXTRACTIONS SPETRES DES LIGNES 2 A 11 DE BDBERAMBADISANSNAN # 

VTEST <- as.data.frame(ExtractSPectres_Test) # CREATION DATAFRAME # 
SpectACPyTest<-t(scale(VTEST, center = FALSE, scale = FALSE)) # CENTRAGE DES VALEURS SPECTRO # 
SpectACPyTest=t(SpectACPyTest)

#ExtractProp_Test=t(ExtractProp_Test) # INVERSEMENT DES LIGNES ET COLONNES #
ProprieteSolTest<-ExtractProp_Test[1,] # CONVERSION EN NUMERIQUE # 
DataTest<-data.frame(cbind(ProprieteSolTest, SpectACPyTest)) # FUSION DES PROPRIETES SOL ET SPECTRO # 
colnames(DataTest)[1] <- "ProprieteSol"

########################################
#   PREPARATION DE TEST / TRAIN DB    #
#######################################

non_selected_values=t(non_selected_values) # INVERSEMENT DES LIGNES ET COLONNES #

ExtractProp_CALVAL<-as.matrix(non_selected_values[1,]) # ETRACTION DES PROPRIETE DU SOL DE LA 1ERE LIGNE ET CONVERSION DANS MATRICE # 
ExtractSPectres_CALVAL<-as.matrix(non_selected_values[2:11,]) # EXTRACTIONS SPETRES DES LIGNES 2 A 11 # 
ExtractSPectres_CALVAL=t(ExtractSPectres_CALVAL)
non_selected_values=t(non_selected_values) # INVERSEMENT DES LIGNES ET COLONNES #

# CALCUL DU NOMBRE DE DONNEES DE VALIDATION ET DE CALIBRATION # 
(NbrDataValid <- trunc(dim(non_selected_values)[1]/4))   # NOMBRE DE DONNEES DE VALIDATION EST CALCULE EN DIVISANT PAR 4 LE NOMBRE TOT DE COLONNES #         
(NbrDataCalib <- dim(non_selected_values)[1] - NbrDataValid) # NOMBRE DE DONNEES DE CALIBRATION EST CALCULE COMME LA DIFF ENTRE LE NOMBRE TOT DE COLONNES ET NBR DE DONNEES DE VALIDATION # 
NbrDataTest <- dim(selected_values)[1] # NOMBRE DE DONNEES TEST #

DF_selected_values <- as.data.frame(selected_values) # CONVERSION DE SELECTED_VALUES EN DATAFRAME #
sf_data <- st_as_sf(DF_selected_values, coords = c("coord.Direct_N", "coord.Direct_E"), remove = FALSE) # CREATION DE OBJET SF EN CONSERVANT LES COORDONNEES  XY #
st_crs(sf_data) <- st_crs(32631) # DEFINITION DU SCR # 
write.csv(sf_data, file = paste0(pathResults, "MLR__Points_Test.csv"), row.names = FALSE) # EXPORT DES POINTS TEST #

# DEPART DE LA BOUCLE #
iteration <- 100   # LE NOMBRE ITERATIONS VOULUES # 
while (iteration > 0) {
  num <- iteration   # INDICATION DU NUMERO DE LA BOUCLE ACTUELLE DANS CERTAINES DONNEES EN SORTIE # 
  
  tirages_CALVAL <- sample(dim(non_selected_values)[1]) # TIRAGE ALEATOIRE DES DONNEES CAL / VAL # 
  
  non_selected_values=t(non_selected_values) # INVERSEMENT DES LIGNES ET COLONNES #
  
  # INITIALISATION DES ENSEMBLES DE DONNEES DE CALIBRATION ET DE VALIDATION # 
  # Réinitialiser la matrice SpectreBerambadiCalib
  SpectreBerambadiCalib <- NULL
  
  SpectreBerambadiCalib <- non_selected_values[, tirages_CALVAL[1:NbrDataCalib]]
  SpectreBerambadiValid <- non_selected_values[, tirages_CALVAL[(NbrDataCalib + 1):length(tirages_CALVAL)]]
  
  # EXTRCATION DES DONNEES SPECTRALES DE CALIB DE LA MATRICE # 
  SpectModifCalib <- SpectreBerambadiCalib[2:11,]      # CHAQUE LIGNE REPRESENTE UN ECHANTILLON ET CHAQUE COLONNE UNE LONGUEUR ONDE OU BANDE SPECTRALE # 
  
  # CONVERSION DE LA MATRICE SPECTMODIFCALIB EN OABJET DE TYPE DATAFRAME (VCAL) 
  vCAL <- as.data.frame(SpectModifCalib)
  
  # SAUVEGARDE DES POINTS DE CALIBRATION EN CSV # 
  SpectreBerambadiCalib=t(SpectreBerambadiCalib)
  
  DF_Calib_values <- as.data.frame(SpectreBerambadiCalib) # CONVERSION DE SELECTED_VALUES EN DATAFRAME #
  sf_data_Calib <- st_as_sf(DF_Calib_values, coords = c("coord.Direct_N", "coord.Direct_E"), remove = FALSE) # CREATION DE OBJET SF EN CONSERVANT LES COORDONNEES  XY #
  st_crs(sf_data_Calib) <- st_crs(32631) # DEFINITION DU SCR # 
  write.csv(sf_data_Calib, file = paste0(pathResults, "MLR_", num, "Points_Calibration.csv"), row.names = FALSE)
  
  #####################################
  #   RECHERCHE OUTLIERS SPECTRAUX    #
  #####################################
  # ANALYSE EN COMPOSANTES PRINCIPALES (PCA) SUR LES DONNEES SPECTRALES [VCAL] 
  dudinir <- dudi.pca(t(vCAL),center = TRUE, scale=TRUE,scannf = F) # DUDI.PCA PREND EN ENTREE UNE MATRICE DE DONNEES OU LES COLONNES REPRESENTENT LES VARIABLES ET LES LIGNES DES OBSERVATIONS #
  
  # CALCUL DU % DE VARIANCE EXPLIQUEE  (PVE) POUR CHAQUE COMPOSANTE PRINCIPALE DE PCA # 
  pve <- 100 * dudinir$eig/sum(dudinir$eig) # DUDINIR$EIG CONTIENT LES VALEURS PROPRES DES COMPOSANTES PRINCIPALES # 
  
  # DETERMINATION DU NBR DE COMPOSANTES ACP # 
  out<-as.double(0)  # INITIALISE UNE VARIABLE OUT A 0 POUR CALCULER LE POURCENTAGE TOT DE LA VARIANCE EXPLIQUEE # 
  nf<-0   # INITALISE LE NBR DE COMPOSANTES A 0 # 
  for (i in 1:100){  # BOUCLE SUR LES COMPOSANTES PRINCIPALES POUR CALCULER LE NBR DE COMPOSANTES NECESSAIRES POUR EXPLIQUER 100% DE LA VARIANCE # 
    out<-as.double(out) + as.double(pve)[i]  # AJOUTE LE PVE DE CHAQUE COMPOSANTE A OUT # 
    nf<-nf+1   # INCREMENTE LE NBR DE COMPOSANTES # 
  }
  out # AFFICHE LE % TOT DE VARIANCE EXPLIQUEE # 
  nf  # AFFICHE LE NBR DE COMPOSANTES NECESSAIRES POUR EXPLIQUER 100% DE LA VARIANCE # 
  
  nf<-3  # FIXE LE NBR DE COMPOSANTES A UTILSER POUR ANALYSE A 3 # 
  
  # ANCIEN CALCUL DE LA DISTANCE DE MAHALANOBIS # 
  S <-  cov(dudinir$tab[,1:nf])  # CALCULE LA MATRICE DE COVARIANCE DES COMPOSANTES PRINCIPALES SELECTIONNEES # 
  det(S)  # CALCULE LE DETERMINANT DE LA MATRICE DE COVARIANCE S # 
  D2<-sqrt(mahalanobis(dudinir$tab[,1:nf],colMeans(dudinir$tab[,1:nf]) , S))   # CALCULE LA DISTANCE DE MAHALANOBIS POUR CHAQUE OBSERVATION / CENTRE DE DISTRIBUTION DES COMPOSANTES PRINCIPALES # 
  
  # CALCULE DU NBR OUTLIERS SPECTRAUX # 
  nbroutlier<-as.double(0) # INITIALISE LE NBR OUTLIERS A 0 # 
  NumEchant<-as.double(0)  # INITALISE UN VECTEUR POUR STOCKER LES INDICES DES OUTLIERS # 
  NomOutlier<-as.double(0)  # INITIALISE UN VECTEUR POUR STOCKER LES NOMS DES OUTLIERS # 
  for (i in 1:dim(vCAL)[2]){   # BOUCLE SUR CHAQUE OBSERVATION POUR IDENTIFIER LES OUTLIERS # 
    if (as.double(D2)[i] > seuilMD)  {   # VERIFICATION SI LA DM DEPASSE LE SEUIL DEFINI # 
      nbroutlier<-nbroutlier+1     # INCREMENTATION DU NBR OUTLIERS # 
      NumEchant[nbroutlier]<-i  # STOCKAGE INDICE OUTLIERS # 
      NomOutlier[nbroutlier]<-rownames(SpectModifCalib)[i]   # STOCKAGE DU NOM DE OUTLIERS # 
      D2[i]<-NA   # MARQUE LA DM DE OUTLIERS COMME NA POUR EXCLURE DE ANALYSE # 
    }
  }
  
  select<-is.na(D2)  # CREATION VECTEUR BOOLEENS POUR SELECTIONNNER OBSERVATIONS QUI NE SONT PAS OUTLIERS # 
  SpectreBerambadiCalib=t(SpectreBerambadiCalib)
  
  SpectreBerambadiCalib2<-SpectreBerambadiCalib[,!select]  
  
  nbroutlier  # AFFICHAGE DU NBR TOT DES OUTLIERS # 
  NumEchant  # AFFICHAGE DES INDICES DES OUTLIERS # 
  D2[NumEchant]  # AFFICHAGE DES DM DES OUTLIERS # 
  
  SavePerfo[2]=nbroutlier  # STOCKAGE DU NBR OUTLIERS IDENTIFIES DANS SAVEPERFO[2] # 
  
  vCAL2<-vCAL[!select]  # FILTRE DES DONNEES DE VCAL EN FONCTION DU VECTEUR SELECT #
  SpectACP<-t(scale(vCAL2, center = FALSE, scale = FALSE)) # CENTRE ET MET A ECHELLE DONNEES DANS VCAL2 # 
  
  NbrDataCalib2<-dim(vCAL2)[2] # CALCULE LE NBR DE VARIABLES DANS VCAL2 EN ACCEDANT A LA DIMENSION 2 DE LA MATRICE VACL2 # 
  
  # MISE DES MESURES MOS DANS LE VECTEUR PROPRIETESOL # 
  ProprieteSolCalib<-SpectreBerambadiCalib2[1,]  # EXTRACTION DES MESURES DE PROPRIETES DU SOL DE SPECTREBERAMBADICALIB2 ET LES ASSIGE A VARIABLE PROPRIETESOLCALIB # 
  
  DataCalib<-data.frame(cbind(ProprieteSolCalib, SpectACP)) # CREATION DATAFRAME EN COMBINANT LES MESURES DE LA PROPRIETE D SOL ET DONNEES SPECTRALES MISES A ECHELLE # 
  colnames(DataCalib)[1]=c("ProprieteSol") # RENOMMAGE DE LA 1ERE COLONNE DE DATACALIB EN PROPRIETESOL # 
  
  ##############################################
  #   PREPARATION DES DONNEES DE VALIDATION    #
  ##############################################
  
  # EXTRACTION DES VARIABLES NIR DE LA MATRIICE SPECTREBERAMBADIVALID EN EXCLUANT LA 1ERE LIGNE DES PROP DU SOL #
  SpectModifValid<-(SpectreBerambadiValid[2:11,])     
  
  vVAL<-as.data.frame(SpectModifValid) # CONVERSION DE LA MATRICE SPECTMODIFVALID EN DATARAME NOMME VVAL # 
  SpectACPy<-t(scale(vVAL, center = FALSE, scale = FALSE)) # CENTRE ET MET A ECHELLE DONNEES NIR DANS VVAL #  
  
  ProprieteSolValid<-SpectreBerambadiValid[1,] # EXTRACTION DES MESIRES DE LA PROP DU SOL DE LA 1ERE LIGNE DE MATRICE SPECTREBERAMBADIVALID ET LES ASSIGNE A LA VARIABLE PROPRIEESOLVALID # 
  
  DataValid<-data.frame(cbind(ProprieteSolValid, SpectACPy)) # CREATION DU DATAFRAME DATAVALID EN COMBINANT LES MESURES DE PROP DU SOL ET DONNEES NIR MISES A ECHELLE # 
  colnames(DataValid)=colnames(DataCalib) # LES NOMS DES COLONNES DE DATAVALID CORRESPONDENT A CEUX DE DATACALIB #
  
  setwd(pathResults) # DEFINITION DE EMPLACEMENT SORTIE #
  tiff(filename = paste(pathResults, "/", "MLR", "_", num, "_", Property, "_", DateImage, "_Distribution_Values.tiff", sep = ""), width = 2000, height = 2000, units = "px", bg = "white", res = 300) 
  
  selected_values=t(selected_values)
  
  # DIFFERENTE DONNEES DE HISTOGRAMME ET COULEUR # 
  hist(SpectreBerambadiCalib[1,], col = "red")
  hist(SpectreBerambadiValid[1,], add = TRUE, col = "blue")
  hist(selected_values[1,], add = TRUE, col = "green")
  
  # DONNEES EN ENTREE HISTOGRAMME # 
  NbTrain <- dim(SpectreBerambadiCalib)[2]
  NbTest <- dim(SpectreBerambadiValid)[2]
  NbT <- dim(selected_values)[2]
  
  legend("topright", c(paste(NbTrain, " Données de calibration"), paste(NbTest, " Données de validation"), paste(NbT, "Données de test")), fill = c("red", "blue", "green"))
  dev.off()
  
  # SAUVEGARDE DES POINTS DE CALIBRATION EN CSV # 
  SpectreBerambadiValid=t(SpectreBerambadiValid)
  write.csv(SpectreBerambadiValid, file = paste0(pathResults, "MLR_", num, "Points_Vaidation.csv"), row.names = FALSE)
  
  
  ###################################
  #   CONSTRUCTION DU MODELE MLR    #
  ###################################
  
  options(digits = 4) # DEFINITION DU NBR DE CHIFFRES SIGNIFICATIFS A AFICHER POUR LES NBR DANS LES SORTIES DE R # 
  
  cemonmodel <- lm(ProprieteSol ~ .,DataCalib) # CREATION DU MODELE OU LA VARIABLE A PREDIRE EST PROPRIETESOL ET VARIABLES DE DAACALIB COMME PREIDCATEURS # 
  print(cemonmodel)  # AFFICHAGE DES DETAILS DU MODELE DANS LA CONSOLE # 
  
  ProprieteSolpredite <- predict(cemonmodel, DataCalib) # UTILISATION SU MODELE CEMONMODEL POUR PREDIRE LES VALEURS DE PROPRIETESOL EN FONCTION DES VARIABLES DANS DATACALIB # 
  
  # CALCULE DU R² # 
  R2C2=cor(ProprieteSolpredite,DataCalib$ProprieteSol)*cor(ProprieteSolpredite,DataCalib$ProprieteSol) # CALCULE LE R² ENTRE VALEURS PREDITES ET VALEURS OBSERVEES # 
  R2C2  # AFFICHAGE DU R² # 
  R2C3<-1-(sum(((ProprieteSolpredite-DataCalib$ProprieteSol)^2)/sum((DataCalib$ProprieteSol-mean(DataCalib$ProprieteSol))^2))) # CALCULE R² AVEC UNE AUTRE FORMULATION # 
  R2C3 # AFFICHAGE DU R² # 
  
  SavePerfo[3]=R2C3
  
  # CALCULE DU RMSE # 
  Sum=0
  RMSEC2=sqrt(sum((ProprieteSolpredite-DataCalib$ProprieteSol)^2)/length(DataCalib$ProprieteSol)) # CALCULE LE RMSE ENTRE VALEURS PREDITES ET VALEURS OBSERVEES 
  RMSEC2 # AFFICHAGE DU RMSE # 
  
  SavePerfo[4]=RMSEC2 # ENREGISTREMENT DANS LA 4EME COLONNE DE SAVEPERFO # 
  propVrai=DataCalib$ProprieteSol 
  
  (biais=mean(propVrai-ProprieteSolpredite)) # CALCULE DU BIAIS COMME MOYENNE DES DIFF ENTRE VALEURS REELLES ET PREDITES # 
  
  (RPIQ<-(quantile(propVrai, 3/4)- quantile(propVrai, 1/4))/as.double(RMSEC2)) # CALCULE DU RPIQ AVEC DIFF ENTRE 3EME QUARTILE ET 1ER QUARTLE DES VALEURS REELLES ET PREDITES DIVISEES PAR RMSE # 
  
  (RPD<-(sd(propVrai))/ (as.double(RMSEC2))) # CALCULE DU RPD EN PRENANT ECART-TYPE DES VALEURS REELES / RMSE # 
  
  
  #############################
  #   VALIDATION DU MODELE    #
  #############################
  
  # PREDICTION DES VALEURS DE TENEUR EN MO POUR DONNEES DE VALIDATION # 
  MOPredictionVal<-predict(cemonmodel, newdata= DataValid)
  
  # CALCUL DU R² POUR DONNEES DE VALIDATION # 
  R2V=cor(MOPredictionVal,ProprieteSolValid)*cor(MOPredictionVal,ProprieteSolValid)
  SavePerfo[5]=R2V
  
  # CALCUL DU RMSE POUR DONNEES DE VALIDATION # 
  RMSEV=sqrt(sum((MOPredictionVal-ProprieteSolValid)^2)/length(ProprieteSolValid))
  SavePerfo[6]=RMSEV
  
  # CALCUL DU RPD POUR DONNEES DE VALIDATION # 
  RPDval<-(sd(ProprieteSolValid))/ (as.double(RMSEV))
  SavePerfo[7]=RPDval
  
  # CALCUL DU RPIQ POUR LES DONNEES DE VALIDATION #  
  RPIQval<-(quantile(ProprieteSolValid, 3/4)- quantile(ProprieteSolValid, 1/4))/as.double(RMSEV)
  SavePerfo[8]=RPIQval
  
  # CALCUL DU BIAIS POUR LES DONNEES DE VALIDATION # 
  biaisVal=mean(ProprieteSolValid)-mean(MOPredictionVal)
  SavePerfo[9]=biaisVal
  
  
  #######################
  #   TEST DU MODELE    #
  #######################
  
  # PREDICTION DES VALEURS DE TENEUR EN MO POUR DONNEES DE TEST # 
  MOPredictionTest<-predict(cemonmodel, newdata= DataTest)
  
  # CALCUL DU R² POUR DONNEES DE TEST # 
  R2T=cor(MOPredictionTest,ProprieteSolTest)*cor(MOPredictionTest,ProprieteSolTest)
  SavePerfo[10]=R2T
  
  # CALCUL DU RMSE POUR DONNEES DE TEST # 
  RMSET=sqrt(sum((MOPredictionTest-ProprieteSolTest)^2)/length(ProprieteSolTest))
  SavePerfo[11]=RMSET
  
  # CALCUL DU RPD POUR DONNEES DE TEST # 
  RPDT<-(sd(ProprieteSolTest))/ (as.double(RMSET))
  SavePerfo[12]=RPDT
  
  # CALCUL DU RPIQ POUR LES DONNEES DE TEST #  
  RPIQT<-(quantile(ProprieteSolTest, 3/4)- quantile(ProprieteSolTest, 1/4))/as.double(RMSET)
  SavePerfo[13]=RPIQT    
  
  # CALCUL DU BIAIS POUR LES DONNEES DE TEST # 
  biaisT=mean(ProprieteSolTest)-mean(MOPredictionTest)
  SavePerfo[14]=biaisT
  
  # AJOUT DES INDICES ERREUR DANS LA MATRICE EN FONCTION DE CHAQUE ITERATIONS # 
  metrics_list[[iteration]] <- c(num, R2C3, RMSEC2, biais, RPIQ, RPD, R2V, RMSEV, biaisVal, RPIQval, RPDval, R2T, RMSET, RPDT, RPIQT, biaisT)
  
  setwd(pathResults)
  tiff(filename = paste(pathResults, "/", "MLR", "_", num, "_", Property, "_", DateImage, "_Cal_Val_MLR.tiff", sep = ""), width = 6000, height = 2000, units = "px", bg = "white", res = 300)
  par(mfrow = c(1, 3))
  
  # Premier graphique avec légende
  plot(propVrai, ProprieteSolpredite, main = "Calibration Data", ylim = c(0, LimMaxscatterplot), xlim = c(0, LimMaxscatterplot), ylab = paste("Predicted Data (", prop_unite, ")"), xlab = paste("Observed Data (", prop_unite, ")"))
  abline(0, 1)
  RMSEC2=round(RMSEC2,digit=1)
  aaCV2<-paste("RMSE_Cal =",RMSEC2)
  R2C2=round(R2C2,digit=2)
  aaCV3<-paste("R2_Cal =",R2C2 )
  RPIQ=round(RPIQ,digit=1)
  aaCAL4<-paste("RPIQ_Cal =",RPIQ)
  RPD=round(RPD,digit=1)
  aaCAL5<-paste("RPD_Cal =",RPD)
  legend("bottomright", c(aaCV2, aaCV3, aaCAL4, aaCAL5), cex = 1, bty = "n")
  
  # Deuxième graphique avec légende
  plot(ProprieteSolValid, MOPredictionVal, main = "Validation Data", ylim = c(0, LimMaxscatterplot), xlim = c(0, LimMaxscatterplot), ylab = paste("Predicted Data (", prop_unite, ")"), xlab = paste("Observed Data (", prop_unite, ")"))
  abline(0, 1)
  biaisVal=round(biaisVal,digit=1)
  aaVAL1<-paste("bias =",biaisVal)
  RMSEV=round(RMSEV,digit=1)
  aaVAL2<-paste("RMSE =",RMSEV)
  R2CV=round(R2V,digit=2)
  aaVAL3<-paste("R2 =",R2CV)
  RPIQval=round(RPIQval,digit=1)
  aaVAL4<-paste("RPIQ =",RPIQval)
  RPDval=round(RPDval,digit=1)
  aaVAL5<-paste("RPD =",RPDval)
  legend("bottomright", c(aaVAL1, aaVAL2, aaVAL3, aaVAL4, aaVAL5), cex = 1, bty = "n")
  
  # Troisième graphique avec légende
  plot(ProprieteSolTest, MOPredictionTest, main = "Test Data", ylim = c(0, LimMaxscatterplot), xlim = c(0, LimMaxscatterplot), ylab = paste("Predicted Data (", prop_unite, ")"), xlab = paste("Observed Data (", prop_unite, ")"))
  abline(0, 1)
  biaisT=round(biaisT,digit=1)
  aaT1<-paste("bias =",biaisT)
  RMSET=round(RMSET,digit=1)
  aaT2<-paste("RMSE =",RMSET)
  R2T=round(R2T,digit=2)
  aaT3<-paste("R2 =",R2T)
  RPIQT=round(RPIQT,digit=1)
  aaT4<-paste("RPIQ =",RPIQT)
  RPDT=round(RPDT,digit=1)
  aaT5<-paste("RPD =",RPDT)
  
  legend("bottomright", c(aaT1, aaT2, aaT3, aaT4, aaT5), cex = 1, bty = "n")
  
  dev.off()
  
  #############################
  #   UNCERTAINTY METRICES    #
  #############################
  
  # PICP DE VALIDATION # 
  MOPrediVal <- predict(cemonmodel, newdata = DataValid)
  SaveMOPredictionVal[,num]=MOPrediVal # SAUVEGARDE DANS LA MATRICE #
  
  # PICP DE TEST #
  MOPrediTest <- predict(cemonmodel, newdata = DataTest)
  SaveMOPredictionTest[,num]=MOPrediTest # SAUVEGARDE DANS LA MATRICE #
  
  
  ##################################################################
  #   APPLICATION DU MODELE MLR A TOTALITE DES PIXELS DE SOL NU    #
  ##################################################################
  
  Mapped_iteration <- predict(S2Band.bare, cemonmodel, fun = predict) # CALCUL DES VALEURS PREDITES MO PAR LE MODELE # 
  
  Mapped_iteration[Mapped_iteration > 100] <- NA # REND EN NA LES VALEURS DEPASSANT 100% MO #
  Mapped_iteration[Mapped_iteration < 0] <- NA # REND NA LES VALEURS INFERIEURES A 0% MO #
  
  FinalMap=stack(FinalMap, Mapped_iteration) # SAUVEGARDE DE LA CARTE PREDITE # 
  
  non_selected_values=t(non_selected_values) # INVERSEMENT DES LIGNES ET COLONNES #
  selected_values=t(selected_values)
  SpectreBerambadiValid=t(SpectreBerambadiValid)
  
  iteration <- (iteration - 1)  #  FIN DE LA BOUCLE #
} # SORTIE DE LA BOUCLE LORSQUE LA VARIABLE "ITERATION" ATTEINT 0 # 

quantiles_5 <- calc(FinalMap, fun = function(x) quantile(x, probs = 0.05, na.rm = TRUE))
quantiles_95 <- calc(FinalMap, fun = function(x) quantile(x, probs = 0.95, na.rm = TRUE))

# CALCUL DE L'INTERVALLE DE PREDICTION
pred_interval_raster_MLR <- quantiles_95 - quantiles_5

# Sauvegarder les résultats en tant que nouvelles rasters
writeRaster(quantiles_5, filename = "quantiles_5.tif", format = "GTiff", overwrite = TRUE)
writeRaster(quantiles_95, filename = "quantiles_95.tif", format = "GTiff", overwrite = TRUE)

# Sauvegarder l'intervalle de prédiction
writeRaster(pred_interval_raster_MLR, filename = "intervalle_prediction.tif", format = "GTiff", overwrite = TRUE)


##########################
#   DONNEES EN SORTIE    #
##########################

# STATISTIQUES DES ERREURS # 
metrics_table <- do.call(rbind, metrics_list) # CREATION DU TABLEAU CONTENANT LES ERREURS # 
colnames(metrics_table) <- c("Itération", "R2_Cal", "RMSE_Cal", "Biais_Cal", "RPIQ_Cal", "RPD_Cal", "R2_Val", "RMSE_Val", "Biais_Val", "RPIQ_Val", "RPD_Val", "R2_Test", "RMSE_Test", "RPD_Test", "RPIQ_Test", "Biais_Test")

write.csv(metrics_table, file = paste0(pathResults, "MLR_tab_erreurs.csv"), row.names = FALSE) # SAUVEGARDE DU FICHIER # 

min_values <- apply(metrics_table[, -1], 2, min, na.rm = TRUE) # CALCUL DES VALEURS MINIMALES PAR COLONNE # 
max_values <- apply(metrics_table[, -1], 2, max, na.rm = TRUE) # CALCUL DES VALEURS MAXIMALES PAR COLONNE # 
mean_values <- apply(metrics_table[, -1], 2, mean, na.rm = TRUE) # CALCUL DES VALEURS MOYENNES PAR COLONNE #
median_values <- apply(metrics_table[, -1], 2, median, na.rm = TRUE) # CALCUL DES VALEURS MEDIANNES PAR COLONNE # 

# CREATION TABLEAU POUR LES STATISTIQUES # 
stats_table <- data.frame(Min = min_values, Max = max_values, Moyenne = mean_values, Médiane = median_values)
write.csv(stats_table, file = paste0(pathResults, "statistiques_erreurs.csv"), row.names = TRUE) # SAUVEGARDE DU FICHIER # 


##############
#   PICP    #
##############

# PICP DE VALIDATION #
lower_boundVal<- apply(SaveMOPredictionVal, 1, quantile, probs = 0.05) # le "1" dans la fonction sert à calculer le quantile pour CHAQUE echantillon, donc sur les lignes.
upper_boundVal <- apply(SaveMOPredictionVal, 1, quantile, probs = 0.95) # le "1" dans la fonction sert à calculer le quantile pour CHAQUE echantillon, donc sur les lignes.

# PICP DE TEST #
lower_boundTest<- apply(SaveMOPredictionTest, 1, quantile, probs = 0.05) # le "1" dans la fonction sert à calculer le quantile pour CHAQUE echantillon, donc sur les lignes.
upper_boundTest <- apply(SaveMOPredictionTest, 1, quantile, probs = 0.95) # le "1" dans la fonction sert à calculer le quantile pour CHAQUE echantillon, donc sur les lignes.

# Fonction pour calculer la PICP
calculate_picp <- function(observed, lower_bound, upper_bound) {
  within_interval <- (observed >= lower_bound) & (observed <= upper_bound)
  picp <- mean(within_interval)*100 #multiplier par 100 pour l'exprimer en %
  return(picp)
}

picp_val <- calculate_picp(ProprieteSolValid, lower_boundVal, upper_boundVal)   
picp_Test <- calculate_picp(ProprieteSolTest, lower_boundTest, upper_boundTest)   

# CHARGEMENT DES DONNEES DEPUIS LE CSV # 
data <- read.csv("MLR_tab_erreurs.csv")

# SELECTION DES VALEURS de biais DEPUIS LES COLONNES DU CSV # 
Biais_Cal <- data$Biais_Cal
Biais_Val <- data$Biais_Val
Biais_Test <- data$Biais_Test

tiff(paste0(pathResults, "MLR_boxplot_Biais.tiff"), width=800, height=600)
boxplot(list(Biais_Cal, Biais_Val, Biais_Test),
        main="Distribution du biais",
        ylab="Valeurs du biais",
        col=c(alpha("lightblue", 0.5), alpha("lightgreen", 0.5), alpha("lightpink", 0.5)),
        border=c("blue", "green", "red"),
        names=c("Biais Cal", "Biais Val", "Biais Test"))

# Ajout de la valeur de nb_total_samples
text(1, min(c(Biais_Cal, Biais_Val, Biais_Test)), paste("Nombre total d'échantillons :", nb_total_samples), pos=1, cex=0.8)
dev.off()

# SELECTION DES VALEURS R² DEPUIS LES COLONNES DU CSV # 
R2_Cal <- data$R2_Cal
R2_Val <- data$R2_Val
R2_Test <- data$R2_Test

tiff(paste0(pathResults, "MLR_boxplot_R2.tiff"), width=800, height=600)
boxplot(list(R2_Cal, R2_Val, R2_Test),
        main="Distribution de R²",
        ylab="Valeurs de R²",
        col=c(alpha("lightblue", 0.5), alpha("lightgreen", 0.5), alpha("lightpink", 0.5)),
        border=c("blue", "green", "red"),
        names=c("R² Cal", "R² Val", "R² Test"))

# Ajout de la valeur de nb_total_samples
text(1, min(c(R2_Cal, R2_Val, R2_Test)), paste("Nombre total d'échantillons :", nb_total_samples), pos=1, cex=0.8)
dev.off()

# SELECTION DES VALEURS RMSE DEPUIS LES COLONNES DU CSV # 
RMSE_Cal <- data$RMSE_Cal
RMSE_Val <- data$RMSE_Val
RMSE_Test <- data$RMSE_Test

# CREATION ET SAUVEGARDE DU BOXPLOT AVEC RMSE # 
tiff(paste0(pathResults, "MLR_boxplot_RMSE.tiff"), width=800, height=600)  
boxplot(list(RMSE_Cal, RMSE_Val, RMSE_Test),
        main="Distribution de RMSE",
        ylab="Valeurs de RMSE",
        col=c(alpha("lightblue", 0.5), alpha("lightgreen", 0.5), alpha("lightpink", 0.5)),
        border=c("blue", "green", "red"),
        names=c("RMSE Cal", "RMSE Val", "RMSE Test"))

# Ajout de la valeur de nb_total_samples
text(1, min(c(RMSE_Cal, RMSE_Val, RMSE_Test)), paste("Nombre total d'échantillons :", nb_total_samples), pos=1, cex=0.8)
dev.off() 


# SELECTION DES VALEURS RPIQ DEPUIS LES COLONNES DU CSV # 
RPIQ_Cal <- data$RPIQ_Cal
RPIQ_Val <- data$RPIQ_Val
RPIQ_Test <- data$RPIQ_Test

# CREATION ET SAUVEGARDE DU BOXPLOT AVEC RPIQ # 
tiff(paste0(pathResults, "MLR_boxplot_RPIQ.tiff"), width=800, height=600)  
boxplot(list(RPIQ_Cal, RPIQ_Val, RPIQ_Test),
        main="Distribution de RPIQ",
        ylab="Valeurs de RPIQ",
        col=c(alpha("lightblue", 0.5), alpha("lightgreen", 0.5), alpha("lightpink", 0.5)),
        border=c("blue", "green", "red"),
        names=c("RPIQ Cal", "RPIQ Val", "RPIQ Test"))

# Ajout de la valeur de nb_total_samples
text(1, min(c(RPIQ_Cal, RPIQ_Val, RPIQ_Test)), paste("Nombre total d'échantillons :", nb_total_samples), pos=1, cex=0.8)
dev.off()

# SELECTION DES VALEURS RPD DEPUIS LES COLONNES DU CSV # 
RPD_Cal <- data$RPD_Cal
RPD_Val <- data$RPD_Val
RPD_Test <- data$RPD_Test

# CREATION ET SAUVEGARDE DU BOXPLOT AVEC RPD # 
tiff(paste0(pathResults, "MLR_boxplot_RPD.tiff"), width=800, height=600)  
boxplot(list(RPD_Cal, RPD_Val, RPD_Test),
        main="Distribution de RPD",
        ylab="Valeurs de RPD",
        col=c(alpha("lightblue", 0.5), alpha("lightgreen", 0.5), alpha("lightpink", 0.5)),
        border=c("blue", "green", "red"),
        names=c("RPD Cal", "RPD Val", "RPD Test"))

# Ajout de la valeur de nb_total_samples
text(1, min(c(RPD_Cal, RPD_Val, RPD_Test)), paste("Nombre total d'échantillons :", nb_total_samples), pos=1, cex=0.8)
dev.off()

# SUAVEGARDE DU STACK DES CARTES DE PREDICTION #
writeRaster(FinalMap, filename = paste0(pathResults, "MLR_Stack_Cartes.tif"), format = "GTiff")

#################################
#  EVALUATION DE INCERTITUDE   #
#################################

# EVALUATION DE INCERTITUDE # 
mean_MO <- calc(FinalMap, mean) # CREATION DE LA CARTE DE LA MOYENNE MO # 
ecart_type_MO <- calc(FinalMap,  sd) # CREATION DE LA CARTE ECART-TYPE  # 
variance_MO <- calc(FinalMap, var) # CREATION DE LA CARTE DE LA VARIANCE # 

# CALCUL DU COEFFICIENT DE VARIATION DES VALEURS PREDITES #
coefficient_variation <- function(FinalMap) {
  cv <- sd(FinalMap) / mean(FinalMap) * 100
  return(cv)
}
cv_raster <- calc(FinalMap, coefficient_variation) # CREATION DE LA CARTE DU COEFFICIENT DE VARIATION # 

#########################################
#  AFFICHAGE ET GENERATION DES CARTES   #
#########################################

levelplot(mean_MO, main = "Carte des valeurs moyennes de matière organique") # AFFICHAGE DE LA CARTE DE LA MOYENNE MO # 
levelplot(ecart_type_MO, main = "Carte de l'incertitude") # AFFICHAGE DE LA CARTE ECART-TYPE # 
levelplot(cv_raster, main = "Carte du coefficient de variation") # AFFICHAGE DE LA CARTE DU COEFFICIENT DE VARIATION # 
levelplot(variance_MO, main = "Variance de la matière organique sur 20 itérations") # AFFICHAGE DE LA CARTE DE LA VARIANCE # 

writeRaster(mean_MO, filename = paste0(pathResults, "MLR_","mean_MO.tif"), format = "GTiff") # SAUVEGARDE DE LA CARTE MOYENNE MO # 
writeRaster(ecart_type_MO, filename = paste0(pathResults, "MLR_","ecart_type_MO.tif"), format = "GTiff") # SAUVEGARDE DE LA CARTE ECART-TYPE # 
writeRaster(cv_raster, filename = paste0(pathResults, "MLR_","coeff_variation_MO.tif"), format = "GTiff") # SAUVEGARDE DE LA CARTE COEFFICIENT DE VARIATION # 
writeRaster(variance_MO, filename = paste0(pathResults, "MLR_", "variance_MO.tif"), format = "GTiff") # SAUVEGARDE DE LA CARTE DE LA VARIANCE # 


#######################################
##### MATRICE DE CONFUSION - QRF ######
#######################################

sample_points<-read.table(paste(path,"164_3Prop_BM_AOBO.csv",sep=""), row.names=1,head=T) # LECTURE DU .CSV #
sample_points <- sample_points[, -c(4, 5)]

# CHARGEMENT DE LA CARTE DE MOYENNE # 
mean_pred_raster_QRF <- raster("//147.99.18.48/shared-data/Maryline_Faviere/RES_STAGE/BERAMBADI/BOOSTRAP/MLR/MLR_mean_argile.tif")
sample_points$observation <- as.numeric(as.character(sample_points$Clay)) # EN FORMAT NUMERIQUE #

# EXTRACTIOND DES COORDONNEES #
coordinates(sample_points) <- sample_points[, 1:2]

# EXTRACTION DES VALEURS PREDICTIVES AUX POINTS ECHANTILLON #
pred_values <- extract(mean_pred_raster_QRF, sample_points)

# AJOUT DES PREDICTIONS AUX POINTS ECHANTILLON #                                                                                                                                        
sample_points$pred <- pred_values

# CLASSER OBSERVATION ET PREDICTIONS EN CLASSES #
classify_values <- function(values, breaks) {
  classified_values <- cut(values, breaks = breaks, labels = FALSE)
  return(classified_values)
}

# DEFINITION DES SEUILS POUR LES CLASSES #
breaks <- c(-Inf, 10, 20, 30, 40, Inf)

# CLASSER LES OBSERVATIONS ET PREDICTIONS #
sample_points$obs_class <- classify_values(sample_points$observation, breaks)
sample_points$pred_class <- classify_values(sample_points$pred, breaks)

# CREATION DU DATA FRAME AVEC CLASSES #
comparison_df <- as.data.frame(sample_points)
comparison_df <- comparison_df[!is.na(comparison_df$obs_class) & !is.na(comparison_df$pred_class), ]

# GENERER LA MATRICE DE CONFUSION #
confusion_matrix <- table(comparison_df$obs_class, comparison_df$pred_class)

rownames(confusion_matrix) <- c("Très Bas", "Bas", "Moyen", "Haut", "Très Haut")
colnames(confusion_matrix) <- c("Très Bas", "Bas", "Moyen", "Haut", "Très Haut")

print(confusion_matrix)

# CALCUL DES METRIQUES DE PERFORMANCES #
confusion_matrix_caret <- confusionMatrix(as.factor(comparison_df$pred_class), as.factor(comparison_df$obs_class))
print(confusion_matrix_caret)

# ENREGISTREMENT HEURE SYSTEME FIN EXECUTION # 
temps_fin <- Sys.time()

# CALCUL DE LA DIFFERENCE HEURE SYSTEME #
temps_ecoule <- temps_fin - temps_debut

# AFFICHAGE DU TEMPS  # 
temps_ecoule_en_minutes <- as.numeric(temps_ecoule, units = "mins")
print(paste("Temps écoulé:", round(temps_ecoule_en_minutes, 2), "minutes"))


################ END  ################ 
################ END  ################ 
################ END  ################ 
