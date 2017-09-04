##### This code produces the parameter estimates of the best model mentioned in Table 7. The code requires installation of the software spAir. 

#########
library(spAir)

###### reading the data files #######################################
spcheck1<-read.table("model/spTimer/spAir_no2_best_model/data/AURN_data_07_11.txt",header=T)

position<-which(is.na(spcheck1$aqum_no2)==TRUE);
spcheck1$aqum_no2[position]=mean(spcheck1$aqum_no2,na.rm=T);
####
spcheck = data.frame(index=spcheck1$index, lon=spcheck1$lon, lat=spcheck1$lat, year=spcheck1$year, month=spcheck1$month, day=spcheck1$day, obs_no2=spcheck1$obs_no2);
spcheck$sqrtaqm = sqrt(spcheck1$aqum_no2);
spcheck$type=spcheck1$type;
spcheck = spcheck[order(spcheck$index, spcheck$year, spcheck$month, spcheck$day),];

###### Step for preparing the data files #########

index_data=unique(spcheck$index);
len_index=length(index_data);
coords_all = as.matrix(unique(cbind(spcheck[,2:3])));
######
index_data=unique(spcheck$index);
len_index=length(index_data);
######

abb=c("Rural","Urban","RKS"); 
name1=abb;

total_length=length(spcheck[,1]); 
addition=matrix(0,nrow=total_length,ncol=(2*length(abb)));

########
name=character(length=(2*length(name1)));
for(i in 1:(length(name1)))
{
  name[(2*i)-1] = paste("fac",i,sep="");
  name[2*i]=name1[i];
}
for(i in 1:len_index)
{
  pos1=c(which(spcheck[,1]==index_data[i]));
  type=spcheck$type[pos1[1]];  
  pos=which(abb==type);
  
    for(j in 1:length(pos1))
    {
      addition[pos1[j],((2*pos)-1)]<-1;
      addition[pos1[j],(2*pos)]<-spcheck$sqrtaqm[pos1[j]];
    }
  
}
colnames(addition)<-name;
deduc=c(1,2);
addition=addition[,-deduc];
spcheck=cbind(spcheck,addition);
spcheck=spcheck[,-9];

####
# Define the coordinates
co_ind<-as.matrix(unique(cbind(spcheck[,1:3])));

######## choosing the validation sites ########
set.seed(11)
SitePred<-c(4,9,28,30,32,44,50,75,102,108,111,122,130,175,194);

#sites for fitting
gh<-co_ind[,1];
sitefit<-gh[!gh %in% SitePred];

new <- sample(sitefit, size=(54-length(SitePred)), replace = FALSE, prob = NULL); 
SitePred <- append(SitePred, new); 

sitefit<-gh[!gh %in% SitePred];

####  prepare the fitting and validation data #########
prediction_file <- merge(data.frame(index=SitePred), spcheck, by="index");
prediction_file<-prediction_file[order(prediction_file$index,prediction_file$year,prediction_file$month,prediction_file$day),];
##########
fitting_file<-merge(data.frame(index=sitefit), spcheck, by="index");
fitting_file<-fitting_file[order(fitting_file$index,fitting_file$year,fitting_file$month,fitting_file$day),];
##############
grid=read.table("data/aurn/grid_pop_eng.txt",header=F);
cand_coords=as.matrix(grid[,1:2]); 
cand_len = length(cand_coords[,1]);
############################
cand_coords1=rep(0,(2*cand_len));
for(i in 1:cand_len)
{
  cand_coords1[(2*(i-1))+1] = cand_coords[i,1];  cand_coords1[(2*(i-1))+2] = cand_coords[i,2];
} 
grid_prob=grid[,3];
#grid_prob=rep(1, length(grid[,1]));

# MCMC via Gibbs not using default choices
# hyper-parameters for the prior distributions
priors<-spT.priors(model="GPP",inv.var.prior=Gamm(2,1),beta.prior=Norm(0,10^10));
# initial values for the model parameters
initials<-spT.initials(model="GPP", sig2eps=0.01,sig2eta=0.5, beta=NULL, phi=0.005, sig2l=0.01);

# input for spatial decay

spatial.decay<-spT.decay(distribution="FIXED",value=0.001);

####
set.seed(100)
coords<-as.matrix(unique(cbind(fitting_file[,2:3])));
new.coords<-as.matrix(unique(cbind(prediction_file[,2:3])));
###########
# Define knots
knots.coords<-spT.grid.coords(Longitude=c(max(coords_all[,1]),min(coords_all[,1])),Latitude=c(max(coords_all[,2]),min(coords_all[,2])), by=c(5,5));

  posT <- c(0, 0)
  XY <- Formula.matrix(formula=obs_no2 ~ sqrtaqm + fac2 + Urban + fac3 + RKS, data=fitting_file);
  Y <- XY[[1]]
  X <- as.matrix(XY[[2]])
  time.data<-list(1,length(Y)/length(coords[,1]));
  
  data_time = unique(fitting_file[,c(4:6)]);
  data_time$loc = 1:time.data[[2]];
  pred_time = unique(prediction_file[,c(4:6)]);
  pos = which(data_time$year==pred_time$year[1] & data_time$month==pred_time$month[1] & data_time$day==pred_time$day[1]);
  posT[1] = pos - 1

  end <- dim(pred_time)[1]
  pos <- which(data_time$year==pred_time$year[end] & data_time$month==pred_time$month[end] & data_time$day==pred_time$day[end])
  posT[2] = pos -1

 

post.ar.fitpred <- spT.Gibbs(formula=obs_no2 ~sqrtaqm + fac2 + Urban + fac3 + RKS, data=fitting_file, model="GPP", coords=coords, knots.coords=knots.coords, newcoords=new.coords, newdata=prediction_file,nItr=10000,nBurn=5000,tol.dist=0.0001,distance.method="geodetic:km",initials=initials,priors=priors,scale.transform="SQRT",spatial.decay=spatial.decay
,cand_coords=cand_coords1,cand_len=cand_len,grid_prob=grid_prob,report=1000, predloc=posT);
  
print(post.ar.fitpred)

##### produces the columns of table 7 ####################################
summary(model.output)
names(model.output)

parameter = model.output$parameters[,c(1,4,5)];
rownames(parameter) = c("gamma_0", "gamma_1", "gamma_02", "gamma_12", "gamma_03", "gamma_13", "rho", "sigma_epsilon", "sigma_eta", "phi");

#### produces the RMSPE and MAPE #########
rmse=spT.validation(prediction_file$obs_no2,c(model.output$prediction[,1]))

bias = mean(prediction_file$obs_no2 - c(model.output$prediction[,1]), na.rm=T);


#####  Producing the row of the Table 3 #############

limits = read.table("Prediction_sites_summary.txt");
limits = data.frame(mean=limits[,1], sd=limits[,2], low95=limits[,3], up95=limits[,4]);
limits = cbind(prediction_file[,c(1:7)], limits);
limits = limits[which(is.na(limits$obs_no2)==F),];
total_length = length(limits[,1]);
prop = 0;
for(ii in 1:total_length){
  if((limits$low95[ii] < limits$obs_no2[ii]) & (limits$obs_no2[ii] < limits$up95[ii])){
      prop = prop + 1;
  }
}
prop = 100*prop/total_length;

####### printing result will give row of Table 3 ############
result = matrix(0, nrow=1, ncol=5);

result[1,1] = rmse[2];  result[1,2] = rmse[3]; result[1,3] = bias; result[1,4]=prop; result[1,5] = (cor(c(model.output$prediction[,1]), prediction_file$obs_no2, use="pairwise.complete.obs"))^2;
rownames(result) = c("Linear");
colnames(result) = c("RMSPE", "MAPE", "Bias", "Coverage (%)", "R2");



## 
# nBurn =  1000 . Iterations =  5000 . 
# Acceptance rate: (phi) =   % 
## 
# Text output is given: 
##
##
# Elapsed time: 2 Hour/s. 3 Min. 0.86 Sec.
##
# >   
# > print(post.ar.fitpred)
# -----------------------------------------------------
# Model: GPP
# Call: obs_no2 ~ sqrtaqm + fac2 + Urban + fac3 + RKS
# Iterations: 5000
# nBurn: 1000
# Acceptance rate for phi (%): 
# -----------------------------------------------------
#         Goodness.of.fit  Penalty     PMCC
# values:        221330.6 177381.3 398711.9
# -----------------------------------------------------
# Computation time: 2  - Hour/s. 3  - Mins. 0.86  - Sec.
# > 
# > ##### produces the columns of table
# > summary(post.ar.fitpred)
# -----------------------------------------------------
# Model: GPP
# Call: obs_no2 ~ sqrtaqm + fac2 + Urban + fac3 + RKS
# Iterations: 5000
# nBurn: 1000
# Acceptance rate for phi (%): 
# -----------------------------------------------------
#         Goodness.of.fit  Penalty     PMCC
# values:        221330.6 177381.3 398711.9
# -----------------------------------------------------
# Computation time: 2  - Hour/s. 3  - Mins. 0.86  - Sec.
# -----------------------------------------------------
# Parameters:
#               Mean Median     SD Low2.5p Up97.5p
# (Intercept) 4.0968 4.0944 0.0643  3.9792  4.2332
# sqrtaqm     0.1868 0.1869 0.0055  0.1758  0.1975
# fac2        1.3738 1.3736 0.0335  1.3078  1.4401
# Urban       0.0504 0.0504 0.0059  0.0389  0.0617
# fac3        1.7879 1.7876 0.0356  1.7181  1.8607
# RKS         0.1237 0.1237 0.0062  0.1114  0.1359
# rho         0.1011 0.1011 0.0066  0.0882  0.1139
# sig2eps     1.7001 1.7001 0.0062  1.6882  1.7127
# sig2eta     0.9160 0.9152 0.0220  0.8762  0.9627
# phi         0.0010 0.0010 0.0000  0.0010  0.0010
# -----------------------------------------------------
# > names(post.ar.fitpred)
#  [1] ""                  "accept"            "call"             
#  [4] "parameters"        "fitted"            "prediction"       
#  [7] "tol.dist"          "distance.method"   "cov.fnc"          
# [10] "scale.transform"   "sampling.sp.decay" "PMCC"             
# [13] "n"                 "knots"             "pred.n"           
# [16] "r"                 "T"                 "initials"         
# [19] "priors"            "iterations"        "nBurn"            
# [22] "computation.time"  "combined.fit.pred" "model"            
# > 
# >

