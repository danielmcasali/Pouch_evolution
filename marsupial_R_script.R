############################################################################################################################################################
#script: marsupium evolution x body mass x litter size
#author: Daniel Casali
############################################################################################################################################################

#Load packages
packages<-c("evobiR","phytools","dplyr","tidyr","ggplot2","ggpubr","phyr","OUwie","car","plotrix","svglite","patchwork","phylosignalDB","dentist")
sapply(packages, function(x) library(x, character.only = TRUE))

#Load data
dat<-read.csv("marsupial_data.csv",row.names=1)

#Load major clade credibility tree (from Upham et al. 2019), prunned to taxomic sample in data
mcc<-read.tree("marsupial_MCC.nwk")
mcc

#Load 100 posterior sample trees (from Upham et al. 2019), prunned to taxomic sample in data
#tree-pruner-79f70b71-c487-403a-94a8-2d5bbe196bb2
pst<-read.tree("marsupial_PST.trees")
pst

##Prepare object "dat_mean"
{
#dat_mean<-dat[c("Pouch","BM_mean","LS_mean")]
#colnames(dat_mean)<-c("pouch","BM","LS")
#dat_mean$pouch[dat_mean$pouch=="absent"]<-0
#dat_mean$pouch[dat_mean$pouch=="present"]<-1
#dat_mean$pouch<-as.numeric(dat_mean$pouch)
#dat_mean$BM<-round(log10(dat_mean$BM),2)
#dat_mean$zmass<-(dat_mean$BM-mean(dat_mean$BM))/sd(dat_mean$BM)
#dat_mean$zlitter<-(dat_mean$LS-mean(dat_mean$LS))/sd(dat_mean$LS)
#dat_mean
#saveRDS(dat_mean,"dat_mean.rds")
}
dat_mean<-readRDS("dat_mean.rds")

##Prepare objects "dat_100" (sampling values from min-max ranges 100 times)
{
#dat_100<-list()
#vec_prov<-rep(0,length(dat_mean$pouch))
#dat_prov<-data.frame(pouch=dat_mean$pouch,BM=vec_prov,LS=vec_prov)
#rownames(dat_prov)<-rownames(dat)
#colnames(dat_prov)<-c("pouch","BM","LS")
#for (i in 1:100)
#{
#	for (j in 1:length(dat$Pouch))
#	{
#		if (is.na(dat$BM_min[j])) {dat_prov$BM[j] <- dat$BM_mean[j]
#		} else {dat_prov$BM[j] <- round(runif(n=1, min=dat$BM_min[j], max=dat$BM_max[j]),2)
#		} # if a range is available, sample from it, otherwise, get the mean.
#		
#		if (is.na(dat$LS_min[j])) {dat_prov$LS[j] <- dat$LS_mean[j]
#		} else {dat_prov$LS[j] <- round(runif(n=1, min=dat$LS_min[j], max=dat$LS_max[j]),2)
#		} # if a range is available, sample from it, otherwise, get the mean.
#	}
#	#log-tranforming, getting Z-scores and saving each dataset
#	dat_prov$BM<-round(log10(dat_prov$BM),2)
#	dat_prov$zmass<-(dat_prov$BM-mean(dat_prov$BM))/sd(dat_prov$BM)
#	dat_prov$zlitter<-(dat_prov$LS-mean(dat_prov$LS))/sd(dat_prov$LS)		
#	dat_100[[i]]<-dat_prov
#}
#saveRDS(dat_100,"dat_100.rds")
}
dat_100<-readRDS("dat_100.rds")

####################################################################Useful functions########################################################################
#Calculate mean and CI (simple and Rubin's)
mean_ci <- function(x, x.se=NULL, alpha = 0.05, dec=2) 
{		
	M <- length(x)  				   								                           #number of replicates
	Q_bar <- mean(x)                   								                           #among-datasets mean
	U_bar <- mean(x.se^2)              								                           #among-datasets variance  
	B_between <- var(x)                								                           #between-imputation variance (B_between)
	T <- U_bar + (1 + 1/M) * B_between 								                           #total variance with the finite-M correction
	SE_total <- sqrt(T)				   								                           #total standard error	
	df_total <- (M - 1) * (1 + (U_bar / ((1 + 1/M) * B_between)))^2                            #degrees of freedom (df) correction
	t_score <- qt(p = 1 - alpha/2, df = df_total)                                              #t-score for the desired alpha level    
	CI <- round(Q_bar + c(-1, 1) * t_score * SE_total,dec)			                           #95% Rubin's pooled confidence interval
	se_mean <- sd(x)/sqrt(M)																   #Std. error of the mean						   
	IN <- ci <- Q_bar + c(-1, 1) * qnorm(1 - alpha/2) * se_mean 							   #95% CI of the meam
	
	if (is.null(x.se)) { 
		result<-paste0(round(Q_bar,dec)," (",paste0(round(IN[1],dec)," – ",round(IN[2],dec)),")")
	} else { 
		result<-paste0(round(Q_bar,dec)," (",paste0(round(IN[1],dec)," – ",round(IN[2],dec)),")"," [", paste0(CI[1]," – ",CI[2]), "]")
	}
	return(result)
}

#BICw function
BICw <- function(BIC_values) 
{
	delta <- BIC_values - min(BIC_values)
	weights <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))
	return(weights)	
}
##########################################################Phylogenetic generalized liner mixed models (PGLMM)###############################################

#Calculate the variance inflation factor (VIF) for body mass and litter size across datasets
VIF<-as.data.frame(matrix(nrow=length(dat_100),ncol=3))
colnames(VIF)<-c("BM","LS","BM:LS")
for (i in 1:length(dat_100))
{
	dat2<-dat_100[[i]]
	reg<-glm(pouch~zmass*zlitter, data = dat2, family = binomial)
	VIF[i,]<- round(vif(reg),2)
}
VIF_summary<-t(sapply(VIF, mean_ci))
VIF_summary
write.csv(VIF_summary,"VIF_summary.csv",quote=FALSE,row.names=TRUE,fileEncoding = "UTF-8")

#Perform phylogenetic logistic regressions (pairing 100 trees with 100 simulated data to conduct regressions for each of the 5 models)
model_names<-c("Null","Body mass","Litter size","Body mass + Litter Size","Body mass * Litter Size")
model_vector<-vector("list",length(model_names))
names(model_vector)<-model_names
PGLM_results<-replicate(length(pst), model_vector, simplify = FALSE)

for (i in 1:length(pst)){
	
	print(paste("Tree and dataset",i))
	TR<-pst[[i]]
	dat3<-ReorderData(TR,dat_100[[i]])
	
	PGLM_results[[i]][[1]]<-tryCatch({pglmm_compare(pouch~1,data=dat3,family="binomial",phy=TR,s2.init=0.01)},error=function(e) e$message)
	PGLM_results[[i]][[2]]<-tryCatch({pglmm_compare(pouch~zmass,data=dat3,family="binomial",phy=TR,s2.init=0.01)},error=function(e) e$message)
	PGLM_results[[i]][[3]]<-tryCatch({pglmm_compare(pouch~zlitter,data=dat3,family="binomial",phy=TR,s2.init=0.01)},error=function(e) e$message)
	PGLM_results[[i]][[4]]<-tryCatch({pglmm_compare(pouch~zmass+zlitter,data=dat3,family="binomial",phy=TR,s2.init=0.01)},error=function(e) e$message)
	PGLM_results[[i]][[5]]<-tryCatch({pglmm_compare(pouch~zmass*zlitter,data=dat3,family="binomial",phy=TR,s2.init=0.01)},error=function(e) e$message)	
}
saveRDS(PGLM_results,"PGLM_results.rds")

#Summarize results
	
	#Extracting results
	model_names<-c("Null","Body mass","Litter size","Body mass + Litter Size","Body mass * Litter Size")
	PGLM_results<-readRDS("PGLM_results.rds")
	PGLM_df<-as.data.frame(matrix(nrow=5,ncol=13))	
	colnames(PGLM_df)<-c("Model","Intercept","Body mass","Litter size","Body mass:Litter Size",
						 "Intercept_se","Body mass_se","Litter size_se","Body mass:Litter Size_se",
						 "s2","BIC","BICw","Convergence")		
	
	PGLM_summary<-replicate(length(PGLM_results),PGLM_df,simplify=FALSE)	
	for (i in 1:length(PGLM_results))
	{
		for (j in 1:length(PGLM_results[[i]]))
		{               
			PGLM_summary[[i]][j,][1]  <- model_names[j]				
			if (class(PGLM_results[[i]][[j]])=="pglmm_compare") {
			PGLM_summary[[i]][j,][2]  <- as.vector(PGLM_results[[i]][[j]]$B[1])
			PGLM_summary[[i]][j,][3]  <- as.vector(PGLM_results[[i]][[j]]$B[2])
			PGLM_summary[[i]][j,][4]  <- as.vector(PGLM_results[[i]][[j]]$B[3])
			PGLM_summary[[i]][j,][5]  <- as.vector(PGLM_results[[i]][[j]]$B[4])
			PGLM_summary[[i]][j,][6]  <- as.vector(PGLM_results[[i]][[j]]$B.se[1])
			PGLM_summary[[i]][j,][7]  <- as.vector(PGLM_results[[i]][[j]]$B.se[2])
			PGLM_summary[[i]][j,][8]  <- as.vector(PGLM_results[[i]][[j]]$B.se[3])
			PGLM_summary[[i]][j,][9]  <- as.vector(PGLM_results[[i]][[j]]$B.se[4])
			PGLM_summary[[i]][j,][10] <- as.vector(PGLM_results[[i]][[j]]$ss)
			PGLM_summary[[i]][j,][11] <- as.vector(PGLM_results[[i]][[j]]$BIC)
			PGLM_summary[[i]][j,][13] <- as.vector(PGLM_results[[i]][[j]]$convcode)
			} else {
				PGLM_summary[[i]][j,][2:12]<-NA
			}
		}		
		PGLM_summary[[i]][12]<-round(BICw(PGLM_summary[[i]][11]),2)	
	}
	
	##Summarize and export
	PGLM_summary2<-do.call(rbind,PGLM_summary) #here check convergence, and if necessary, replace BICw of non-convergent to NA (not the case here)
	PGLM_summary3<-PGLM_summary2[!is.na(PGLM_summary2$BICw),]#remove results for tree:data comb. that at least one model presented issues
	PGLM_summary4<-split(PGLM_summary3, PGLM_summary3$Model)[model_names]
	
	PGLM_summary_df<-as.data.frame(matrix(nrow=1,ncol=8))
	colnames(PGLM_summary_df)<-c("Model","Intercept","Body mass","Litter size","Body mass:Litter Size","s2","BIC","BICw")
	
	PGLM_summary5<-replicate(length(PGLM_summary4), PGLM_summary_df, simplify=FALSE)
	for (i in 1:length(PGLM_summary5))
	{		
		PGLM_summary5[[i]][1] <- model_names[i]
		PGLM_summary5[[i]][2] <- mean_ci(PGLM_summary4[[i]][[2]],PGLM_summary4[[i]][[6]])#intercept
		PGLM_summary5[[i]][3] <- mean_ci(PGLM_summary4[[i]][[3]],PGLM_summary4[[i]][[7]])#bm
		PGLM_summary5[[i]][4] <- mean_ci(PGLM_summary4[[i]][[4]],PGLM_summary4[[i]][[8]])#ls
		PGLM_summary5[[i]][5] <- mean_ci(PGLM_summary4[[i]][[5]],PGLM_summary4[[i]][[9]])#bm:ls
		PGLM_summary5[[i]][6] <- mean_ci(PGLM_summary4[[i]][[10]])
		PGLM_summary5[[i]][7] <- mean_ci(PGLM_summary4[[i]][[11]])
		PGLM_summary5[[i]][8] <- mean_ci(PGLM_summary4[[i]][[12]])
	}	
	PGLM_summary_final<-do.call(rbind,PGLM_summary5)
	PGLM_summary_final[PGLM_summary_final == "NA (NA – NA) [NA – NA]"] <- NA
	PGLM_summary_final[3,4]<-PGLM_summary_final[3,3]
	PGLM_summary_final[3,3]<-NA
	write.csv(PGLM_summary_final,"Summary_PGLM.csv",row.names=FALSE,fileEncoding = "UTF-8")

#Plots

##Biplot
dat_plot<-dat
dat_plot$BM_mean<-log10(dat_plot$BM_mean)
dat_plot$BM_min<-log10(dat_plot$BM_min)
dat_plot$BM_max<-log10(dat_plot$BM_max)
biplot<-ggplot(dat=dat_plot,aes(x=BM_mean,y=LS_mean,color=Pouch,shape=Pouch)) +
geom_errorbarh(aes(xmin = BM_min, xmax = BM_max), height = 0, alpha=0.4) +
geom_errorbar(aes(ymin = LS_min, ymax = LS_max), width = 0, alpha=0.4) +
geom_point(alpha=1, size=2) +
theme_classic() + theme(legend.position=c(0.12,0.9),legend.title = element_text(face = "bold")) + 
labs (x="Body mass (log10)", y="Litter size", title="") +
scale_color_viridis_d()+
scale_shape_manual(values=c(19,19))
biplot

##PGLM plots
PGLM_results<-readRDS("PGLM_results.rds")
pglm_res<-list()
for (i in 1:length(PGLM_results))
{
	pglm_res[[i]]<- data.frame(species=PGLM_results[[i]][[4]]$phy$tip.label,
							   pouch=PGLM_results[[i]][[4]]$data$pouch,
							   bm=PGLM_results[[i]][[4]]$data$BM,
							   ls=PGLM_results[[i]][[4]]$data$LS,
							   mu=rowMeans(PGLM_results[[i]][[4]]$mu,PGLM_results[[i]][[5]]$mu))
}
pglm_res2<-do.call(rbind,pglm_res)
pglm_res3<-aggregate(. ~ species, data = pglm_res2, FUN = mean)
pglm_res3$Pouch<-pglm_res3$pouch
pglm_res3$Pouch[pglm_res3$Pouch==0]<-"absent"
pglm_res3$Pouch[pglm_res3$Pouch==1]<-"present"

BM_S<-ggplot(pglm_res3, aes(x = bm, y = pouch, color=Pouch)) + geom_point(alpha = 0.5, size=2, show.legend = FALSE) + scale_color_viridis_d()+
labs(x = "Body mass (log10)", y = "Probability of pouch") + theme_classic() +
geom_smooth(aes(x = bm, y = mu), method = "glm", method.args = list(family = "binomial"), color = "black", fill = "grey80", size=0.5)

LS_S<-ggplot(pglm_res3, aes(x = ls, y = pouch, color=Pouch)) + geom_point(alpha = 0.5, size=2, show.legend = FALSE) + scale_color_viridis_d()+
labs(x = "Litter size", y = "Probability of pouch") + theme_classic()+
geom_smooth(aes(x = ls, y = mu), method = "glm", method.args = list(family = "binomial"), color = "black", fill = "grey80", size=0.5)

#Get plots together and export
SS<-ggarrange(BM_S,LS_S, ncol = 1, nrow = 2, align = "h", labels = c("B","C"), hjust=-0.3)
FF<-ggarrange(biplot,NULL,SS, ncol = 3, nrow = 1, align = "v", labels = c("A",""), hjust=-0.3, widths=c(1.5,-0.2,1))
ggsave("Figure_1.svg",width=7.5,height=5)

#####################################################DISCRETE + CONTINUOUS MODELS######################################################################

#BODY MASS - Prepare dataset and set the number of stochastic maps (nSim) and random starts (n_st).
dataset<-data.frame(taxon=rownames(dat_mean),pouch=dat_mean$pouch,mass=dat_mean$BM)
nSim<-50
n_st=10
model_names<-c("CD_ER_BMV", "CD_ER_OUV", "CD_ER_OUM", "CD_ER_OUMV","CD_ARD_BMV", "CD_ARD_OUV", "CD_ARD_OUM", "CD_ARD_OUMV",
			   "CID_ER_BM1", "CID_ER_OU1", "CID_ARD_BM1", "CID_ARD_OU1",
			   "CIDP_ER_BMV", "CIDP_ER_OUV", "CIDP_ER_OUM", "CIDP_ER_OUMV","CIDP_ARD_BMV", "CIDP_ARD_OUV", "CIDP_ARD_OUM", "CIDP_ARD_OUMV",
			   "HYB_ER_BMV", "HYB_ER_OUV", "HYB_ER_OUM", "HYB_ER_OUMV","HYB_ARD_BMV", "HYB_ARD_OUV", "HYB_ARD_OUM", "HYB_ARD_OUMV")
model_results<-vector("list",length(model_names))
names(model_results)<-model_names   

	##Character-dependent models (CD)	
	model_results$CD_ER_BMV      <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "BMV",  null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$CD_ER_OUV      <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "OUV",  null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$CD_ER_OUM      <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "OUM",  null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$CD_ER_OUMV     <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "OUMV", null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$CD_ARD_BMV     <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "BMV",  null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$CD_ARD_OUV     <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "OUV",  null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$CD_ARD_OUM     <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "OUM",  null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$CD_ARD_OUMV    <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "OUMV", null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	
	##Character-independent models (CID)	
	model_results$CID_ER_BM1     <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "BM1", null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = TRUE, adaptive_sampling = TRUE)
	model_results$CID_ER_OU1     <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "OU1", null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = TRUE, adaptive_sampling = TRUE)
	model_results$CID_ARD_BM1    <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "BM1", null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = TRUE, adaptive_sampling = TRUE)
	model_results$CID_ARD_OU1    <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "OU1", null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = TRUE, adaptive_sampling = TRUE)
	
	##Character-independent plus models (CID+)
	model_results$CIDP_ER_BMV    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "BMV",  null.model=TRUE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = TRUE, adaptive_sampling = TRUE)
	model_results$CIDP_ER_OUV    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUV",  null.model=TRUE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = TRUE, adaptive_sampling = TRUE)
	model_results$CIDP_ER_OUM    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUM",  null.model=TRUE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = TRUE, adaptive_sampling = TRUE)
	model_results$CIDP_ER_OUMV   <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUMV", null.model=TRUE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = TRUE, adaptive_sampling = TRUE)
	model_results$CIDP_ARD_BMV   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "BMV",  null.model=TRUE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = TRUE, adaptive_sampling = TRUE)
	model_results$CIDP_ARD_OUV   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUV",  null.model=TRUE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = TRUE, adaptive_sampling = TRUE)
	model_results$CIDP_ARD_OUM   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUM",  null.model=TRUE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = TRUE, adaptive_sampling = TRUE)
	model_results$CIDP_ARD_OUMV  <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUMV", null.model=TRUE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = TRUE, adaptive_sampling = TRUE)
	
	##Hybrid models (HYB)
	model_results$HYB_ER_BMV     <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "BMV",  null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$HYB_ER_OUV     <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUV",  null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$HYB_ER_OUM     <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUM",  null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$HYB_ER_OUMV    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUMV", null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$HYB_ARD_BMV    <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "BMV",  null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$HYB_ARD_OUV    <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUV",  null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$HYB_ARD_OUM    <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUM",  null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)
	model_results$HYB_ARD_OUMV   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUMV", null.model=FALSE, nSim=nSim, n_starts=n_st, ncores=n_st, sample_nodes = FALSE, adaptive_sampling = FALSE)

	##Save model results
	saveRDS(model_results,"model_results.rds")

	##Sumarize BIC
	model_results<-readRDS("model_results.rds")
	models_BIC<-as.data.frame(matrix(nrow=length(model_results),ncol=2))
	colnames(models_BIC)<-c("Model","BIC")
	for (i in 1:length(model_results))
	{
		models_BIC[i,1]<-names(model_results)[i]
		models_BIC[i,2]<-round(model_results[[i]]$BIC,2)
	}
	models_BIC$BICw<-round(BICw(models_BIC$BIC),2)
	models_BIC
	write.csv(models_BIC,"models_BIC.csv")

	
#Ancestral state estimations	
	
	model_results<-readRDS("model_results.rds")
	
	mA<-hOUwie.recon(model_results$CD_ER_BMV, nodes = "internal")
	saveRDS(mA,"mA_ase.rds")	
	
	mB<-hOUwie.recon(model_results$CD_ER_OUMV, nodes = "internal")
	saveRDS(mB,"mB_ase.rds")
	
	#PLOTS	
	mA<-readRDS("mA_ase.rds")
	mB<-readRDS("mB_ase.rds")
	tip_data<-as.factor(setNames(dat$Pouch,rownames(dat)))
	cols<-make.transparent(c("#440154FF","#FDE725FF"),alpha=0.9)
	svg("ASR.svg",width=8,height=9.2)		
	
		par(mfrow=c(2,1))
		arc_height<-0.01
		states <- c("absent","present")
		
		plotTree(mcc, type="arc", fsize=0.3, ftype="i", offset = 3, lwd=0.3, arc_height=arc_height)
		title(main="", adj=0.35, line = -1.0, cex.main=1.5)
		h=max(nodeHeights(mcc))
		labs<-seq(0,h,by=10)
		a1<-axis(1,pos=-0.02*h,at=h-labs+arc_height*h,
		labels=labs,cex.axis=0.7,lwd=1,lend=1,padj=-1.5)
		text(mean(a1),-0.13*h,"million years ago",font=1, cex=0.8)
		a2<-axis(1,pos=-0.02*h,at=-h+labs-arc_height*h,
		labels=labs,cex.axis=0.7,lwd=1,lend=1,padj=-1.5)
		text(mean(a2),-0.13*h,"million years ago",font=1, cex=0.8)
		draw.arc(0,0,radius=h-labs[2:length(labs)]+arc_height*h,
		angle1=,angle2=pi,col=make.transparent("grey",1),lty=1,lwd=0.5)
		par(lty="solid",fg="transparent")
		nodelabels(node=1:mcc$Nnode+Ntip(mcc), pie=mA, piecol=cols, cex=0.20)
		tiplabels(pie=to.matrix(tip_data[mcc$tip.label],levels(tip_data)),piecol=cols,cex=0.15)
		par(lty="solid",fg="black")	
		legend(x=-103,y=98,states,cex=0.80, pt.cex=1.3, pch=21, pt.bg = cols, col = "white", bty="n",
		title="Pouch", title.adj = 0.25, title.cex = 1, title.font = 2)
		title(paste("A"),adj=0.03,line=-1.5,cex.main=1.2)
		
		plotTree(mcc, type="arc", fsize=0.3, ftype="i", offset = 3, lwd=0.3, arc_height=arc_height)
		title(main="", adj=0.35, line = -1.0, cex.main=1.5)
		h=max(nodeHeights(mcc))
		labs<-seq(0,h,by=10)
		a1<-axis(1,pos=-0.02*h,at=h-labs+arc_height*h,
		labels=labs,cex.axis=0.7,lwd=1,lend=1,padj=-1.5)
		text(mean(a1),-0.13*h,"million years ago",font=1, cex=0.8)
		a2<-axis(1,pos=-0.02*h,at=-h+labs-arc_height*h,
		labels=labs,cex.axis=0.7,lwd=1,lend=1,padj=-1.5)
		text(mean(a2),-0.13*h,"million years ago",font=1, cex=0.8)
		draw.arc(0,0,radius=h-labs[2:length(labs)]+arc_height*h,
		angle1=,angle2=pi,col=make.transparent("grey",1),lty=1,lwd=0.5)
		par(lty="solid",fg="transparent")
		nodelabels(node=1:mcc$Nnode+Ntip(mcc), pie=mB, piecol=cols, cex=0.20)
		tiplabels(pie=to.matrix(tip_data[mcc$tip.label],levels(tip_data)),piecol=cols,cex=0.15)
		par(lty="solid",fg="black")	
		legend(x=-103,y=98,states,cex=0.80, pt.cex=1.3, pch=21, pt.bg = cols, col = "white", bty="n",
		title="Pouch", title.adj = 0.25, title.cex = 1, title.font = 2)
		title(paste("B"),adj=0.03,line=-1.8,cex.main=1.2)		
		
	dev.off()
	
#Calculate 95% CI for hOUwie models with AICw >= 0.1
fn_hOUwie <- function(par, phy, data, rate.cat, discrete_model, continuous_model, null.model, nSim)
{
	hOUwie_fit <-hOUwie(phy=phy, data=data, rate.cat=rate.cat, discrete_model=discrete_model, continuous_model=continuous_model, 
	null.model=null.model, nSim=nSim, p=par)
	loglik <- hOUwie_fit$loglik
	neg_loglik <- -loglik
	return(neg_loglik)
}

 ##body mass
 model_results<-readRDS("model_results.rds")
 BEST<-model_results$CD_ER_BMV
 par<-BEST$p
 names(par)<-c("Rate 0<->1", "Sigma2_0", "Sigma2_1", "Theta")
 dent_A <- dent_walk(par=par, fn=fn_hOUwie, best_neglnL=-BEST$loglik, phy=BEST$phy, data=BEST$data, rate.cat=BEST$rate.cat,
 discrete_model=BEST$discrete_model, continuous_model=BEST$continuous_model, null.model=FALSE, nSim=BEST$nSim, nsteps = 2000)
 saveRDS(dent_A,"dent_A.rds")
 pdf("dent_A.pdf",height=10,width=10)
 plot(dent_A)
 dev.off()
 
 BEST<-model_results$CD_ER_OUMV
 par<-BEST$p
 names(par)<-c("Rate 0<->1", "Alpha", "Sigma2_0", "Sigma2_1", "Theta_0", "Theta_1")
 dent_B <- dent_walk(par=par, fn=fn_hOUwie, best_neglnL=-BEST$loglik, phy=BEST$phy, data=BEST$data, rate.cat=BEST$rate.cat,
 discrete_model=BEST$discrete_model, continuous_model=BEST$continuous_model, null.model=FALSE, nSim=BEST$nSim, nsteps = 2000)
 saveRDS(dent_B,"dent_B.rds")
 pdf("dent_B.pdf",height=10,width=12)
 plot(dent_B)
 dev.off()
 
 ##export param tables
	format_ci_table <- function(tab, digits = 2) 
	{
	  
	  stopifnot(nrow(tab) == 3)
	  formatted <- sapply(1:ncol(tab), function(j) 
	  {
		mean  <- tab[1, j]
		lower <- tab[2, j]
		upper <- tab[3, j]
		sprintf(paste0("%.", digits, "f (%.", digits, "f – %.", digits, "f)"), 
				mean, lower, upper)
	  })
	  names(formatted) <- colnames(tab)
	  as.data.frame(as.list(formatted), stringsAsFactors = FALSE, check.names = FALSE)
	}
 
	dent_A<-readRDS("dent_A.rds")
	dent_B<-readRDS("dent_B.rds")
	D1A<-format_ci_table(dent_A$all_ranges[1:3,],digits = 3)
	D1B<-format_ci_table(dent_B$all_ranges[1:3,],digits = 3)
		
	options(scipen=1000)
	write.table(D1A, "hOUwie_params.csv", sep = ",", row.names = FALSE, col.names = TRUE)
	write("\n---\n", file = "hOUwie_params.csv", append = TRUE)
	write.table(D1B, "hOUwie_params.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = TRUE)

#####################################################PHYLOGENETIC SIGNAL######################################################################

#calculating signal
data_name<-c("Pouch","Body mass","Litter size","Body mass x litter size",
			   "Pouch x body mass","Pouch x litter size","Pouch x body mass x litter size")
M_list<-vector("list",length(data_name))
names(M_list)<-data_name
M_results<-replicate(length(pst),M_list,simplify=FALSE)
for (i in 1:length(pst))
{	
	print(paste("dataset_tree",i))
	#preparing data and tree
	dat_M<-dat_100[[i]]
	dat_M$pouch<-as.factor(dat_M$pouch)
	TR<-pst[[i]]
	
	#obtaining Gower distances
	PC        <-gower_dist(dat_M["pouch"],type=list(factor="pouch"))
	BM        <-gower_dist(dat_M["BM"])
	LS        <-gower_dist(dat_M["LS"])
	BM_LS     <-gower_dist(dat_M[c("BM","LS")])
	PC_BM     <-gower_dist(dat_M[c("pouch","BM")],type=list(factor="pouch"))
	PC_LS     <-gower_dist(dat_M[c("pouch","LS")],type=list(factor="pouch"))
	PC_BM_LS  <-gower_dist(dat_M[c("pouch","BM","LS")],type=list(factor="pouch"))	

	#calculating phylogenetic signal
	M_results[[i]][[1]]<-phylosignal_M(trait_dist = PC,       phy = TR, reps = 999)
	M_results[[i]][[2]]<-phylosignal_M(trait_dist = BM,       phy = TR, reps = 999)
	M_results[[i]][[3]]<-phylosignal_M(trait_dist = LS,       phy = TR, reps = 999)
	M_results[[i]][[4]]<-phylosignal_M(trait_dist = BM_LS,    phy = TR, reps = 999)
	M_results[[i]][[5]]<-phylosignal_M(trait_dist = PC_BM,    phy = TR, reps = 999)
	M_results[[i]][[6]]<-phylosignal_M(trait_dist = PC_LS,    phy = TR, reps = 999)
	M_results[[i]][[7]]<-phylosignal_M(trait_dist = PC_BM_LS, phy = TR, reps = 999)
}
saveRDS(M_results,"M_results.rds")

for (i in 1:length(M_results))
	
#summarizing	
M_df<-as.data.frame(matrix(ncol=3,nrow=length(M_results[[1]])))
colnames(M_df)<-c("Data set","M statistic","p-value")
M_df[,1]<-data_name
M_summary<-replicate(length(M_results),M_df,simplify=FALSE)

for (i in 1:length(M_results))
{
	for (j in 1:length(data_name))
	{
		M_summary[[i]][j,][[2]]<-M_results[[i]][[j]]$stat
		M_summary[[i]][j,][[3]]<-M_results[[i]][[j]]$pvalue
	}
}
M_summary2<-do.call(rbind,M_summary)
M_summary3<-split(M_summary2,M_summary2[,1])
M_summary_final<-M_df[c(1:length(data_name)),]
for (i in 1:length(data_name))
{
	M_summary_final[i,2]<-mean_ci(M_summary3[[i]][,2],dec=2)
	M_summary_final[i,3]<-mean_ci(M_summary3[[i]][,3],dec=3)
}
M_summary_final
write.csv(M_summary_final,"Summary_M.csv",row.names=FALSE,fileEncoding = "UTF-8")


