############################################################################################################################################################
#script: marsupium evolution x body mass x litter size
#author: Daniel Casali
############################################################################################################################################################

#Load packages
packages<-c("evobiR","phytools","dplyr","tidyr","openxlsx","ggplot2","ggpubr","phyr","OUwie","car","DHARMa","plotrix","svglite","caper")
sapply(packages, function(x) library(x, character.only = TRUE))

#Load data and filter to keep only taxa with both body mass and litter size information
dat<-read.csv("marsupial_data.csv",row.names=1)

#Load major clade credibility tree (from Upham et al. 2019), prunned to taxomic sample in data
mcc<-read.tree("marsupial_MCC.nwk")
mcc

#Load 100 posterior sample trees (from Upham et al. 2019), prunned to taxomic sample in data
#tree-pruner-79f70b71-c487-403a-94a8-2d5bbe196bb2
pst<-read.tree("marsupial_PST.trees")
pst

###################################Phylogenetic generalized liner mixed models (PGLMM)######################################################################
#BICw function
BICw <- function(BIC_values) 
{
  delta <- BIC_values - min(BIC_values)
  weights <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))
  return(weights)
}

#Calculate the variance inflation factor (VIF) for body mass and litter size
dat2<-dat
dat2$pouch[dat2$pouch=="absent"]<-0
dat2$pouch[dat2$pouch=="present"]<-1
dat2$pouch<-as.numeric(dat2$pouch)
dat2$mass<-log10(dat2$mass)
dat2$zmass<-(dat2$mass-mean(dat2$mass))/sd(dat2$mass)
dat2$zlitter<-(dat2$litter-mean(dat2$litter))/sd(dat2$litter)
reg <- glm(pouch~zmass*zlitter, data = dat2, family = binomial)
VIF<-round(vif(reg),2)
VIF

#Biplot
dat2<-dat
dat2$mass<-log10(dat2$mass)
dat2$zmass<-(dat2$mass-mean(dat2$mass))/sd(dat2$mass)
dat2$zlitter<-(dat2$litter-mean(dat2$litter))/sd(dat2$litter)		
ggplot(data=dat2,aes(x=zmass,y=zlitter,color=pouch,shape=pouch)) + geom_point(alpha=0.75, size=4) +
theme_classic() + theme(legend.position=c(0.85,0.9),legend.title = element_text(face = "bold")) + 
labs (x="Body mass (Z-scores)", y="Litter size (Z-scores)", title="") +
scale_color_viridis_d()+
scale_shape_manual(values=c(19,19))+
guides(color = guide_legend("Pouch"), shape = guide_legend("Pouch")) +
#annotate("text", hjust=0, x = 1.15, y = 2, label = "VIF", color = "black", size = 4, fontface="bold")+
#annotate("text", hjust=0, x = 1.15, y = 1.8, label = paste("mass -", VIF[[1]]), color = "black", size = 3)+
#annotate("text", hjust=0, x = 1.15, y = 1.6, label = paste("litter -", VIF[[2]]), color = "black", size = 3)+
#annotate("text", hjust=0, x = 1.15, y = 1.4, label = paste("mass:litter -", VIF[[3]]), color = "black", size = 3)
ggsave("Biplot.svg",width=5,height=5)
	
#Perform logistic regression - body mass and litter size as fixed effects and phylogeny as a random effect

dat2$pouch[dat2$pouch=="absent"]<-0
dat2$pouch[dat2$pouch=="present"]<-1
dat2$pouch<-as.numeric(dat2$pouch)

PGLMM_results<-list()	
for (i in 1:length(pst))
{
	print(paste("Tree",i,"of",length(pst)))
	TR<-pst[[i]]
	dat3<-ReorderData(TR,dat2)
	null<-pglmm_compare(pouch~1,data=dat3,family="binomial",phy=TR)
	mass<-pglmm_compare(pouch~zmass,data=dat3,family="binomial",phy=TR)
	litt<-pglmm_compare(pouch~zlitter,data=dat3,family="binomial",phy=TR)
	both<-pglmm_compare(pouch~zmass+zlitter,data=dat3,family="binomial",phy=TR)
	inte<-pglmm_compare(pouch~zmass*zlitter,data=dat3,family="binomial",phy=TR)
	PGLMM_results[[i]]<-list(null,mass,litt,both,inte)		
}
saveRDS(PGLMM_results,"PGLMM_results.rds")

	#Summarize results	
	BIC_PGLMM<-as.data.frame(matrix(nrow=length(PGLMM_results[[1]]),ncol=length(PGLMM_results)))
	rownames(BIC_PGLMM)<-c("Null","Body mass","Litter size","Body mass + Litter Size","Body mass * Litter Size")
	dd<-dim(BIC_PGLMM)
	BICw_PGLMM<-BIC_PGLMM
	s2_PGLMM<-BIC_PGLMM
	par_PGLMM<-BIC_PGLMM[,c(1:4)]
	colnames(par_PGLMM)<-c("Intercept","Body mass","Litter size","Body mass:Litter Size")
	PAR_PGLMM<-replicate(length(pst),par_PGLMM,simplify=FALSE)
	for (i in 1:dd[2])
	{
	for (j in 1:dd[1])
	{
		BIC_PGLMM[j,i]<-PGLMM_results[[i]][[j]]$BIC
		s2_PGLMM[j,i]<-PGLMM_results[[i]][[j]]$ss
		params<-as.vector(PGLMM_results[[i]][[j]]$B)
		PAR_PGLMM[[i]][j,]<-c(params, rep(NA,(4-length(params))))
	}		
	BICw_PGLMM[,i]<-round(BICw(BIC_PGLMM[,i]),2)
	}	
	
	##BIC
	BIC1=data.frame(median=round(apply(BIC_PGLMM, 1, median), 2),
				lower=round(apply(BIC_PGLMM, 1, quantile, 0.025), 2), 
				upper=round(apply(BIC_PGLMM, 1, quantile, 0.975), 2))
	BIC2=unite(BIC1, "CI", c(lower,upper), sep=" - ")
	BIC2$prov1<-" ("
	BIC2$prov2<-") "
	BIC3=unite(BIC2, "BIC", c(median,prov1,CI,prov2), sep="")
	BIC3
	
	##BICw
	BICw1=data.frame(median=round(apply(BICw_PGLMM, 1, median), 2),
				lower=round(apply(BICw_PGLMM, 1, quantile, 0.025), 2), 
				upper=round(apply(BICw_PGLMM, 1, quantile, 0.975), 2))
	BICw2=unite(BICw1, "CI", c(lower,upper), sep=" - ")
	BICw2$prov1<-" ("
	BICw2$prov2<-") "
	BICw3=unite(BICw2, "BICw", c(median,prov1,CI,prov2), sep="")
	BICw3
	
	##s2
	s2_phy=data.frame(median=round(apply(s2_PGLMM, 1, median), 2),
				lower=round(apply(s2_PGLMM, 1, quantile, 0.025), 2), 
				upper=round(apply(s2_PGLMM, 1, quantile, 0.975), 2))
	s2_phy2=unite(s2_phy, "CI", c(lower,upper), sep=" - ")
	s2_phy2$prov1<-" ("
	s2_phy2$prov2<-") "
	s2_phy3=unite(s2_phy2, "s2.phy", c(median,prov1,CI,prov2), sep="")
	s2_phy3
	
	##Model parameters (intercept and slope)
	ss<-dim(par_PGLMM)
	PAR<-par_PGLMM
	for (i in 1:ss[[1]])
	{
		for (j in 1:ss[[2]])
		{
			par_md<-round(quantile(sapply(PAR_PGLMM, function(df) df[i, j]),0.500, na.rm=TRUE)[[1]],2)
			par_q1<-round(quantile(sapply(PAR_PGLMM, function(df) df[i, j]),0.025, na.rm=TRUE)[[1]],2)
			par_q2<-round(quantile(sapply(PAR_PGLMM, function(df) df[i, j]),0.975, na.rm=TRUE)[[1]],2)
			
			if(!is.na(par_md)) { PAR[i,j]<-paste0(par_md," (",par_q1," - ",par_q2,") ")
			} else { PAR[i,j]<-"NA"
			}
		}
	}
	PAR
	
	##Joining everything...
	PGLMM_summary<-data.frame(PAR,s2_phy3,BIC3,BICw3)
	PGLMM_summary	
	write.csv(PGLMM_summary,"Summary_PGLMM.csv")

#Check residuals
pdf("Quality_check.pdf")
for (i in 1:length(PGLMM_results))
{
	for (j in 1:length(PGLMM_results[[i]]))
	{
		print(paste0("Tree - ",i,", Model - ",rownames(BIC_PGLMM)[j]))
		class(PGLMM_results[[i]][[j]])<-"communityPGLMM"
		plot(simulateResiduals(PGLMM_results[[i]][[j]]))
		title(paste0("Tree - ",i,", Model - ",rownames(BIC_PGLMM)[j]), outer=TRUE, line=-1)
	}
}
dev.off()

#####################################################DISCRETE + CONTINUOUS MODELS######################################################################

#BODY MASS - Prepare dataset and set the number of stochastic maps (nSim)
dataset<-data.frame(taxon=rownames(dat),pouch=dat$pouch,mass=log10(dat$mass))
nSim<-100

	##Character-dependent models (CD, Boyko et al. 2023)	
	CD_ER_BMV    <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "BMV",   null.model=FALSE, nSim=nSim)
	CD_ER_OUV    <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "OUV",   null.model=FALSE, nSim=nSim)
	CD_ER_OUM    <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "OUM",   null.model=FALSE, nSim=nSim)
	CD_ER_OUMV   <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "OUMV",  null.model=FALSE, nSim=nSim)
	CD_ARD_BMV   <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "BMV",   null.model=FALSE, nSim=nSim)
	CD_ARD_OUV   <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "OUV",   null.model=FALSE, nSim=nSim)
	CD_ARD_OUM   <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "OUM",   null.model=FALSE, nSim=nSim)
	CD_ARD_OUMV  <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "OUMV",  null.model=FALSE, nSim=nSim)
	
	##Character-independent models (CID,  Boyko et al. 2023)	
	CID_ER_BM1    <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "BM1",   null.model=FALSE, nSim=nSim)
	CID_ER_OU1    <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "OU1",   null.model=FALSE, nSim=nSim)
	CID_ARD_BM1   <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "BM1",   null.model=FALSE, nSim=nSim)
	CID_ARD_OU1   <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "OU1",   null.model=FALSE, nSim=nSim)	
	
	##Character-independent plus models (CID+, Boyko et al. 2023)
	CIDP_ER_BM1    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "BM1",   null.model=TRUE, nSim=nSim)
	CIDP_ER_BMV    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "BMV",   null.model=TRUE, nSim=nSim)
	CIDP_ER_OU1    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OU1",   null.model=TRUE, nSim=nSim)
	CIDP_ER_OUV    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUV",   null.model=TRUE, nSim=nSim)
	CIDP_ER_OUM    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUM",   null.model=TRUE, nSim=nSim)
	CIDP_ER_OUMV   <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUMV",  null.model=TRUE, nSim=nSim)
	CIDP_ARD_BM1   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "BM1",   null.model=TRUE, nSim=nSim)
	CIDP_ARD_BMV   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "BMV",   null.model=TRUE, nSim=nSim)
	CIDP_ARD_OU1   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OU1",   null.model=TRUE, nSim=nSim)
	CIDP_ARD_OUV   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUV",   null.model=TRUE, nSim=nSim)
	CIDP_ARD_OUM   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUM",   null.model=TRUE, nSim=nSim)
	CIDP_ARD_OUMV  <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUMV",  null.model=TRUE, nSim=nSim)
	
	##Hybrid models (HYB, Boyko et al. 2023)
	HYB_ER_BMV    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "BMV",   null.model=FALSE, nSim=nSim)
	HYB_ER_OUV    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUV",   null.model=FALSE, nSim=nSim)
	HYB_ER_OUM    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUM",   null.model=FALSE, nSim=nSim)
	HYB_ER_OUMV   <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUMV",  null.model=FALSE, nSim=nSim)
	HYB_ARD_BMV   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "BMV",   null.model=FALSE, nSim=nSim)
	HYB_ARD_OUV   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUV",   null.model=FALSE, nSim=nSim)
	HYB_ARD_OUM   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUM",   null.model=FALSE, nSim=nSim)
	HYB_ARD_OUMV  <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUMV",  null.model=FALSE, nSim=nSim)

	##Save model results
	model_results<-list(CD_ER_BMV, CD_ER_OUV, CD_ER_OUM, CD_ER_OUMV,  
						CD_ARD_BMV, CD_ARD_OUV, CD_ARD_OUM, CD_ARD_OUMV,  
						CID_ER_BM1, CID_ER_OU1, CID_ARD_BM1, CID_ARD_OU1, 
						CIDP_ER_BM1, CIDP_ER_BMV, CIDP_ER_OU1, CIDP_ER_OUV, CIDP_ER_OUM, CIDP_ER_OUMV, 
						CIDP_ARD_BM1, CIDP_ARD_BMV, CIDP_ARD_OU1, CIDP_ARD_OUV, CIDP_ARD_OUM, CIDP_ARD_OUMV,
						HYB_ER_BMV, HYB_ER_OUV, HYB_ER_OUM, HYB_ER_OUMV, 
						HYB_ARD_BMV, HYB_ARD_OUV, HYB_ARD_OUM, HYB_ARD_OUMV)
						
	model_names<-c("CD_ER_BMV", "CD_ER_OUV", "CD_ER_OUM", "CD_ER_OUMV", 
				   "CD_ARD_BMV", "CD_ARD_OUV", "CD_ARD_OUM", "CD_ARD_OUMV",
				   "CID_ER_BM1", "CID_ER_OU1", "CID_ARD_BM1", "CID_ARD_OU1",
				   "CIDP_ER_BM1", "CIDP_ER_BMV", "CIDP_ER_OU1", "CIDP_ER_OUV", "CIDP_ER_OUM", "CIDP_ER_OUMV",
				   "CIDP_ARD_BM1", "CIDP_ARD_BMV", "CIDP_ARD_OU1", "CIDP_ARD_OUV", "CIDP_ARD_OUM", "CIDP_ARD",
				   "HYB_ER_BMV", "HYB_ER_OUV", "HYB_ER_OUM", "HYB_ER_OUMV",
				   "HYB_ARD_BMV", "HYB_ARD_OUV", "HYB_ARD_OUM", "HYB_ARD_OUMV")
	names(model_results)<-model_names
	saveRDS(model_results,"model_results1.rds")

	##Sumarize BIC
	model_results<-readRDS("model_results1.rds")
	models_BIC<-as.data.frame(matrix(nrow=length(model_results),ncol=2))
	colnames(models_BIC)<-c("Model","BIC")
	for (i in 1:length(model_results))
	{
		models_BIC[i,1]<-names(model_results)[i]
		models_BIC[i,2]<-round(model_results[[i]]$BIC,2)
	}
	min_BIC <- min(models_BIC$BIC)
	delta_BIC <- models_BIC$BIC - min_BIC
	BIC_weights <- exp(-0.5 * delta_BIC) / sum(exp(-0.5 * delta_BIC))
	models_BIC$BICw<-round(BIC_weights,2)
	models_BIC
	write.csv(models_BIC,"models_BIC1.csv")
	
#LITTER SIZE - Prepare dataset and set the number of stochastic maps (nSim)
dataset<-data.frame(taxon=rownames(dat),pouch=dat$pouch,litter=log10(dat$litter))
nSim<-100

	##Character-dependent models (CD, Boyko et al. 2023)	
	CD_ER_BMV    <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "BMV",   null.model=FALSE, nSim=nSim)
	CD_ER_OUV    <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "OUV",   null.model=FALSE, nSim=nSim)
	CD_ER_OUM    <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "OUM",   null.model=FALSE, nSim=nSim)
	CD_ER_OUMV   <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "OUMV",  null.model=FALSE, nSim=nSim)
	CD_ARD_BMV   <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "BMV",   null.model=FALSE, nSim=nSim)
	CD_ARD_OUV   <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "OUV",   null.model=FALSE, nSim=nSim)
	CD_ARD_OUM   <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "OUM",   null.model=FALSE, nSim=nSim)
	CD_ARD_OUMV  <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "OUMV",  null.model=FALSE, nSim=nSim)
	
	##Character-independent models (CID,  Boyko et al. 2023)	
	CID_ER_BM1    <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "BM1",   null.model=FALSE, nSim=nSim)
	CID_ER_OU1    <-hOUwie(mcc, dataset, rate.cat=1, "ER",  "OU1",   null.model=FALSE, nSim=nSim)
	CID_ARD_BM1   <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "BM1",   null.model=FALSE, nSim=nSim)
	CID_ARD_OU1   <-hOUwie(mcc, dataset, rate.cat=1, "ARD", "OU1",   null.model=FALSE, nSim=nSim)	
	
	##Character-independent plus models (CID+, Boyko et al. 2023)
	CIDP_ER_BM1    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "BM1",   null.model=TRUE, nSim=nSim)
	CIDP_ER_BMV    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "BMV",   null.model=TRUE, nSim=nSim)
	CIDP_ER_OU1    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OU1",   null.model=TRUE, nSim=nSim)
	CIDP_ER_OUV    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUV",   null.model=TRUE, nSim=nSim)
	CIDP_ER_OUM    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUM",   null.model=TRUE, nSim=nSim)
	CIDP_ER_OUMV   <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUMV",  null.model=TRUE, nSim=nSim)
	CIDP_ARD_BM1   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "BM1",   null.model=TRUE, nSim=nSim)
	CIDP_ARD_BMV   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "BMV",   null.model=TRUE, nSim=nSim)
	CIDP_ARD_OU1   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OU1",   null.model=TRUE, nSim=nSim)
	CIDP_ARD_OUV   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUV",   null.model=TRUE, nSim=nSim)
	CIDP_ARD_OUM   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUM",   null.model=TRUE, nSim=nSim)
	CIDP_ARD_OUMV  <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUMV",  null.model=TRUE, nSim=nSim)
	
	##Hybrid models (HYB, Boyko et al. 2023)
	HYB_ER_BMV    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "BMV",   null.model=FALSE, nSim=nSim)
	HYB_ER_OUV    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUV",   null.model=FALSE, nSim=nSim)
	HYB_ER_OUM    <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUM",   null.model=FALSE, nSim=nSim)
	HYB_ER_OUMV   <-hOUwie(mcc, dataset, rate.cat=2, "ER",  "OUMV",  null.model=FALSE, nSim=nSim)
	HYB_ARD_BMV   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "BMV",   null.model=FALSE, nSim=nSim)
	HYB_ARD_OUV   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUV",   null.model=FALSE, nSim=nSim)
	HYB_ARD_OUM   <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUM",   null.model=FALSE, nSim=nSim)
	HYB_ARD_OUMV  <-hOUwie(mcc, dataset, rate.cat=2, "ARD", "OUMV",  null.model=FALSE, nSim=nSim)

	##Save model results
	model_results<-list(CD_ER_BMV, CD_ER_OUV, CD_ER_OUM, CD_ER_OUMV,  
						CD_ARD_BMV, CD_ARD_OUV, CD_ARD_OUM, CD_ARD_OUMV,  
						CID_ER_BM1, CID_ER_OU1, CID_ARD_BM1, CID_ARD_OU1, 
						CIDP_ER_BM1, CIDP_ER_BMV, CIDP_ER_OU1, CIDP_ER_OUV, CIDP_ER_OUM, CIDP_ER_OUMV, 
						CIDP_ARD_BM1, CIDP_ARD_BMV, CIDP_ARD_OU1, CIDP_ARD_OUV, CIDP_ARD_OUM, CIDP_ARD_OUMV,
						HYB_ER_BMV, HYB_ER_OUV, HYB_ER_OUM, HYB_ER_OUMV, 
						HYB_ARD_BMV, HYB_ARD_OUV, HYB_ARD_OUM, HYB_ARD_OUMV)
						
	model_names<-c("CD_ER_BMV", "CD_ER_OUV", "CD_ER_OUM", "CD_ER_OUMV", 
				   "CD_ARD_BMV", "CD_ARD_OUV", "CD_ARD_OUM", "CD_ARD_OUMV",
				   "CID_ER_BM1", "CID_ER_OU1", "CID_ARD_BM1", "CID_ARD_OU1",
				   "CIDP_ER_BM1", "CIDP_ER_BMV", "CIDP_ER_OU1", "CIDP_ER_OUV", "CIDP_ER_OUM", "CIDP_ER_OUMV",
				   "CIDP_ARD_BM1", "CIDP_ARD_BMV", "CIDP_ARD_OU1", "CIDP_ARD_OUV", "CIDP_ARD_OUM", "CIDP_ARD",
				   "HYB_ER_BMV", "HYB_ER_OUV", "HYB_ER_OUM", "HYB_ER_OUMV",
				   "HYB_ARD_BMV", "HYB_ARD_OUV", "HYB_ARD_OUM", "HYB_ARD_OUMV")
	names(model_results)<-model_names
	saveRDS(model_results,"model_results2.rds")

	##Sumarize BIC
	model_results2<-readRDS("model_results2.rds")
	models_BIC2<-as.data.frame(matrix(nrow=length(model_results2),ncol=2))
	colnames(models_BIC2)<-c("Model","BIC")
	for (i in 1:length(model_results2))
	{
		models_BIC2[i,1]<-names(model_results2)[i]
		models_BIC2[i,2]<-round(model_results2[[i]]$BIC,2)
	}
	min_BIC <- min(models_BIC2$BIC)
	delta_BIC <- models_BIC2$BIC - min_BIC
	BIC_weights <- exp(-0.5 * delta_BIC) / sum(exp(-0.5 * delta_BIC))
	models_BIC2$BICw<-round(BIC_weights,2)
	models_BIC2
	write.csv(models_BIC2,"models_BIC2.csv")
	
#Ancestral state estimations	
	
	##Model set 1 - pouch x bodymass
	model_results1<-readRDS("model_results1.rds")
	m1<-hOUwie.recon(model_results1$CD_ER_BMV,nodes = "internal")
	saveRDS(m1,"m1_ase.rds")
	
	##Model set 2 - pouch x litter size
	model_results2<-readRDS("model_results2.rds")
	m2<-hOUwie.recon(model_results2$CD_ER_BMV,nodes = "internal")
	saveRDS(m2,"m2_ase.rds")
	
	#PLOTS	
	m1<-readRDS("m1_ase.rds")
	m2<-readRDS("m2_ase.rds")	
	tip_data<-as.factor(setNames(dat$pouch,rownames(dat)))
	cols<-make.transparent(c("#440154FF","#FDE725FF"),alpha=0.9)	
	pdf("ASR.pdf",width=8,height=9.5)		
		par(mfrow=c(2,1))
		arc_height<-0.5
		plotTree(mcc, type="arc", fsize=0.35, ftype="i", offset = 3, lwd=0.3, arc_height=arc_height)
		title(main="", adj=0.35, line = -1.0, cex.main=1.5)
		par(lty="solid",fg="transparent")
		nodelabels(node=1:mcc$Nnode+Ntip(mcc), pie=m1, piecol=cols, cex=0.20)
		tiplabels(pie=to.matrix(tip_data[mcc$tip.label],levels(tip_data)),piecol=cols,cex=0.15)
		par(lty="solid",fg="black")	
		legend(x=-20,y=25,c("absent","present"),cex=0.75, pt.cex=1.3, pch=21, pt.bg = cols, col = "white", bty="n",title="Pouch")
		h=max(nodeHeights(mcc))
		labs<-seq(0,h,by=10)
		a1<-axis(1,pos=-0.02*h,at=h-labs+arc_height*h,
		labels=labs,cex.axis=0.7,lwd=1,lend=1,padj=-2)
		text(mean(a1),-0.20*h,"million years ago",font=1, cex=0.8)
		a2<-axis(1,pos=-0.02*h,at=-h+labs-arc_height*h,
		labels=labs,cex.axis=0.7,lwd=1,lend=1,padj=-2)
		text(mean(a2),-0.20*h,"million years ago",font=1, cex=0.8)
		draw.arc(0,0,radius=h-labs[2:length(labs)]+arc_height*h,
		angle1=,angle2=pi,col=make.transparent("grey",1),lty="dotted")
		title(paste("a)","Correlated with body mass"),adj=0.1,line=-2.5,cex.main=0.9)
		
		plotTree(mcc, type="arc", fsize=0.35, ftype="i", offset = 3, lwd=0.3, arc_height=arc_height)
		title(main="", adj=0.35, line = -1.0, cex.main=1.5)
		par(lty="solid",fg="transparent")
		nodelabels(node=1:mcc$Nnode+Ntip(mcc), pie=m2, piecol=cols, cex=0.20)
		tiplabels(pie=to.matrix(tip_data[mcc$tip.label],levels(tip_data)),piecol=cols,cex=0.15)
		par(lty="solid",fg="black")	
		legend(x=-20,y=25,c("absent","present"),cex=0.75, pt.cex=1.3, pch=21, pt.bg = cols, col = "white", bty="n",title="Pouch")
		h=max(nodeHeights(mcc))
		labs<-seq(0,h,by=10)
		a1<-axis(1,pos=-0.02*h,at=h-labs+arc_height*h,
		labels=labs,cex.axis=0.7,lwd=1,lend=1,padj=-2)
		text(mean(a1),-0.20*h,"million years ago",font=1, cex=0.8)
		a2<-axis(1,pos=-0.02*h,at=-h+labs-arc_height*h,
		labels=labs,cex.axis=0.7,lwd=1,lend=1,padj=-2)
		text(mean(a2),-0.20*h,"million years ago",font=1, cex=0.8)
		draw.arc(0,0,radius=h-labs[2:length(labs)]+arc_height*h,
		angle1=,angle2=pi,col=make.transparent("grey",1),lty="dotted")
		title(paste("b)","Correlated with litter size"),adj=0.1,line=-2.5,cex.main=0.9)
	dev.off()
	
	
#####################################################PHYLOGENETIC SIGNAL######################################################################

#Prepare data and calculate phylogenetic signal (Blombergs's K for continuous predictors and Fritz and Purvis's D for pouch presence)
bodymass<-setNames(log10(dat$mass),rownames(dat))
litter<-setNames(dat$litter,rownames(dat))
pouch_data<-data.frame(taxa=rownames(dat),dat$pouch)
BM_K<-LS_K<-P_D<-list()
for (i in 1:length(pst))
{
	print(paste("tree",i))
	BM_K[[i]]<-phylosig(pst[[i]], bodymass, method="K", test=TRUE, nsim=1000)
	LS_K[[i]]<-phylosig(pst[[i]], litter, method="K", test=TRUE, nsim=1000)
		
	comp_data<-comparative.data(pst[[i]], pouch_data, names.col = taxa)
	P_D[[i]]<-phylo.d(data=comp_data, binvar = dat.pouch, permut = 1000)
}
signal_results<-list(BM_K=BM_K, LS_K=LS_K, P_D=P_D)
saveRDS(signal_results,"signal_results.rds")

#Summarize stats and p-values and export
phylosig_summary<-list(BM_K=c(),BM_Pvalue=c(),LS_K=c(),LS_Pvalue=c(),P_D=c(),P_Random=c(),P_Brownian=c())
for (i in 1:length(pst))
{
	phylosig_summary[[1]][i]<-signal_results[[1]][[i]][[1]]
	phylosig_summary[[2]][i]<-signal_results[[1]][[i]][[2]]
	phylosig_summary[[3]][i]<-signal_results[[2]][[i]][[1]]
	phylosig_summary[[4]][i]<-signal_results[[2]][[i]][[2]]
	phylosig_summary[[5]][i]<-signal_results[[3]][[i]]$DEstimate[[1]]
	phylosig_summary[[6]][i]<-signal_results[[3]][[i]]$Pval1
	phylosig_summary[[7]][i]<-signal_results[[3]][[i]]$Pval0
}

summary_median_quantiles <- function(x) {
  q <- quantile(x, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
  sprintf("%.2f (%.2f–%.2f)", q[2], q[1], q[3])
}
result <- sapply(phylosig_summary, summary_median_quantiles)
result
write.table(result,"physignal_summary.txt")