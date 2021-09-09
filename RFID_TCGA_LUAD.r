/Users/jli19/Documents/G2P_Project/Cell_Type_Decomposition/NBF_analysis/BayCount_1

##preprocess

library(R.utils)
sourceDirectory("./")

Rcpp::sourceCpp('src/likelihood_NBFA_RF_cpp.cpp')
Rcpp::sourceCpp('src/update_pos_NBFA_RF_cpp.cpp')
Rcpp::sourceCpp('src/update_pos_cube_NBFA_RF_cpp.cpp')

hyper=NULL
hyper=NULL
hyper$eta = 0.1
hyper$a0 = 0.01
hyper$b0 = 0.01
hyper$e0 = 1
hyper$f0 = 1
hyper$g0 = 1
hyper$h0 = 1
hyper$u0 = hyper$v0 = 100
lab_switch.check = FALSE
B = 1000
nmc = 1000
K_range = 1:10
lab_switch.check = FALSE


flow_data<-read.delim("cd8_cd3p.txt",sep="\t",header=T)

signature_gene<-read.delim("signature_gene5260.txt",sep="\t",header=T)
getSimID1<-function(id)
{
	tmp<-unlist(strsplit(id,"\\."))
	return(tmp[1])
	
}



signature_ensemble<-as.vector(signature_gene[,2])
simid<-unlist(lapply(signature_ensemble,getSimID1))


ref_data<-read.delim("rma_20091gene_1194sample_742human_sample_profile_sortbyshark_gene.txt",sep="\t",header=T)
ref_sample<-as.vector(colnames(ref_data))
##ref_data_tumor<-read.delim("rma_20091gene_742_human_sample_profile_sortbyshark_gene.txt",sep="\t",header=T)

rma_normal_tissue<-read.delim("rma_human_hg133_158samples.txt",sep="\t",header=T)

f3<-read.delim("TCGA-LUAD.htseq_counts.tsv_526_tumor.txt",sep="\t",row.names=1)
tcga_luad_data<-2^f3-1
tcga_luad_data<-tcga_luad_data[rownames(tcga_luad_data)%in%signature_gene[,2],]
tcga_gene<-signature_gene[match(rownames(tcga_luad_data),signature_gene[,2]),1]
exp_tcga<-cbind(tcga_luad_data,tcga_gene)
 

rna_gene<-tcga_gene
ref_gene<-ref_data$Gene.Symbol


normal_gene<-as.vector(rownames(rma_normal_tissue))
shared_gene<-intersect(rna_gene,ref_gene)
shared_gene<-intersect(shared_gene,normal_gene)
 


getSimID<-function(id)
{
	tmp<-unlist(strsplit(id,"\\_"))
	return(tmp[1])
	
}

getIndex<-function(subexp,exp_rpkm)
{
	cor_vec<-vector()
	for(i in 1:length(subexp[1,]))
	{
		tmp<-vector()
		for(j in 1:length(exp_rpkm[1,]))
		{
			tmp<-c(tmp,cor(subexp[,i],exp_rpkm[,j]))
		}
		cor_vec<-c(cor_vec,sum(tmp))
	}
	return(which(cor_vec==max(cor_vec)))
}

simid<-unlist(lapply(ref_sample,getSimID))
simid<-gsub(".CEL.gz","",simid)
colnames(ref_data)<-simid

ref_data_sort<-ref_data[match(shared_gene,ref_data$Gene.Symbol),]
rma_normal_tissue_sort<-rma_normal_tissue[match(shared_gene,rownames(rma_normal_tissue)),]
exp_tcga_sort<-exp_tcga[match(shared_gene,exp_tcga$tcga_gene),]

gene_ref<-as.vector(ref_data_sort$Gene.Symbol)
tid<-as.vector(flow_data[,1])
tid<-paste(tid,"_T",sep="")


ref_exp<-ref_data_sort[,-1]
ref_exp<-ref_exp[,-1]
ref_exp<-ref_exp[,-1]
K=length(ref_exp[1,])

oo<-ref_exp

for(j in 1:length(ref_exp[,1]))
{
	## ref_exp[j, ] = ref_exp[j, ]/sum(ref_exp[j,])
}
 for (k in 1:K){
     ref_exp[, k] = ref_exp[, k]/sum(ref_exp[, k])
   }

rownames(ref_exp)<-as.vector(ref_data_sort[,1])


anno<-read.delim("sample_info_full.txt",sep="\t",header=F)


populations.names=c("T cells","CD4 T cells","NK cells","B cells",
"Monocytes","Myeloid Dendritic Cells","Neutrophils",
"Endothelium","Fibroblast","Lung")

populations.names=c("T cells","CD4 T cells","CD8 T cells","NK cells","B cells",
"Monocytes","Myeloid Dendritic Cells","Neutrophils",
"Endothelium","Fibroblast","Lung")

populations.names=c("T cells","CD4 T cells","CD8 T cells","NK cells","B cells",
"Monocytes","Myeloid Dendritic Cells","Neutrophils",
"Endothelium","Fibroblast","Lung")


features1=(anno[anno[,3]%in%populations.names,c(1,3)])
features2=split(as.vector(features1[,1]),as.vector(features1[,2]))

Phi<-matrix(0,nrow=length(ref_exp[,1]),ncol=length(features2))

for(i in 1:length(features2))
{
	subexp<-ref_exp[,colnames(ref_exp)%in%as.vector(unlist(features2[i]))]
	print(paste(length(subexp[1,]),names(features2[i])))
	##Phi[,i]<-rowMeans(subexp)
	index<-sample(1:length(subexp),1)
	Phi[,i]<-subexp[,index]
}


rownames(Phi)<-rownames(ref_exp)
colnames(Phi)<-c(names(features2))




Y<-exp_tcga_sort[,1:(length(exp_tcga_sort[1,])-1)]
rownames(Y)<-as.vector(exp_tcga_sort[,length(exp_tcga_sort[1,])])
Y<-as.matrix(Y)
G<-length(Y[,1])
Seeded_Phi=Phi
K_tmp=4
Init = list(Phi=matrix(1/G,nrow=G,ncol=K_tmp),Seeded_Phi=Seeded_Phi,Theta=matrix(15,nrow=K_tmp,ncol=length(Y[1,])))

B = 1000
nmc = 1000
d<- BayCount_MCMC(Y, Init, hyper, B, nmc, K_tmp, Print_Iter = FALSE, lab_switch.check)



purity_file<-read.delim("TCGA_mastercalls.abs_tables_JSedit.fixed.txt",sep="\t",header=T)
 

phi<-d$Phi[,,1999]
o<-d$Theta_norm[,,1000:1999]
pp_vec<-vector()
for(i in 1:length(Y[1,]))
{
	q<-(o[,i,1:1000])
	pp<-(rowMeans(q))
	pp_vec<-rbind(pp_vec,pp)
	
}

rownames(pp_vec)<-gsub("01A","01",colnames(Y))
rownames(pp_vec)<-gsub("\\.","\\-",rownames(pp_vec))
p<-purity_file[match(as.vector(rownames(pp_vec)),purity_file$array),]
sample_hit<-as.vector(p[!is.na(p$purity),1])


pdf("TCGA_GBM_purity.pdf")


r1<-pp_vec[match(sample_hit,rownames(pp_vec)),2]
r2<-p[!is.na(p$purity),]$purity
r1<-as.vector(unlist(r1))
cor(r1,r2)

test<-cbind(r1,r2)
##test<-test[!(test[,1]<0.4&(test[,2]>0.8)),]
##test<-test[!(test[,1]>0.6&(test[,2]<0.2)),]
##test<-test[!((test[,1]<0.05)),]
summary(lm(test[,1]~test[,2]))


##r1<-as.vector(unlist(r1))
##ddd<-cbind(r1,r2)

##test<-ddd

reg1<-lm(test[,1]~test[,2])
plot(test[,2],test[,1],xlab="RNA-based purity",ylab="Absolute Predicted Purity",main="GBM",xlim=c(0,1),ylim=c(0,1))
 reg1<-lm(test[,1]~test[,2])
 abline(reg1,col="red")
 anno<-summary(reg1)
 	 
 	 r2=anno$adj.r.squared
 	 my.p=anno$coefficients[2,4]
 	 rp = vector('expression',3)
 	 val<-cor(test[,1],test[,2])+0.05
 	  rp[1] =substitute(expression(pearson_correlation == MYOTHERVALUE), 
		list(MYOTHERVALUE = format(val, digits = 2)))[2]
	 rp[2] = substitute(expression(italic(R)^2 == MYVALUE), 
		list(MYVALUE = format(r2,dig=3)))[2]
	rp[3] = substitute(expression(italic(p) == MYOTHERVALUE), 
		list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
 	 legend('bottomright', legend = rp, bty = 'n')
 	
 	
 	dev.off()
 	
 	

gettumorsimPID<-function(id_all)
{
				l<-unlist(strsplit(id_all,"\\-"))
				patient_id<-paste(l[1],l[2],l[3],sep="-")
				return((l[4]))
}	
fm<-list.files("./",pattern="htseq_counts.tsv$")
for(k in 1:length(fm))
{
fm_k<-fm[k]



		f<-read.delim(fm_k,sep="\t",header=T,row.names=1,check.names=F)
		nm<-as.vector(colnames(f))
		tumor_sample<-vector()
		normal_sample<-vector()
		
		for(j in 1:length(nm))
		{
		
			name1=gettumorsimPID(nm[j])
			name1<-substr(name1,1,2)
			num<-as.integer(name1)
			if(num<10)
			{
				tumor_sample<-c(tumor_sample,nm[j])
			}else if(num<20)
			{
				normal_sample<-c(normal_sample,nm[j])
			}
		}
		
		if(length(tumor_sample)>0)
		{
			tumor_file=paste(fm_k,"_",length(tumor_sample),"_tumor.txt",sep="")
			data<-f[,tumor_sample]
			Ensembl_ID<-as.vector(rownames(f))
			data<-cbind(Ensembl_ID,data)
			write.table(data,file=tumor_file,col.names=T,row.names=F,sep="\t",quote=F)
		}
		if(length(normal_sample)>0)
		{
			normal_file=paste(fm_k,"_",length(normal_sample),"_normal.txt",sep="")
			data<-f[,normal_sample]
			Ensembl_ID<-as.vector(rownames(f))
			data<-cbind(Ensembl_ID,data)
			write.table(data,file=normal_file,col.names=T,row.names=F,sep="\t",quote=F)
		}



}

