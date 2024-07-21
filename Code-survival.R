

#install.packages('survival')
#install.packages("survminer")



library(survival)
library(survminer)

coxPfilter=0.05                   
inputFile="TCGA.expTime.txt"      
setwd("C:\\Users\\chen\\Desktop\\单细胞亚群分析\\免疫")    


rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)


outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
	#cox分析
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	
	if(coxP<coxPfilter){
	    sigGenes=c(sigGenes,i)
		outTab=rbind(outTab,
			         cbind(id=i,
			         HR=coxSummary$conf.int[,"exp(coef)"],
			         HR.95L=coxSummary$conf.int[,"lower .95"],
			         HR.95H=coxSummary$conf.int[,"upper .95"],
			         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
			        )
		
	
		data=rt[,c("futime", "fustat", i)]
		colnames(data)=c("futime", "fustat", "gene")
		
		res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
		res.cat=surv_categorize(res.cut)
	
		diff=survdiff(Surv(futime, fustat) ~gene,data =res.cat)
		pValue=1-pchisq(diff$chisq, df=1)
		#print(pValue)
		if(pValue<0.001){
			pValue="p<0.001"
		}else{
			pValue=paste0("p=",sprintf("%.03f",pValue))
		}
	
		fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)
		surPlot=ggsurvplot(fit,
						data=res.cat,
						pval=pValue,
						pval.size=6,
						legend.title=i,
						legend.labs=c("high","low"),
						xlab="Time(years)",
						palette=c("#FF34B3", "#AB82FF"),
						break.time.by=2,
						conf.int=F,
						risk.table=TRUE,
						risk.table.title="",
						risk.table.height=.25)
		pdf(file=paste0("survival.",i,".pdf"), onefile=FALSE, width=6, height =7)
		print(surPlot)
		dev.off()
	}
}
