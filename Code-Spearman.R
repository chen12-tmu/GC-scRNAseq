

#install.packages("corrplot")

library(corrplot)


M=cor(data,method = 'spearman')
res1=cor.mtest(data, conf.level = 0.95)


pdf(file="cor.pdf", width=8, height=8)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=0.9, pch=T,
         insig = "label_sig",
         pch.cex = 0,
         sig.level=0.05,
         number.cex = 0.1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()


