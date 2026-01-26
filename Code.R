#####Phenotypic data visualization and statistical analysis
rm(list=ls())
data<-read.csv("data environment.csv")
group_data<-read.csv("Sample_Treatment.csv",header = T)
a<- merge(data,group_data, by="sample")
library("ggplot2")
library("ggpubr")
a$Treatments<- factor(a$Treatments,levels=c("T1","T2","S1","S2"))
p<-ggplot(data =a,aes(x = Treatments,y =Dry.weight,color = Treatments)) +
  stat_boxplot(geom="errorbar",size=1,width=0.3,linetype="solid")+
  geom_boxplot()+
  geom_jitter(width=0.2,size=2)+
  scale_color_manual(values=c("#293A6A","#3C8A84","#8B8F9E","#344D8B"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="",y="Dry weight of the plants (g/plant)")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),axis.text.x=element_text(angle=45,vjust = 0.5,hjust = 0.5,colour='black',size=12))
p
ggsave("Dry weight.pdf",p,height=5,width=8)

fit <- aov(Dry.weight~soil_s, a)
summary(fit)

library(agricolae)
LSD1<-LSD.test(fit,"soil_s",p.adj = "none")
LSD1




###PCoA+PERMANOVA
library(vegan)
library(ape)
library(ggplot2)
library(grid)
library(dplyr)
library(multcomp)
library(patchwork)

data<-read.csv("OTU.csv",row.names = 1)
groups<-read.csv("Sample_Treatment.csv")

cbbPalette <- c("#293A6A","#748BBA","#3C8A84","#CBE1E3")

#PCoA
data <- vegdist(data)
pcoa<- pcoa(data, correction = "none", rn = NULL)
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2,groups$soil_s)
colnames(plotdata) <-c("sample","PC1","PC2","soil_s")
pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)
plotdata$soil_s <- factor(plotdata$soil_s,levels = c("T1","T2","S1","S2"))

#PC1PC2
yf <- plotdata
yd1 <- yf %>% group_by(soil_s) %>% summarise(Max = max(PC1))
yd2 <- yf %>% group_by(soil_s) %>% summarise(Max = max(PC2))
yd1$Max <- yd1$Max + max(yd1$Max)*0.1
yd2$Max <- yd2$Max + max(yd2$Max)*0.1
fit1 <- aov(PC1~soil_s,data = plotdata)
tuk1<-glht(fit1,linfct=mcp(soil_s="Tukey"))
res1 <- cld(tuk1,alpah=0.05)
fit2 <- aov(PC2~soil_s,data = plotdata)
tuk2<-glht(fit2,linfct=mcp(soil_s="Tukey"))
res2 <- cld(tuk2,alpah=0.05)
test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,
                   yd1 = yd1$Max,yd2 = yd2$Max,soil_s = yd1$soil_s)
test$soil_s <- factor(test$soil_s,levels = c("T1","T2","S1","S2"))

#Plot
p<-ggplot(plotdata, aes(PC1, PC2)) +
  scale_color_manual(values=cbbPalette,name = "soil_s")+
  geom_point(size=8,aes(color=groups$soil_s))+
  stat_ellipse(aes(fill = groups$soil_s),level = 0.7,geom = "polygon", alpha = 0.05, show.legend = F) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank())
p

#PERMANOVA
otu.adonis=adonis(data~soil_s,data = groups,distance = "bray")
pvalue<-otu.adonis$aov.tab



####Microbial difference analysis
library(ggplot2)
library(vegan)
library(limma)
data<- read.csv("OTU.csv",header = T,row.names = 1)
group<-read.csv("Sample_Treatment.csv")
group <- factor(group$soil_s,levels = unique(group$soil_s))

#_______________________________________________________________________________________________design
#===========Group
design <- model.matrix(~0 + group)
rownames(design) <- colnames(data)
colnames(design) <- levels(group)
# Matrix
cont.matrix <- makeContrasts("R1 vs S1" = R1-S1, 
                             "R1 vs S2" = R1-S2, 
                             "R2 vs S1" = R2-S1,
                             "R2 vs S2" = R2-S2,
                             levels = colnames(design)
)
cont.matrix

#_______________________________________________________________________________________________difference analysis
logFCcutoff <- log2(0)

adjPvalueCutoff <- 0.05
fit <- lmFit(data, design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
plotSA(fit2)
summary(decideTests(fit2, p.value = adjPvalueCutoff,adjust.method = "none", lfc = logFCcutoff))


dt <- decideTests(fit2,p.value = adjPvalueCutoff,adjust.method = "none", lfc = logFCcutoff)
de.common <- which(dt[,1]!=0 & dt[,2]!=0 & dt[,3]!=0& dt[,4]!=0)
length(de.common)
vennDiagram(dt[,1:4], circle.col=c("#3C5488E5", "#91D1C2E5", "#4DBBD5E5","#8B8F9E"))

colnames(fit2)
R1_vs_S1 <- topTreat(fit2, coef=1, n=Inf)
R1_vs_S2 <- topTreat(fit2, coef=2, n=Inf)
R2_vs_S1 <- topTreat(fit2, coef=3, n=Inf)
R2_vs_S2 <- topTreat(fit2, coef=4, n=Inf)


#_______________________________________________________________________________________________Plot
rm(list=ls())
setwd("C:/Users/NG/Desktop/Wheat/work")
library(ggplot2)
library(tidyverse)
library(ggrepel)
data<- read.csv("G2_Fun_limma.csv",header = T)

data$label <- ifelse(data$pvalue<0.01,"pvalue<0.01","pvalue>=0.01")

top10sigR1_VS_S1 <- filter(data,cluster=="R1_VS_S1") %>% distinct(ASV,.keep_all = T) %>% top_n(10,abs(log2FC))
head(top10sigR1_VS_S1)
top10sigR1_VS_S2 <- filter(data,cluster=="R1_VS_S2") %>% distinct(ASV,.keep_all = T) %>% top_n(10,abs(log2FC))
head(top10sigR1_VS_S2)
top10sigR2_VS_S1 <- filter(data,cluster=="R2_VS_S1") %>% distinct(ASV,.keep_all = T) %>% top_n(10,abs(log2FC))
head(top10sigR2_VS_S1)
top10sigR2_VS_S2 <- filter(data,cluster=="R2_VS_S2") %>% distinct(ASV,.keep_all = T) %>% top_n(10,abs(log2FC))
head(top10sigR2_VS_S2)


top10sig <- rbind(top10sigR1_VS_S1,top10sigR1_VS_S2,top10sigR2_VS_S1,top10sigR2_VS_S2)

data$size <- case_when(!(data$ASV %in% top10sig$ASV)~ 1,data$ASV %in% top10sig$ASV ~ 2)

dt <- filter(data,size==1)
head(dt)


p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = log2FC, color = label),
              size = 0.85,
              width =0.4)
p


p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = log2FC, color = label),
              size = 1,
              width =0.4)
p


dfbar<-data.frame(x=c("R1_VS_S1","R1_VS_S2","R2_VS_S1","R2_VS_S2"),
                  y=c(8,8,8,8))
dfbar1<-data.frame(x=c("R1_VS_S1","R1_VS_S2","R2_VS_S1","R2_VS_S2"),
                   y=c(-7,-7,-7,-7))

p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1


p2 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = cluster, y = log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = log2FC, color = label),
              size = 1,
              width =0.4)
p2


dfcol<-data.frame(x=c("R1_VS_S1","R1_VS_S2","R2_VS_S1","R2_VS_S2"),
                  y=0,
                  label=c("R1_VS_S1","R1_VS_S2","R2_VS_S1","R2_VS_S2"))
mycol <- c("#E64B357F","#4DBBD57F","#00A0877F","#91D1C27F")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=1.0,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)
p3

#Cluster
p4 <- p3+
  geom_text_repel(
    data=top10sig,
    aes(x=cluster,y=log2FC,label=ASV),
    force = 1.2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last")
  )
p4


p5 <- p4 +
  scale_color_manual(name=NULL,
                     values = c("#293A6A","#CBE1E3"))
p5


p6 <- p5+
  labs(x="Cluster",y="average logFC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =4,
            color ="white")
p6


p7 <- p6+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )
p7



#_______________________________________________________________________________________________Relative abundance
spe <-read.csv("spe_RA.csv",row.names = 1,header = T)
tax <-read.csv("RA_taxonomy.csv",row.names = 1,header = T)

group <- read.csv("Sample_Treatment.csv",header = T)

rel_spe <- spe
head(rel_spe)

library("tidyr")

rel_Order <- rel_spe %>%
  group_by(Order = tax$Genus) %>% 
  summarise_all(sum)
write.csv(rel_Order,"RA_Genus.csv")


dat<-read.csv("RA_Genus.csv", header = TRUE, row.names = 1)

dat$RowSum <- rowSums(dat)
dat<- dat[order(dat$RowSum, decreasing = TRUE), ]
dat$RowSum <- NULL

##top 20
dd<-rbind(colSums(dat[22:as.numeric(length(rownames(dat))),]),dat[21:1,])
rownames(dd)[1]<-"Others2"
dd<-dd[-c(1,15),]
rel_spe<-dd
write.csv(rel_spe,"RA_Genus.csv")

library(tidyverse)
library(ggprism)
library(vegan)
library(amplicon)
dd <- read.csv('RA_Genus.csv', header = TRUE, row.names = 1)
metadata<- read.csv('Sample_Treatment.csv', header = T)

bb<-merge(t(dd),dplyr::select(metadata,sample,soil_s),
          by.x = "row.names",by.y ="sample")

kk<-tidyr::gather(bb,Order,Abundance,-c(soil_s,Row.names))

hh <- kk %>%
  group_by(soil_s,Order) %>%
  dplyr :: summarise(Abundance=mean(Abundance))

hh$Order = factor(hh$Order,order = T,levels = row.names(dd))

library(ggsci)


p1 <- ggplot(hh,aes(x = soil_s,y = Abundance,fill = Order)) + 
  geom_bar(stat = "identity",color="black",width = 0.5)
  scale_fill_manual(values = yanse) +
  scale_color_manual(values = yanse) +
  labs(x="", y="Relative abundance (%)")+
  scale_x_discrete(limits = c("R1","R2","S1","S2"))+
  guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  scale_y_continuous(expand = c(0,0))
p1


cc =dplyr::select(bb,Row.names,soil_s,everything())
dat <- gather(cc,Order,value,-(Row.names:soil_s))
p2<-ggplot(dat,aes(x = soil_s,y =value,color = soil_s)) +
  stat_boxplot(geom="errorbar",size=1.0,width=0.3,linetype="solid")+
  geom_boxplot()+
  geom_jitter(width=0.2,size=2)+
  scale_color_manual(values=c("#293A6A","#748BBA","#3C8A84","#CBE1E3","#293A6A","#8B8F9E","#8B8F9E","#344D8B","#9FA1B6","#9FA1B6","#3C8A84","#C8C8C8","#C8C8C8","#2F9E98","#2F9E98","#CBE1E3"))+
  theme_bw()+
  labs(x="", y="Relative abundance (%)")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_blank(data = dat,aes(x = soil_s,y =value))+
  facet_wrap(~Order,scales = "free_y",nrow = 4)+
  theme(panel.border=element_rect(color = "black", fill = NA, size = 1),axis.text.x=element_text(angle=45,vjust = 0.5,hjust = 0.5,colour='black',size=12))
p2

### Kruskal-Wallis
library(stats)
library(tidyverse)

boxplerk <- function(X, Y, main = NULL, xlab = NULL, ylab = NULL, 
                     bcol = c("#293A6A", "#293A6A", "#8B8F9E", "#8B8F9E"), 
                     p.adj = "none", cexy = 1, varwidth = TRUE, las = 1, paired = FALSE) {
  aa <- levels(as.factor(Y))
  an <- as.character(c(1:length(aa)))
  tt1 <- matrix(nrow = length(aa), ncol = 7)    
  for (i in 1:length(aa)) {
    temp <- X[Y == aa[i]]
    tt1[i, 1] <- mean(temp, na.rm = TRUE)
    tt1[i, 2] <- sd(temp, na.rm = TRUE) / sqrt(length(temp))
    tt1[i, 3] <- sd(temp, na.rm = TRUE)
    tt1[i, 4] <- min(temp, na.rm = TRUE)
    tt1[i, 5] <- max(temp, na.rm = TRUE)
    tt1[i, 6] <- median(temp, na.rm = TRUE)
    tt1[i, 7] <- length(temp)
  }
  tt1 <- as.data.frame(tt1)
  row.names(tt1) <- aa
  colnames(tt1) <- c("mean", "se", "sd", "min", "max", "median", "n")
  boxplot(X ~ Y, main = main, xlab = xlab, ylab = ylab, las = las, col = bcol, 
          cex.axis = cexy, cex.lab = cexy, varwidth = varwidth)    
  require(agricolae)
  Yn <- factor(Y, labels = an)
  comp <- kruskal(X, Yn, p.adj = p.adj)
  sig <- "ns"
  if (paired == TRUE & length(aa) == 2) {
    coms <- wilcox.test(X ~ Yn, paired = TRUE)
    pp <- coms$p.value
  } else {
    pp <- comp$statistics$p.chisq
  }
  if(pp <= 0.1) sig <- "."
  if(pp <= 0.05) sig <- "*"
  if(pp <= 0.01) sig <- "**"
  if(pp <= 0.001) sig <- "***"
  gror <- comp$groups[order(rownames(comp$groups)), ]
  tt1$rank <- gror$X
  tt1$group <- gror$groups
  mtext(sig, side = 3, line = 0.5, adj = 0, cex = 2, font = 1)   
  if(pp <= 0.1) {
    mtext(tt1$group, side = 3, at = c(1:length(aa)), line = 0.5, cex = 1, font = 4)
  }
  list(comparison = tt1, p.value = pp)
}


dat <- dat %>%rename(Group = soil_s,)

genera <- unique(dat$Order)


if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

results <- list()
for (Order in genera) {
  subset_dat <- dat[dat$Order == Order, ]
  png_filename <- paste(output_dir, paste0("Kruskal_Wallis_", Order, ".png"), sep = "/")
  png(png_filename, width = 800, height = 600)
  result <- boxplerk(subset_dat$value, subset_dat$Group, 
                     main = paste("Kruskal-Wallis Test for", Order), 
                     xlab = "", ylab = "Relative abundance %", 
                     bcol = c("#293A6A", "#293A6A", "#8B8F9E", "#8B8F9E"), p.adj = "fdr", las = 1)
  dev.off()
  results[[Order]] <- result
}




####PSTI LM
library(ggplot2)
df<-read.csv("PSTI.csv")

LM_results<- lm(PSTI~RDPWC,data =df)
summary(LM_results)

#se置信区间
p<-ggplot(df, aes(x=RDPWC, y=PSTI))+
  geom_point(pch=21,fill='#3C8A84',stroke = NA,size=3.5)+#点大小
  geom_smooth(method = 'lm',se = T, fill="#CBE1E3",color = '#3C8A84') +
  theme_bw()+
  theme(axis.text=element_text(colour='black',size=10))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  annotate('text', label = 'R2 = 0.160; p < 0.001', x =0.4, y = 0.6, size =4,color="#748BBA")
p




####RFM
library(tidyverse)
library(randomForest)
library(export)
data<-as.data.frame(data)
testData2 <- data %>% 
  mutate(Treatments = as.factor(Treatments))

set.seed(123)
otu_rf= randomForest(Treatments ~ ., data = testData2, importance=TRUE, proximity=TRUE)
otu_rf#1-estimate of  error rate



imp_otu <-  as_tibble(round(importance(otu_rf), 2),rownames = "OTUid") %>% 
  arrange(desc(MeanDecreaseAccuracy))
imp_otu#MeanDecreaseAccuracy


ncol(testData2)
myotu= testData2[,-1]
set.seed(123)
result<-rfcv(myotu, testData2$Treatments, cv.fold=5, scale = "log", step = 0.9)
result1<-result
result1$n.var
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))
error.cv <- data.frame(num = result$n.var, error.1 =  result$error.cv)



for (i in 124:127){
  print(i)
  set.seed(i)
  result= rfcv(myotu, testData2$Treatments, cv.fold=5, scale = "log", step = 0.9)
  error.cv = cbind(error.cv, result$error.cv)
}

n.var <- error.cv$num
error.cv <- error.cv[,-1]
colnames(error.cv)<- paste('err',1:5,sep='.')
err.mean <-apply(error.cv,1,mean)
allerr<-data.frame(num=n.var,err.mean=err.mean,error.cv) 
head(allerr[,1:6])

#plot
optimal = 621
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=0.5, colour="black"),
                   axis.line.y=element_line(size=0.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=12),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=12),
                   text=element_text(family="sans", size=12))
max(allerr$num)

ggplot() + 
  geom_line(aes(x = allerr$num, y = allerr$err.1), colour = '#CBE1E3',size=0.5) + 
  geom_line(aes(x = allerr$num, y = allerr$err.2), colour = '#CBE1E3',size=0.5) + 
  geom_line(aes(x = allerr$num, y = allerr$err.3), colour = '#CBE1E3',size=0.5) + 
  geom_line(aes(x = allerr$num, y = allerr$err.4), colour = '#CBE1E3',size=0.5) + 
  geom_line(aes(x = allerr$num, y = allerr$err.5), colour = '#CBE1E3',size=0.5) + 
  geom_line(aes(x = allerr$num, y = allerr$err.mean), colour = '#3C8A84',size=0.5) + 
  geom_vline(xintercept = optimal, colour='#CBE1E3', lwd=0.36, linetype=2) + 
  coord_trans(x = "log2") +
  scale_x_continuous(breaks = c(1, 10, 20, 100, 300,600, 1604)) + # , max(allerr$num)
  labs( x='Number of OTUs ', y='Cross-validation error rate') + 
  annotate("text", x = optimal, y = max(allerr$err.mean), label=paste("Optimal = ", optimal, sep=""))+
  main_theme ->p1
p1


#top 621 otu
write.csv(imp_otu[1:621,], "Biomarker.csv")
imp_otu[1:30,] %>%
  dplyr::select(OTUid,MeanDecreaseAccuracy) %>%
  arrange(MeanDecreaseAccuracy) %>% 
  mutate(OTUid =  forcats::fct_inorder(OTUid)) %>% 
  ggplot(aes(x = OTUid, y = MeanDecreaseAccuracy))+
  geom_bar(stat = "identity")+
  labs(x = "", y = "Mean decrease accuracy")+
  geom_col(fill = '#3C8A84') +
  coord_flip()+
  main_theme ->p2
p2


library(devtools)
library(customLayout)
library(magrittr)
lay1 = lay_new(mat = matrix(1, ncol = 1))
lay2 =lay_bind_col(lay1, lay1, widths = c(1, 1))
lay_show(lay2)
lay_grid(list(p1,p2), lay2)
