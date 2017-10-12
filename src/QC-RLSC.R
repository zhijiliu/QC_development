rm(list=ls())
gc()
setwd("/Users/zhijiliu/Google Drive/Mac Rstudio/Pathway analysis/PE/QC_DEV")
create.project("../QC_DEV", 
               merge.strategy = "allow.non.conflict")

# source('http://bioconductor.org/biocLite.R')
# biocLite("statTarget")
# install.packages("ProjectTemplate")
library("ProjectTemplate")
library(statTarget)
library(xlsx)

a1=read.xlsx(file='PE plan 6-27-17.xlsx',sheetIndex = 1)
a1=a1[c(2:6,10:14),2:9]

x=c()
for(i in 1:10){
  x=c(x,as.character(unlist(a1[i,])))
}
x=x[!is.na(x)]
x=x[-grep('QCB',x)]
x=x[-1]

a2=read.xlsx(file='PE plan 6-27-17.xlsx',sheetIndex = 2)
a2=a2[c(2:6,10:14,17:21,24:28),2:9]

y=c()
for(i in 1:20){
  y=c(y,as.character(unlist(a2[i,])))
  
}

y=y[!is.na(y)]
y=y[-1]

ph=data.frame(sample=c(x,y),batch=c(rep(1,length(x)),rep(2,length(y))),
              class=rep(1,length(x)+length(y)),order=1:(length(x)+length(y)))

ph$class[grep('QC',ph$sample)]='NA'
ph$class[grep('P',ph$sample)]=2

write.csv(ph,file='ph.csv',row.names = F)


load('/Users/dna/Dropbox/ling.manuscripts/PE/code/input.matrix.promedex.rdata')

pd=conc_PE_Promedex.t
pd=pd[-grep('QCB',pd$sample_ID),]
matrix.pd=pd[,4:2249]
matrix.pd=t(matrix.pd)
colnames(matrix.pd)=pd$sample_ID
matrix.pd=data.frame(matrix.pd)
matrix.pd$name=row.names(matrix.pd)
matrix.pd=matrix.pd[,c(69,1:68)]

load('/Users/dna/Dropbox/ling.manuscripts/PE/code/input.matrix.rdata')
st=conc_PE_stanford.t
matrix.st=st[,12:2257]
matrix.st=t(matrix.st)
colnames(matrix.st)=st$sample_ID
matrix.st=data.frame(matrix.st)
matrix.st$name=row.names(matrix.st)
matrix.st=matrix.st[,c(155,1:154)]

r.order=row.names(matrix.st)
matrix.pd=matrix.pd[r.order,]

matrix.all=cbind(matrix.pd,matrix.st)
matrix.all=matrix.all[,-70]

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
matrix.all$name=trim(matrix.all$name)
matrix.all$name=gsub('/','-to-',matrix.all$name)

write.csv(matrix.all,file='pr.csv',row.names = F)

shiftCor('ph.csv','pr.csv')

#shiftCor('PhenoFile.csv','ProfileFile.csv')
rm(list=ls())
gc()

setwd("/Users/dna/Dropbox/ling.manuscripts/PE/code/code.new")

load('/Users/dna/Dropbox/ling.manuscripts/PE/code/input.matrix.promedex.rdata')
load('/Users/dna/Dropbox/ling.manuscripts/PE/code/input.matrix.rdata')
matrix.all.cr=read.csv('statTarget/shiftCor/After_shiftCor/shift_sample_cor.csv')

st.cr=conc_PE_stanford.t[,1:11]
st.cr=merge(st.cr,matrix.all.cr[,c(-1,-3)],by.x='sample_ID',by.y='sample',all.x=T)
st.cr=st.cr[1:138,]
conc_PE_stanford.t.cr=st.cr
conc_PE_stanford.t.cr=conc_PE_stanford.t.cr[order(conc_PE_stanford.t.cr$time,conc_PE_stanford.t.cr$subject),]

control=conc_PE_stanford.t.cr[conc_PE_stanford.t.cr$Class==0,]
control.matrix=control[,12:1560]
row.names(control.matrix)=control$ID

save(conc_PE_stanford.t,conc_PE_stanford.t.cr,control.matrix,control,file='input.matrix.st.rdata')

pd.cr=conc_PE_Promedex.t[1:64,1:3]
pd.cr=merge(pd.cr,matrix.all.cr[,c(-1,-3)],by.x='sample_ID',by.y='sample',all.x=T)
conc_PE_Promedex.t.cr=pd.cr
conc_PE_Promedex.t.cr=conc_PE_Promedex.t.cr[order(conc_PE_Promedex.t.cr$sample_ID),]

control=conc_PE_Promedex.t.cr[conc_PE_Promedex.t.cr$Class==0,]
control.matrix=control[,4:ncol(control)]
row.names(control.matrix)=control$sample_ID

save(conc_PE_Promedex.t,conc_PE_Promedex.t.cr,control.matrix,control,file='input.matrix.pd.rdata')

