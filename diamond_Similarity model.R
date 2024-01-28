################################################################################
## The main functions used in cratons similarity measurement
library(readxl)
library(tidyverse)
library(corrplot)
library(matrixcalc)
library(fmsb)
library(ggplot2)
library(ggrepel)

### A function to normalize a vector in to interval [0,1] 
norm01 = function(x){
  return((x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T)))
}

### A function to calculate the mean and covariance structure of two cratons
## Input: two cratons' dataframe
## Output: the mean and correlation matrix of each craton
total_Matrix = function(df1,df2){
  df11 = apply(df1, 2, norm01)
  df22 = apply(df2, 2, norm01)
  mudf1 = colMeans(df11,na.rm = T)
  mudf2 = colMeans(df22,na.rm = T)
  
  cdf1 = cor(df1,use = "pairwise.complete.obs")
  cdf2 = cor(df2,use = "pairwise.complete.obs")

  return(list(mudf1 = mudf1,
              mudf2 = mudf2,
              cdf1=cdf1,
              cdf2=cdf2))
}

### A function to calculate the similarity measurement between two cratons 
## Input: 
## tM: the output from function total_Matrix
## p: the order of vector norm
## q: the order of matrix norm
## w1,w2: the weights vectors of each craton for elements
calculateM = function(tM,p=1,q='F',w1=NULL,w2=NULL){
  mudf1 = tM$mudf1
  mudf2 = tM$mudf2
  cdf1 = tM$cdf1
  cdf2 = tM$cdf2
  
  if(is.null(w1)){w1=rep(1,length(mudf1))}
  if(is.null(w2)){w2=rep(1,length(mudf2))}
  
  cdf1 = w1%*%t(w1)*cdf1
  cdf2 = w2%*%t(w2)*cdf2
  
  mudf1[is.na(mudf1)] = mean(na.omit(mudf1))
  mudf2[is.na(mudf2)] = mean(na.omit(mudf2))
  cdf1[is.na(cdf1)] = mean(na.omit(as.numeric(cdf1)))
  cdf2[is.na(cdf2)] = mean(na.omit(as.numeric(cdf2)))
  
  mid1 = mudf1-mudf2
  dmu_craton = sum(abs(mid1)^p)^(1/p)
  mid2 = cdf1-cdf2
  dinner_craton = norm(mid2,type = q)
  
  cor1 = 1 - dmu_craton/(sum(abs(mudf1)^p)^(1/p)+sum(abs(mudf2)^p)^(1/p))
  cor2 = 1 - dinner_craton/(norm(cdf1,type = q)+norm(cdf2,type = q))
  
  return(list(dmu_craton = dmu_craton,
              dinner_craton = dinner_craton,
              cor = c(cor1,cor2),
              sim = sqrt((cor1^2+cor2^2)/2)))
}

### The total function 
## Input the dataframe of cratons, and the output path
## Output correlation matrixs and graphs
diamond_project = function(diamond,outpath){
  craton1 = names(table(diamond$craton))
  cra = substr(craton1,1,7)
  ncra = length(cra)
  X_cra = list()
  for (i in 1:ncra) {
    X_cra[[cra[i]]] = diamond %>% filter(craton==craton1[i])
    X_cra[[cra[i]]] = X_cra[[cra[i]]][,-1]
  }
  
  X_tM = list()
  X_tM1 = list()
  for (i in 1:(ncra-1)) {
    for (j in 1:(ncra-i)) {
      tm_name = paste(cra[i],cra[i+j],sep = '')
      tm_name1 = paste(cra[i+j],cra[i],sep = '')
      X_tM[[tm_name]] = total_Matrix(X_cra[[cra[i]]],X_cra[[cra[i+j]]])
      X_tM1[[tm_name1]] = total_Matrix(X_cra[[cra[i+j]]],X_cra[[cra[i]]])
    }
  }
  
  X_W = list()
  for (i in 1:ncra) {
    X_W[[cra[i]]] = colSums(!is.na(X_cra[[cra[i]]]))/nrow(X_cra[[cra[i]]])
  }
  
  X_cor = list()
  for (i in names(X_tM)) {
    X_cor[[i]] = calculateM(X_tM[[i]],p=1,q='2',w1=X_W[[substr(i,1,3)]],w2=X_W[[substr(i,4,6)]])
  }
  X_cor1 = list()
  for (i in names(X_tM1)) {
    X_cor1[[i]] = calculateM(X_tM1[[i]],p=1,q='2',w1=X_W[[substr(i,1,3)]],w2=X_W[[substr(i,4,6)]])
  }
  X_cor2 = c(X_cor,X_cor1)
  
  gtitle = c('center','inner','total')
  for (k in 1:2) {
    cormat = diag(ncra)
    l=1
    for (j in 1:(ncra-1)) {
      for (i in (j+1):ncra) {
        cormat[i,j] = X_cor[[l]]$cor[k]
        l = l+1
      }
    }
    colnames(cormat) = cra
    rownames(cormat) = cra
    corrplot(cormat,method = 'number',type = 'lower',title = gtitle[k], number.digits = 3)
  }
  
  cormat = diag(ncra)
  l=1
  for (j in 1:(ncra-1)) {
    for (i in (j+1):ncra) {
      cormat[i,j] = X_cor[[l]]$sim
      l = l+1
    }
  }
  colnames(cormat) = cra
  rownames(cormat) = cra
  corrplot(cormat,method = 'number',type = 'lower',title = gtitle[3], number.digits = 3)
  
  X_tri = list()
  for (i in 1:ncra) {
    cra_diff = setdiff(cra,cra[i])
    X_tri[[cra[i]]] = data.frame(center = c(1,0),inner = c(1,0),outer = c(1,0))
    for (j in 1:(ncra-1)) {
      mc = paste(cra[i],cra_diff[j],sep = '')
      X_tri[[cra[i]]] = rbind(X_tri[[cra[i]]],X_cor2[[mc]]$cor)
    }
    rownames(X_tri[[cra[i]]]) = c('max','min',cra_diff)
  }
  capture.output(X_tri,file=outpath)
  
  for (i in 1:ncra) {
    triX = X_tri[[i]]
    mintri = min(triX[3:nrow(triX),1:2])  
    mintri = ((mintri*10)%/%1)/10
    
    G1 = ggplot(triX[-c(2),], aes(x = center, y = inner )) +
      geom_point(color = "firebrick") +
      labs(x = "center_distance", y = "inner_structure") +
      geom_text_repel(aes(center, inner, label = c(cra[i],rownames(triX[-c(1,2),])))) +
      ylim(c(mintri, 1))+
      xlim(c(mintri,1))+
      coord_fixed(ratio = 1)+
      ggtitle(paste("similarity_with_",craton1[i],sep = ''))+
      geom_vline(xintercept = c(1))+
      geom_hline(yintercept = c(1))+
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold.italic"))
    print(G1)
  }
}

################################################################################
## An example for using function diamond_project
# read a dataframe
diamond <- read_excel("E:/paper4_diamond/lithosphere/LDSD.xlsx")
# The data is seem as:
head(diamond)
# indicate an output path
# outpath = '/diamond_out.txt'
diamond_project(diamond,outpath)
