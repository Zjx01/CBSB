#single folding mutation

#loading necessary packages
library(ggplot2)
library(rootSolve)

wt<-function(exp_mut,deltadelta_Gf){#exp_mut: mutation affect the expression level
  pr1=1.1*exp_mut # protein expression level
  deltag1_wt = -2 # WT folding energy, kcal/mol
  R= 1.98*10^(-3) # gas constant
  Temp= 310.15 # abs temperature for 37 degrees, Kelvin
  kwt= exp(-deltag1_wt/(R*Temp))
  
  deltaG1= deltag1_wt + deltadelta_Gf
  k=exp(-deltaG1/(R*Temp))
  
  dslnex<-function(x){
    y<-numeric(2)
    y[1]<-x[1]+x[2]-pr1
    y[2]<-k*x[2]-x[1] #to notice 不要用除式
    y
  }
  xstart<-c(0.9,0.1)
  components= multiroot(dslnex, start=xstart, positive = TRUE,rtol=1e-6,atol= 1e-10, ctol= 1e-10)$root
  return(c(components[1]/pr1,components[1],components[2]))
}

#返回的数量级有点问题,我觉得这很有问题

?multiroot()
#Simulation for the gradual mutation 
folding_fraction<-c()
for (deltadelta_Gf in seq(-2,5.5,0.1)) {
    folding_fraction<-c(folding_fraction,wt(1,deltadelta_Gf)[1])
}
deltadelta_Gf<-c(seq(-2,5.5,0.1))
structure<-data.frame(deltadelta_Gf,folding_fraction)
colnames(structure)<-c("deltadelta_Gf","folding_fraction") 
#ΔΔG f(kcal/mol),folding_fraction (%)

ggplot(data=structure,aes(x=deltadelta_Gf,y=folding_fraction*100))+
  geom_line(col="darkred")+theme(panel.border = element_rect(colour = "black",
        fill=NA),panel.background = element_blank())+geom_vline(aes(xintercept=0),linetype="dashed")+
  labs(x="ΔΔG f(kcal/mol)",y="Pf fraction (%)")


#single expression mutation, we set the deltadeltaG_f=0, to represent the mutation does not 
#interfere with the folding energy mutation-------------------
library(ggplot2)
sexp_folding_fraction<-c()
folded_amount<-c()
unfolded_amount<-c()

for (exp in seq(0.1,10,0.1)) {
    sexp_folding_fraction<-c(sexp_folding_fraction,wt(exp,0)[1])
    folded_amount<-c(folded_amount,wt(exp,0)[2])
    unfolded_amount<-c(unfolded_amount,wt(exp,0)[3])
}
expdata=data.frame(c(seq(0.1,10,0.1)*1.1),sexp_folding_fraction,folded_amount,unfolded_amount)
colnames(expdata)<-c("total_protein_amount","Pf_fraction","Pf_amount","Puf_amount")
head(expdata)


ep_1<-ggplot(data = expdata,aes(x=total_protein_amount,y=Pf_fraction))+
  geom_line()+theme(panel.border = element_rect(colour = "black",
                fill=NA),panel.background = element_blank())+labs(x="total protein amount",y="Pf fraction")+
            ylim(0.96,1)

#加上legend!!!!
ep_2<-ggplot(data = expdata,aes(x=total_protein_amount))+
  geom_line(aes(y=Pf_amount),col='red')+
  geom_line(aes(y=Puf_amount),col='blue')+
  theme(panel.border = element_rect(colour = "black",
    fill=NA),panel.background = element_blank())+
  labs(x="total protein amount",y="protein amount")+
  theme(legend.position = "right",legend.title = element_text(colour="blue", size=10, 
                                                              face="bold"))
?theme()  
#legend( "right", c("Pf amount","Puf amount"),col = c("red","blue"))
 
###为什么phenotype没有改变却致病？推测：
ep_3<-ggplot(data = expdata,aes(x=total_protein_amount,y=Pf_amount/Puf_amount))+
  geom_line(col="orange")+ylim(25.8,26)+
  theme(panel.border = element_rect(colour = "black",
      fill=NA),panel.background = element_blank())+labs(x="total protein amount",y="Pf amount/ Puf amount")
  
ggarrange(ep_2,ep_3,ep_1,ncol=3)




#expression level mutation combination with folding energy mutation-----------------
#here we gradually change the expression level to the impact of expression on phenotype
ex_folding_fraction_0.8<-c()
ex_folding_fraction_1<-c()
ex_folding_fraction_1.2<-c()
for (deltadelta_Gf in seq(-2,2,0.01)) {
  ex_folding_fraction_0.8<-c(ex_folding_fraction_0.8,wt(0.8,deltadelta_Gf))
  ex_folding_fraction_1<-c(ex_folding_fraction_1,wt(1,deltadelta_Gf))
  ex_folding_fraction_1.2<-c(ex_folding_fraction_1.2,wt(1.2,deltadelta_Gf))
}

expmu_data<-data.frame(c(seq(-2,2,0.01)),ex_folding_fraction_0.8,ex_folding_fraction_1,ex_folding_fraction_1.2)
colnames(expmu_data)[1]<-"deltadelta_Gf"
colnames(expmu_data)[2]<-"ex_folding_fraction_0.8"
colnames(expmu_data)[3]<-"ex_folding_fraction_1"
colnames(expmu_data)[4]<-"ex_folding_fraction_1.2"
dim(expmu_data)

ggplot(data=expmu_data,aes(x=deltadelta_Gf))+
  geom_line(aes(y=ex_folding_fraction_0.8),col="orange")+
  geom_line(aes(y=ex_folding_fraction_1),col="green")+
  geom_line(aes(y=ex_folding_fraction_1.2),col="blue")+
  theme(panel.border = element_rect(colour = "black",fill=NA),
        panel.background = element_blank())+
  labs(x="ΔΔG f(kcal/mol)",y="folding_fraction (%)")






#binding and folding combination----------------------------
comb_bf<- function(a1, deltadelta_G1,  deltadelta_GB1) { # a1 can be seen as a mutation that affects expression of the protein
  a1 = 1
  pr1=1.1 # protein expression level
  deltag1_wt = -2 # WT folding energy, kcal/mol
  R= 1.98*10^(-3) # gas constant
  Temp= 310.15 # abs temperature for 37 degrees, Kelvin
  deltaG1= deltag1_wt + deltadelta_G1  # Mut folding energy
  
  k4 = exp(-deltaG1/(R*Temp))  # Folding kinetic
  
  pr2=1.1  # ligand concentration, this can be a DNA molecule or small molecules as well. 
  deltaGB_ligand_wt = -5 # WT binding energy
  deltaGB_ligand1 = deltaGB_ligand_wt + deltadelta_GB1 # Mut binding energy
  k1= exp(-deltaGB_ligand1/(R*Temp))
  
  require(rootSolve)
  dslnex <- function(x) {
    y <- numeric(4) # x[1]= bound and fold, x[2]= folded,x[3]= unbound ligand(or DNA), x[4] = unfolded. 
    y[1] <- x[1]+ x[2] + x[4]- pr1*a1  # equation for total protein amount
    y[2] <- x[2]*x[3]*k1-x[1] # equation for binding-unbinding kinetics
    y[3] <- x[4]*k4-x[2] # equation for folding-unfolding kinetics
    y[4] <- x[3]+x[1]-pr2 # equation for total ligand (DNA) concentration 
    y
  } 
  xstart <- c(0.9, 0.05, 0.05, 0.01)  # give a starting position to search for the solution
  components= multiroot( dslnex, start=xstart, positive = TRUE,rtol=1e-6,atol= 1e-10, ctol= 1e-10)$root
  #print(deltadelta_G1,deltadelta_GB1)
  return(components[1]/pr1)
}



#we did not involve expression level here,so the expression level a is set to 1
data=data.frame(rep(0,76))

i=0
while(i<=76){
  for (deltadelta_G1 in seq(-2,5.5,0.1)) {
    i=i+1
    j=0
    for(deltadelta_GB1 in seq(-2,5.5,0.1)){
      j=j+1
      temp= comb_bf(1,deltadelta_G1,deltadelta_GB1)
      data[i,j]=temp
    }
  }
}
head(data)
dim(data)

#data = as.data.frame(t(data.frame(rep(0,76))))
xx = seq(-2,5.5,0.1)
yy = seq(-2,5.5,0.1)
# for (i in 1:length(xx)) {
#   tempList = c()
#   for (j in 1:length(yy)) {
#     temp = comb_bf(1,xx[i],yy[j])
#     tempList = c(tempList, temp)
#   }
#   data[i,]= tempList
# }

rownames(data)<-seq(-2,5.5,0.1)
colnames(data)<-seq(-2,5.5,0.1)

# mutation combination和 phenotype组合在一起
matrix_combfb<-matrix(ncol = 3)
for (i in 1:length(xx)) {
  loc_x <- xx[i]
  for (j in 1:length(yy)) {
    loc_y <- yy[j]
    fb_phenotype <- data[i,j]
    matrix_combfb<-rbind(matrix_combfb,c(loc_x,loc_y,fb_phenotype))
  }
}
dim(matrix_combfb)
colnames(matrix_combfb)<-c("ΔΔGf","ΔΔGb","folded and bound protein fraction")
write.csv(matrix_combfb,'matrix.csv')
getwd()
#plotting,画一下folded and bounded states 随着 folding 和binding mutation的变化

#我们可以考虑一下具体的实例，例如血红蛋白的结合
library(scatterplot3d)

p1<-scatterplot3d(x=c(matrix_combfb[1:5776,2]),y=c(matrix_combfb[1:5776,1]),z=c(matrix_combfb[1:5776,3]*100),
              main = "Folded and Bounded Protein Phenotype Distribution",
              size = 1,
              xlab = "ΔΔGb",
              ylab = "ΔΔGf",
              xlim = c(-2,5.5),
              ylim = c(-2,5.5),
              zlab ="P_fb fraction",              
              type = 'p',
              highlight.3d = T,   # choose a colorscale
              #wt：c(0,0,0.9835605), 
              opacity=0.8)

p2<-scatterplot3d(x=c(matrix_combfb[1:5776,1]),y=c(matrix_combfb[1:5776,2]),z=c(matrix_combfb[1:5776,3]*100),
              main = "Folded and Bounded Protein Phenotype Distribution",
              size = 1,
              xlab = "ΔΔGb",
              ylab = "ΔΔGf",
              xlim = c(-2,5.5),
              ylim = c(-2,5.5),
              zlab ="P_fb fraction",              
              type = 'p',
              highlight.3d = T,   # choose a colorscale
              #wt：c(0,0,0.9835605), 怎么标注上wt
              opacity=0.8)

library(ggpubr)  
ggarrange(p1,p2,ncol=2,nrow = 1)


#我们可以看出在 folded and binded protein 对于相同能量大小的binding mutation
#更加的sensitive




#在binding和folding突变蛋白的基础上，我们给予相同的能量改变ΔΔGf和ΔΔGb，观察性状的绝对值改变大小
#null hypotheis: there is no difference on phenotype when given the same alteration of ΔΔGf and ΔΔGb


com_b=c()
com_f=c()
for(deltadelta_GB1 in seq(-2,5.5,0.1)){
  com_b=c(com_b,comb_bf(1,0,deltadelta_GB1))
}

for (deltadelta_G1 in seq(-2,5.5,0.1)) {
  com_f=c(com_f,comb_bf(1,deltadelta_G1,0))
}
length(com_f)

###legend！！！！
phenotype_fb=data.frame(ΔΔG=seq(-2,5.5,0.1),ΔΔGf_p=com_f,ΔΔGb_p=com_b)
ggplot(data=phenotype_fb,aes(x=ΔΔG))+
  geom_line(aes(y=ΔΔGf_p),col='blue')+
  geom_line(aes(y=ΔΔGb_p),col='red')+
  theme(panel.border = element_rect(colour = "black",fill=NA),
        panel.background = element_blank())+
  labs(x="ΔΔG (kcal/mol)",y="P_fb fraction")+
  geom_vline(aes(xintercept=0),linetype="dashed")
 

#statistical analysis 
# probability_com=c()
# for (i in 1:100) {
#   #randomly set 100 points each time 
#   G1=sample(seq(-1,4.5,0.1),100,replace = T)
#   GB=sample(seq(-1,4.5,0.1),100,replace= T)
#   #ΔΔGf=G1, ΔΔGb=GB, wt_phenotype=comb_bf(1,0,0)=0.9835605
#   energy=sample(seq(-1,1,0.1),1)
#   count=0
#   for(j in 1:100){
#     current_phenotype=comb_bf(1,G1[j],GB[j])
#     folding_change=abs(comb_bf(1,G1[j]+energy,GB)-current_phenotype)
#     binding_change=abs(comb_bf(1,G1,GB[j]+energy)-current_phenotype)
#     if (binding_change>folding_change) {
#       count=count+1
#     }           
#   }
#   #the probability of P_GB>P_Gf
#   p=count/100
#   probability_com=c(probability_com,p)
# }
# sum(probability_com>0.5)
# 
# 
# 
# 


#Codes for Fig4 visualized in python
# X = np.arange(-2,5.6,0.1)       
# Y = np.arange(-2,5.6,0.1)                      
# for i in range(0,len(df),76): 
#   Z.append(list(df["folded and bound protein fraction"])[i: i+ 76]) 
# Z = np.array(Z) 
# plt.contourf(X,Y,Z,levels=10)