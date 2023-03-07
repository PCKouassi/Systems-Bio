#install and load required packages
#install.packages("dplyr")
library(dplyr)

#Read csv file into R and name as noise.data
#
noise.data <- read.csv(file="Documents/Grad 2021-2022/IBIS 404/HW 3/noise dataset wt and mutant.csv", header=T)
#convert to dplyr readable format i.e. table and dataframe
noise.data <- tibble::as_tibble(noise.data)

#calculate total protein fluorescence per nucleus and add as a column
noise.data <- mutate(noise.data, total.protein = GFP + mCherry)

#subset by genotype
wt <- collect(filter(noise.data, Genotype == "wt"))
mut <- collect(filter(noise.data, Genotype == "mutant"))

#density plot of protein fluorescence
plot(density(wt$total.protein), lwd=3,col="blue", main = "Density Plot of Fluorescence", xlab = "Fluorescence Intensity")
lines(density(mut$total.protein), lwd=3,col="orange")
legend("topright", col = c("blue", "orange"), lwd = c(3), box.lwd = c(1), text.col = c("blue","orange"), legend=c("Wildtype","Mutant"))

#nuclear volume comparision
plot(ecdf(wt$Volume),lwd=3,col="blue", main = "Cumulative Frequency Distribution of Nuclear Volume", cex.main = 0.9, xlab = "Nuclear Volume", ylab = "Cumulative Frequency")
lines(ecdf(mut$Volume),lwd=3,col="orange")
legend("bottomright", col = c("blue", "orange"), lwd = c(3), box.lwd = c(1), text.col = c("blue","orange"), legend=c("Wildtype","Mutant"))

#scatter plot of GFP vs mCherry
plot(mut$GFP, mut$mCherry,pch=20,col="orange", xlab ="GFP", ylab="mCherry", main = "GFP vs mCherry Scatterplot")
points(wt$GFP, wt$mCherry, pch=20, col="blue")
legend("topleft", col = c("blue", "orange"), pch=20, text.col = c("blue","orange"), legend=c("Wildtype","Mutant"))

#Code to plot the regression line
#mut$GFP, mut$mCherry

Fit<-lm(noise.data$GFP ~ noise.data$mCherry)
abline(Fit)

## rounded coefficients for better output
cf <- round(coef(Fit),4)

## sign check to avoid having plus followed by minus for negative coefficients
eq1 <- paste0("mCherry = ", cf[1],
             ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " GFP ")

## printing of the equation
mtext(eq1,3, line=-2)

#sort cells by ascending total protein level
wt <- arrange(wt, total.protein)
mut <- arrange(mut, total.protein)

#check number of datapoints/rows
nrow(wt)
nrow(mut)

#bin cells into groups of 100 cells each
N = 100
binning.vector <- seq(1,16501,N) 

#create empty lists to store cells by bin 
wt.bin = mut.bin = list(NULL)

#bin cells(rows) 1-100 in bin 1, cells(rows) 101-200 in bin 2 etc
for(i in 1:(length(binning.vector)-1)){
  wt.bin[[i]] <- wt[binning.vector[i]:binning.vector[i+1]-1,]  
  mut.bin[[i]] <- mut[binning.vector[i]:binning.vector[i+1]-1,] 
}

int_noise <- function (data, indices){
  d <- data[indices]
  d$numerator <- (d$GFP - d$mCherry)^2
  a <- mean(d$numerator)
  b <- mean(d$GFP)
  c <- mean(d$mCherry)
  product <- 2*b*c
  int_noise.sq = a/product
  return(int_noise.sq)
}

ext_noise <- function (data, indices){
  d <- data[indices]
  d$numerator1 <- mean(d$GFP*d$mCherry)-mean(d$GFP)*mean(d$mCherry)
  a1 <- mean(d$numerator1)
  b1 <- mean(d$GFP)
  c1 <- mean(d$mCherry)
  product1 <- b1*c1
  ext_noise.sq = a1/product1
  return(ext_noise.sq)
}

#create vectors to store the group-wise population characteristics mean protein and intrinsic noise level
wt.mean.prot = mut.mean.prot = c(NULL)
wt.int.noise.sq = mut.int.noise.sq = c(NULL)
wt.ext.noise.sq = mut.ext.noise.sq = c(NULL)

#calculate mean protein conc, and intrinsic noise squared for each bin of 100 cells
for(i in 1:(length(binning.vector)-1)){
  wt.int.noise.sq[i] = int_noise(wt.bin[[i]])
  wt.mean.prot[i] = mean(wt.bin[[i]]$total.protein)
  
  mut.int.noise.sq[i] = int_noise(mut.bin[[i]])
  mut.mean.prot[i] = mean(mut.bin[[i]]$total.protein)
}

#calculate mean protein conc, and extrinsic noise squared for each bin of 100 cells
for(i in 1:(length(binning.vector)-1)){
  wt.ext.noise.sq[i] = -ext_noise(wt.bin[[i]])
  #wt.mean.prot[i] = mean(wt.bin[[i]]$total.protein)
  
  mut.ext.noise.sq[i] = -ext_noise(mut.bin[[i]])
  #mut.mean.prot[i] = mean(mut.bin[[i]]$total.protein)
}

  
#plot intrinsic noise vs mean protein conc
plot(wt.mean.prot, sqrt(wt.int.noise.sq),pch=20,col="blue",xlab="Protein Concentration", ylab = "Intrinsic Noise")
points(mut.mean.prot, sqrt(mut.int.noise.sq), pch=20, col = "orange")
legend("topright", col = c("blue", "orange"), pch = 20, text.col = c("blue","orange"), legend=c("Wildtype","Mutant"))

#plot extrinsic noise vs mean protein conc
plot(wt.mean.prot, sqrt(wt.ext.noise.sq),pch=20,col="blue",xlab="Protein Concentration", ylab = "Extrinsic Noise")
points(mut.mean.prot, sqrt(mut.ext.noise.sq), pch=20, col = "orange")
legend("topright", col = c("blue", "orange"), pch = 20, text.col = c("blue","orange"), legend=c("Wildtype","Mutant"))


#Calculate inverse protein protein level 
inv_mut_prot <- 1/mut.mean.prot
inv_wt_prot<- 1/wt.mean.prot


#plot square of intrinsic noise vs mean protein conc
plot(inv_wt_prot, (wt.int.noise.sq)^2,pch=20,col="blue",xlab="Inverse Protein Concentration", ylab = "Square of Intrinsic Noise")
points(inv_mut_prot, (mut.int.noise.sq)^2, pch=20, col = "orange")
legend("topright", col = c("blue", "orange"), pch = 20, text.col = c("blue","orange"), legend=c("Wildtype","Mutant"))




noise.data$TotProtPerNuc <- (noise.data$GFP + noise.data$mCherry)*noise.data$Volume
d <- density(noise.data$TotProtPerNuc)

wt2 <- subset(noise.data, Genotype == "wt")
mut2 <- subset(noise.data, Genotype == "mutant")

wt_densities<-density(wt2$TotProtPerNuc)
mut_densities<-density(mut2$TotProtPerNuc)

ks.test(wt2$TotProtPerNuc, mut2$TotProtPerNuc)



plot(wt_densities, pch=20,col="blue", main = "Kernal Density")
points(mut_densities, pch=20, col = "orange")
points(wt_densities,pch=20,col="blue")
legend("topright", col = c("blue", "orange"), pch = 20, text.col = c("blue","orange"), legend=c("Wildtype","Mutant"))


