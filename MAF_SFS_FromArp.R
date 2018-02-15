#Read a modified arp  file, just the rows with genotypes

data<-read.table(file="/Users/tpyhajar/Documents/Marja/PurePops29607Loci_NewTP.arp", sep="\t")
#Remove unencessary columns, if needed
data<-data[,-1]

#Naming table
rivi<-as.character(data[1,3])
numerorivi<-as.integer(unlist(strsplit(rivi,split=" ")))
rivinr<-length(test[,1])

#make a matrix for genotype data
new<-matrix(nrow=rivinr, ncol=length(numerorivi))


#make a new table using the third column where all the genotype data is as list
for(i in 1:rivinr){
  rivi<-as.character(data[i,3])
  numerorivi<-as.integer(unlist(strsplit(rivi,split=" ")))
  new[i,]<-numerorivi
}

#check for missing data
hist(colSums(apply(new, 2, is.na))/length(new[,1]))
sum(colSums(apply(new, 2, is.na))==0)

#id data without missing data
nomissing<-colSums(apply(new, 2, is.na))==0

#remove missing data
new_nomis<-new[,nomissing]

#double check
table(apply(new_nomis, 2, is.na))

#funciton to count minor allele frequency for a single population
frequency1pop<-function(alleles){
	als<-unique(alleles)
	als<-als[!is.na(als)]
	al1<-als[1]
	al2<-als[2]
	#which is more commom
	al1sum<-sum(alleles==al1, na.rm=T)
	al2sum<-sum(alleles==al2, na.rm=T)
	if(al1sum>al2sum){al1<-al2}
	
	frekv<-sum(alleles==al1, na.rm=T)
	return(frekv)

	
}

#test the function for one column/locus
frequency1pop(new_nomis[,1])

#applky to whole table
allfreqs<-apply(new_nomis, 2, frequency1pop)

hist(allfreqs)
hist(allfreqs, breaks=200)

#Use poplabels to controls pops
poplabels<-c(9,9,9,9,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,17,17,17,17,17,17,17,17,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,15,15,15,15,16,16,16,16,18,18,18,18,18,18,18,18,18,18,19,19,19,19,3,3,3,3,3,3,3,3)

poplabels<-as.factor(poplabels)


#testing
frekvenssit(wild_nomis[,1])

#Apply
allfreqs<-apply(wild_nomis, 2, frekvenssit)

#inspecting results
hist(allfreqs)
hist(allfreqs, breaks=48)
sum(allfreqs==0)
sum(allfreqs==1)
sum(allfreqs==2)
barplot(table(allfreqs))

#without monomorphic
barplot(table(allfreqs)[-1])

#this can be provided to fastsimcoal2
table(allfreqs)



#2d SFS



#function to calculate allele frquencies for two populations
frequencies2pop<-function(alleles, n1){
  als<-unique(alleles)
  als<-als[!is.na(als)]
  al1<-als[1]
  al2<-als[2]
  #which is more commom
  al1sum<-sum(alleles==al1, na.rm=T)
  al2sum<-sum(alleles==al2, na.rm=T)
  if(al1sum>al2sum){al1<-al2}
  
  frekv1<-sum(alleles[1:npop1]==al1, na.rm=T)
  frekv2<-sum(alleles[(npop1+1):length(alleles)]==al1, na.rm=T)
  
  #returns the frequency of minor allele (counted for joint sample of 2 populations) in two pops
  return(c(frekv1,frekv2))
}

#choose 2 populations
#like this or by using poplabels
wild<-c(1:98)
dom<-(99:200)

#just wilds
pop1<-new_nomis[wild,]
pop2<-new_nomis[dom,]

npop1<-length(pop1[,1])
npop2<-length(pop2[,1])


#combine the two pops by rbind
pop2test<-rbind(pop1,pop2)

#just checking
class(pop2test[1,1])

#Testing that monomorphics are correctly labeled
for (i in 1:100) {
  print(frequencies2pop(pop2test[,i],npop1))
}
#looking good

#apply, also provide the nr of samples in the first pop
allelecounts<-apply(pop2test, 2, frequencies2pop, n1=npop1)

#make a matrix to classify allele counts
sfs2d<-matrix(nrow=npop1+1, ncol=npop2+1)
rownames(sfs2d)<-c(0:npop1)
colnames(sfs2d)<-c(0:npop2)


#stupid for loop to populate based on result table
for(i in 0:npop1){
  for (j in 0:npop2){
    sfs2d[(i+1),(j+1)]<-sum(allelecounts[1,]==i&allelecounts[2,]==j)
  }
}

#this can be provided to fastsimcoal
#note that this is bases on minor allele frequency and has to be treated as folded 2D SFS
sfs2d
heatmap(sfs2d, Rowv=NA, Colv = NA)
write.table(sfs2d,file="/Users/tpyhajar/Documents/Marja/2dsfs.txt")
