### Script for running PGLMM models using MCMCglmm
### Made by: Alex Kurilshikov
### Date: 01 Feb 2024
### ----------------------------------------------------------------------------

## Loading libraries 
library(ape)
library(foreach)
library(stringr)
library(MCMCglmm)
library(phangorn)
library(pbapply)
library(dispRity)

setwd("~/Documents/LocalData/UMCG/LifeLines-NEXT/CS_Baby_biome/")

# models ------------------------------------------------------------------


metadata2 = read.table("METADATA_INFANTS_EARLY_CS_BABY_BIOME_09_06_2023_UPDATED_FEEDING.txt",
                       header=T,as.is = T,sep="\t")



metaphlan = read.table("TAXA_CLEAN_CS_BABY_BIOME_14_03_2023.txt",
                       header=T,as.is = T,sep="\t")

set.seed(1214)
## Making trees from scratch
tree_files = dir("updated_trees/raxml_trees/")
tree_files = paste0("updated_trees/raxml_trees/",tree_files)


trees = lapply(tree_files,read.tree)
names(trees) = sub("^updated_trees/raxml_trees/RAxML_bestTree.","",
                   sub("[.]StrainPhlAn4.tre$","",tree_files))



trees = lapply(trees,function(tree){
  labels = tree$tip.label
  labels.short = sub("_.*","",labels)
  which.duplicated = which(duplicated(labels.short))
  duplicated.labels = labels.short[which.duplicated]
  for_removing = c()
  for (i in duplicated.labels) {
    both_labels = labels[which(labels.short == i)]
    strlen = str_length(both_labels)
    for_removing = c(for_removing,both_labels[which(strlen > min(strlen))])
  }
  tree2 = drop.tip(tree,for_removing)
  tree2
}
)
trees = lapply(trees,function(x){x$tip.label = sub("_.*","",x$tip.label);x})

#MCMC model functions
execute_mcmc_feeding.standard = function(name_tree,tree_object,scale = F) {
  print(name_tree)
  tree = tree_object[[name_tree]]
  metadata.subset = metadata2[match(tree$tip.label,metadata2$NG_ID),]
  metadata.subset = metadata.subset[
    metadata.subset$Timepoint_categorical!="MOM",]
  metadata.4analysis = data.frame(sample_id = metadata.subset$NG_ID,
                                  baby_id = metadata.subset$CS_BABY_BIOME_ID,
                                  outcome = metadata.subset$feeding_mode)
  metadata.4analysis = metadata.4analysis[!is.na(metadata.4analysis$outcome),]
  metadata.4analysis$outcome = as.integer(metadata.4analysis$outcome=="breast_feeding")
  tree = drop.tip(tree,setdiff(tree$tip.label,metadata.4analysis$sample_id))
  
  
  if (sum(unlist(lapply(split(metadata.4analysis,metadata.4analysis$outcome),function(x) length(unique(x$baby_id))))>=2)>=2) {
    tree.ultra1 = tree
    
    Ainv.1<-try(inverseA(tree.ultra1,nodes = "TIPS",scale = scale)$Ainv)
    if(class(Ainv.1)[1]=="try-error") Ainv.1<-try(inverseA(remove.zero.brlen(tree.ultra1),nodes = "TIPS",scale = scale)$Ainv)

    
    prior.Fix <- list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)),R=list(V=1,fix =1))
    prior.chisq<-list(G=list(G1=list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                             G2=list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)
                             ),
                      R=list(V=1,fix=1))
    #defined as Ainv
    model.fix <- try(MCMCglmm(outcome ~ 1, 
                            random = ~sample_id + baby_id, 
                            ginverse = list(sample_id = Ainv.1),
                            prior = prior.Fix,
                            data = metadata.4analysis,
                            family = "threshold",nitt = 600000,thin = 100,burnin=10000))
    model.chisq <- try(MCMCglmm(outcome ~ 1, 
                                random = ~sample_id + baby_id, 
                                ginverse = list(sample_id = Ainv.1),
                                prior = prior.chisq,
                                data = metadata.4analysis,
                                family = "threshold",nitt = 600000,thin=100,burnin=10000))
    
    both_models = list(model.fix,model.chisq,N = nrow(metadata.4analysis),
                       Nbabies = length(unique(metadata.4analysis$baby_id)),bug = name_tree)
    if(any(c(class(Ainv.1),
             class(tree.ultra1),
             class(model.fix),
             class(model.chisq))=="try-error")) {both_models[[1]] <- NA;both_models[[2]] <-NA} 
    #defined as covariance matrix
    
    
  } else {
    both_models <- list(
      NA,NA,N = nrow(metadata.4analysis),
      Nbabies = length(unique(metadata.4analysis$baby_id)),bug = name_tree
    )
  }
  
  both_models
}

## AB models
execute_mcmc_AB.standard = function(name_tree,tree_object,scale = F) {
  print(name_tree)
  tree = tree_object[[name_tree]]
  metadata.subset = metadata2[match(tree$tip.label,metadata2$NG_ID),]
  metadata.subset = metadata.subset[
    metadata.subset$Timepoint_categorical!="MOM",]
  metadata.4analysis = data.frame(sample_id = metadata.subset$NG_ID,
                                  baby_id = metadata.subset$CS_BABY_BIOME_ID,
                                  outcome = metadata.subset$Randomization_AB_all_numeric)
  metadata.4analysis = metadata.4analysis[!is.na(metadata.4analysis$outcome),]
  tree = drop.tip(tree,setdiff(tree$tip.label,metadata.4analysis$sample_id))
  
  
  if (sum(unlist(lapply(split(metadata.4analysis,metadata.4analysis$outcome),function(x) length(unique(x$baby_id))))>=2)>=2) {
    tree.ultra1 = tree
    
    Ainv.1<-try(inverseA(tree.ultra1,nodes = "TIPS",scale = scale)$Ainv)
    if(class(Ainv.1)[1]=="try-error") Ainv.1<-try(inverseA(remove.zero.brlen(tree.ultra1),nodes = "TIPS",scale = scale)$Ainv)
    
    
    prior.Fix <- list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)),R=list(V=1,fix =1))
    prior.chisq<-list(G=list(G1=list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                             G2=list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)
    ),
    R=list(V=1,fix=1))
    #defined as Ainv
    model.fix <- try(MCMCglmm(outcome ~ 1, 
                              random = ~sample_id + baby_id, 
                              ginverse = list(sample_id = Ainv.1),
                              prior = prior.Fix,
                              data = metadata.4analysis,
                              family = "threshold",nitt = 600000,thin = 100,burnin=10000))
    model.chisq <- try(MCMCglmm(outcome ~ 1, 
                                random = ~sample_id + baby_id, 
                                ginverse = list(sample_id = Ainv.1),
                                prior = prior.chisq,
                                data = metadata.4analysis,
                                family = "threshold",nitt = 600000,thin=100,burnin=10000))
    
    both_models = list(model.fix,model.chisq,N = nrow(metadata.4analysis),
                       Nbabies = length(unique(metadata.4analysis$baby_id)),bug = name_tree)
    if(any(c(class(Ainv.1),
             class(tree.ultra1),
             class(model.fix),
             class(model.chisq))=="try-error")) {both_models[[1]] <- NA;both_models[[2]] <-NA} 
    #defined as covariance matrix
    
    
  } else {
    both_models <- list(
      NA,NA,N = nrow(metadata.4analysis),
      Nbabies = length(unique(metadata.4analysis$baby_id)),bug = name_tree
    )
  }
  
  both_models
}

mcmc_parse = function(x) {do.call(rbind, lapply(x,function(x){
  if (class(x[[1]])=="logical"|class(x)[[1]] == "try-error") {
    result = data.frame(
      bug = x[["bug"]],
      N=x[["N"]],
      Nbabies=x[["Nbabies"]],
      
     
      effectiveSize.fix = NA,
      H2mode.fix = NA,
      H2mean.fix = NA,
      H2median.fix = NA,
      post.HPD.05.fix = NA,
      post.HPD.95.fix = NA,
      effectiveSize.Chisq = NA,
      H2mode.Chisq = NA,
      H2mean.Chisq = NA,
      H2median.Chisq = NA,
      post.HPD.05.Chisq = NA,
      post.HPD.95.Chisq = NA
      )
  } else {
    post.herit1 = x[[1]]$VCV[,"sample_id"] / (rowSums(x[[1]]$VCV[,c(1,3)]))
    post.h21 = posterior.mode(post.herit1)
    
    post.int1 = HPDinterval(post.herit1, 0.95)
    
    post.herit2 = x[[2]]$VCV[,"sample_id"] / (rowSums(x[[2]]$VCV[,c(1,3)]))
    post.h22 = posterior.mode(post.herit2)
    post.int2 = HPDinterval(post.herit2, 0.95)
    
   
    
    result = data.frame(
      bug = x[["bug"]],
      N = x[["N"]],
      Nbabies = x[["Nbabies"]],
      
     
      effectiveSize.fix = effectiveSize(post.herit1)[1],
      H2mode.fix = post.h21,
      H2mean.fix = mean(post.herit1),
      H2median.fix = median(post.herit1),
      post.HPD.05.fix = post.int1[1],
      post.HPD.95.fix = post.int1[2],
      effectiveSize.Chisq = effectiveSize(post.herit2)[1],
      H2mode.Chisq = post.h22,
      H2mean.Chisq = mean(post.herit2),
      H2median.Chisq = median(post.herit2),
      post.HPD.05.Chisq = post.int2[1],
      post.HPD.95.Chisq = post.int2[2])
  }
  result
}))
}


# Pictures ----------------------------------------------------------------

make_tree.AB = function(x,tree_object,  chr = F){
  metadata.subset = metadata2[match(tree_object[[x]]$tip.label,metadata2$NG_ID),]
  subtree = tree_object[[x]]
  metadata.subset$Randomization_AB_all_numeric[metadata.subset$Randomization_AB_all_numeric ==0] = "-"
  metadata.subset$Randomization_AB_all_numeric[metadata.subset$Randomization_AB_all_numeric ==1] = "+"
  
  subtree$tip.label = paste0("AB",metadata.subset$Randomization_AB_all_numeric,
                             "_",
                             metadata.subset$CS_BABY_BIOME_ID,
                             "_",
                             metadata.subset$Timepoint_categorical
                             
  )
  subtree = drop.tip(subtree,"ABNA_NA_NA")
  if(chr == T) subtree = chronos(subtree,lambda = 1)
  
  plot(subtree,
       cex.main = 0.9,
       main = sub(".t__.*","",sub(".*s__","",grep(x,colnames(metaphlan),value = T)),),
       tip.color = as.integer(as.factor(sub("^[^_]*_","",sub("_[^_]*$","",subtree$tip.label)))),
       font = 1+as.integer(grepl("[+]",subtree$tip.label))+2,
       cex = 0.7 - as.integer(grepl("[-]",subtree$tip.label))*0.1)
  table(metadata.subset$Randomization_AB_all_numeric,metadata.subset$CS_BABY_BIOME_ID)
}


make_tree.feed = function(x,tree_object,chr = F){
  metadata.subset = metadata2[match(tree_object[[x]]$tip.label,metadata2$NG_ID),]
  metadata.subset$BF = c("MF/FF","BF")[as.integer(metadata.subset$feeding_mode == "breast_feeding")+1]
  subtree = tree_object[[x]]
  
  subtree$tip.label = paste0(metadata.subset$BF,
                             "_",
                             metadata.subset$CS_BABY_BIOME_ID,
                             "_",
                             metadata.subset$Timepoint_categorical
                             
  )
  subtree = drop.tip(subtree,"NA_NA_NA")
  if(chr == T) subtree = chronos(subtree,lambda = 1)
  
  plot(subtree,
       cex.main = 0.9,
       main = sub(".t__.*","",sub(".*s__","",grep(x,colnames(metaphlan),value = T)),),
       tip.color = as.integer(as.factor(sub("^[^_]*_","",sub("_[^_]*$","",subtree$tip.label)))),
       font = 1+as.integer(grepl("[+]",subtree$tip.label))+2,
       cex = 0.7 - as.integer(grepl("[-]",subtree$tip.label))*0.1)
  table(metadata.subset$BF,metadata.subset$CS_BABY_BIOME_ID)
}



# phangorn trees ----------------------------------------------------------
trees_good1 = names(trees[AB.parsed[!is.na(AB.parsed$H2.fix),1]])
trees_good2 = trees[feeding.parsed[!is.na(feeding.parsed$H2.fix),1]]
trees_good = unique(c(names(trees_good1),names(trees_good2)))


aln_files = dir("updated_trees/aln/")
aln_files = paste0("updated_trees/aln/",aln_files)
alns = lapply(aln_files,\(x) {read.phyDat(x,format = "fasta")})
names(alns) = sub("updated_trees/aln/","",sub(".StrainPhlAn4_concatenated.aln","",aln_files))
alns = alns[names(alns) %in% trees_good]

phangorn_trees = pblapply(alns,\(x) {
 aln.tmp = x
 names(aln.tmp) = sub("_.*","",names(aln.tmp))
 aln.tmp1 = aln.tmp[names(aln.tmp) %in% metadata2$NG_ID]
 d1 = pml_bb(aln.tmp1,model = "GTR",method = "ultrametric",
             rearrangement="NNI", control = pml.control(trace = 0))
 d1$tree
})
phangorn.feeding.threshold = lapply(names(phangorn_trees), \(x) execute_mcmc_feeding.standard(x,phangorn_trees,scale = T))
phangorn.AB.threshold = lapply(names(phangorn_trees), \(x) execute_mcmc_AB.standard(x,phangorn_trees,scale = T))

phangorn.feed.parsed.threshold = mcmc_parse(phangorn.feeding.threshold)
phangorn.AB.parsed.threshold = mcmc_parse(phangorn.AB.threshold)

source("trees.R")
pdf(onefile = F,width =4,height = 7)
make_tree.AB("t__SGB6952",phangorn_trees)
make_tree.AB("t__SGB17234",phangorn_trees)
make_tree.AB("t__SGB7962",phangorn_trees)
make_tree.feed("t__SGB6191",phangorn_trees)
make_tree.feed("t__SGB4037",phangorn_trees)
dev.off()


