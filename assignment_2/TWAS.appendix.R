require(tidyverse)

# Merge "typical analysis and output" results - GCTA reference
setwd("C:/Users/Johnathan He/fusion_twas-master/OUT/GCTA")

top.snp <- list()
for(i in 1:22) {
    top.snp[[i]] <- read.table(paste("BC17.", i, ".top", sep = ""), header = TRUE)
}
top.snp.merge <- top.snp[[1]]
for(i in 2:22) {
    if(nrow(top.snp[[i]]) > 0) {
        top.snp.merge <- top.snp.merge %>% bind_rows(top.snp[[i]])
    }
}

top.snp.merge %>% 
    mutate(TWAS.indi = TWAS.Z > 0) %>% 
    pull(TWAS.indi) %>% 
    table()

# Merge "joint/conditional tests and plots" results - GCTA reference
joint.snp <- list()
for(i in 1:22) {
    if(paste("BC17.", i, ".top.analysis.joint_included.dat", sep = "") %in% list.files()) {
        joint.snp[[i]] <- read.table(paste("BC17.", i, ".top.analysis.joint_included.dat", sep = ""), header = TRUE)
    }
}
joint.snp.merge <- joint.snp[[1]]
for(i in 2:22) {
    if(paste("BC17.", i, ".top.analysis.joint_included.dat", sep = "") %in% list.files()) {
        joint.snp.merge <- joint.snp.merge %>% bind_rows(joint.snp[[i]])
    }
}

joint.snp.merge <- joint.snp.merge %>% filter(JOINT.P < 2.5e-6)

joint.snp.merge %>% 
    mutate(joint.indi = JOINT.Z > 0) %>% 
    pull(joint.indi) %>% 
    table()

top.snp.merge %>% 
    select(CHR, ID, TWAS.Z, TWAS.P, EQTL.ID, EQTL.Z, EQTL.GWAS.Z) %>% 
    left_join(joint.snp.merge %>% select(ID, JOINT.Z, JOINT.P), by = "ID") %>% 
    write_csv(file = "D:/OneDrive - Johns Hopkins/Course/140.686.01 - Advanced Methods for Statistical Genetics and Genomics/assignment/homework_2/TWAS output/GCTA reference.csv")

# Merge "typical analysis and output" results - GTEx reference
setwd("C:/Users/Johnathan He/fusion_twas-master/OUT/GTEx")

top.snp <- list()
for(i in 1:22) {
    top.snp[[i]] <- read.table(paste("BC17.", i, ".top", sep = ""), header = TRUE)
}
top.snp.merge <- top.snp[[1]]
for(i in 2:22) {
    if(nrow(top.snp[[i]]) > 0) {
        top.snp.merge <- top.snp.merge %>% bind_rows(top.snp[[i]])
    }
}

top.snp.merge %>% 
    mutate(TWAS.indi = TWAS.Z > 0) %>% 
    pull(TWAS.indi) %>% 
    table()

# Merge "joint/conditional tests and plots" results - GTEx reference
joint.snp <- list()
for(i in 1:22) {
    if(paste("BC17.", i, ".top.analysis.joint_included.dat", sep = "") %in% list.files()) {
        joint.snp[[i]] <- read.table(paste("BC17.", i, ".top.analysis.joint_included.dat", sep = ""), header = TRUE)
    }
}
joint.snp.merge <- joint.snp[[1]]
for(i in 2:22) {
    if(paste("BC17.", i, ".top.analysis.joint_included.dat", sep = "") %in% list.files()) {
        joint.snp.merge <- joint.snp.merge %>% bind_rows(joint.snp[[i]])
    }
}

joint.snp.merge <- joint.snp.merge %>% filter(JOINT.P < 2.5e-6)

joint.snp.merge %>% 
    mutate(joint.indi = JOINT.Z > 0) %>% 
    pull(joint.indi) %>% 
    table()

top.snp.merge %>% 
    select(CHR, ID, TWAS.Z, TWAS.P, EQTL.ID, EQTL.Z, EQTL.GWAS.Z) %>% 
    left_join(joint.snp.merge %>% select(ID, JOINT.Z, JOINT.P), by = "ID") %>% 
    write_csv(file = "D:/OneDrive - Johns Hopkins/Course/140.686.01 - Advanced Methods for Statistical Genetics and Genomics/assignment/homework_2/TWAS output/GTEx reference.csv")
