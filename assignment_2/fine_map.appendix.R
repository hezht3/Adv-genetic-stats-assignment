require(tidyverse)

setwd("C:/Users/Johnathan He/fusion_twas-master/LDREF")
data <- read.table("1000G.EUR.22.bim")

data %>% filter(V4 %in% (28374004-500000):(28374004+500000)) %>% nrow()

# Find sentinel SNP

setwd("D:/OneDrive - Johns Hopkins/Course/140.686.01 - Advanced Methods for Statistical Genetics and Genomics/assignment/homework_2")
Z_score_region <- read.table("./files/Z.txt")   # Z-score in TTC28 +/- 500kb region
colnames(Z_score_region) <- c("rsid", "Z")

Z_score_gw <- read.table("./files/data_BC17_1.txt", header = TRUE)   # Summary statistics genome wide
Z_score_region <- Z_score_region %>% left_join(Z_score_gw, by = "rsid")   # Merge summary statistics

sentinel <- Z_score_region %>% filter(P == min(P)) %>% pull(rsid)   # Sentinel SNP

# LD matrix

LD_matrix <- read.table("./files/LD.txt") %>% as.matrix()
rownames(LD_matrix) <- Z_score_region$rsid
colnames(LD_matrix) <- Z_score_region$rsid

# LD of causal variants with sentinel SNP

PIP_1 <- read.table("./Fine map output/out_causal_1_95_post", header = TRUE)
PIP_2 <- read.table("./Fine map output/out_causal_2_95_post", header = TRUE)
PIP_3 <- read.table("./Fine map output/out_causal_3_95_post", header = TRUE)

causal_1 <- PIP_1 %>% arrange(desc(`Causal_Post._Prob.`)) %>% slice(1) %>% pull(SNP_ID)
causal_2 <- PIP_2 %>% arrange(desc(`Causal_Post._Prob.`)) %>% slice(1:2) %>% pull(SNP_ID)
causal_3 <- PIP_3 %>% arrange(desc(`Causal_Post._Prob.`)) %>% slice(1:3) %>% pull(SNP_ID)

LD_matrix[sentinel, causal_1]
LD_matrix[sentinel, causal_2]
LD_matrix[sentinel, causal_3]

# Credible set

CS_95 <- list()   # 95% credible set
for(i in 1:3) {
    CS_95[[i]] <- read.table(paste0("./Fine map output/out_causal_", i, "_95_set"))
}
CS_99 <- list()   # 99% credible set
for(i in 1:3) {
    CS_99[[i]] <- read.table(paste0("./Fine map output/out_causal_", i, "_99_set"))
}

PIP_1 %>% filter(SNP_ID %in% CS_95[[1]]$V1) %>% 
    mutate(Causal_Post._Prob. = scales::scientific(Causal_Post._Prob.))
PIP_1 %>% filter(SNP_ID %in% CS_95[[2]]$V1) %>% 
    mutate(Causal_Post._Prob. = scales::scientific(Causal_Post._Prob.))
PIP_1 %>% filter(SNP_ID %in% CS_95[[3]]$V1) %>% 
    mutate(Causal_Post._Prob. = scales::scientific(Causal_Post._Prob.))

PIP_1 %>% filter(SNP_ID %in% CS_99[[1]]$V1) %>% 
    mutate(Causal_Post._Prob. = scales::scientific(Causal_Post._Prob.))
PIP_1 %>% filter(SNP_ID %in% CS_99[[2]]$V1) %>% 
    mutate(Causal_Post._Prob. = scales::scientific(Causal_Post._Prob.))
PIP_1 %>% filter(SNP_ID %in% CS_99[[3]]$V1) %>% 
    mutate(Causal_Post._Prob. = scales::scientific(Causal_Post._Prob.))

tiff("./Fine map output/Figures/LD.tiff",
     width = 4000, height = 3000, pointsize = 12, res = 600)
LD_matrix[CS_99[[3]]$V1, CS_99[[3]]$V1] %>% 
    as.data.frame() %>% 
    mutate(snp2 = rownames(LD_matrix[CS_99[[3]]$V1, CS_99[[3]]$V1])) %>% 
    pivot_longer(cols = -snp2, names_to = "snp", values_to = "ld") %>% 
    ggplot(aes(x = snp, y = snp2, fill = ld)) +
    geom_tile() +
    scale_fill_gradient(limits = c(0, 1), low = "white", high = "red") +
    xlab("") +
    ylab("") +
    theme_minimal()
dev.off()

# Plot credible set

tiff("./Fine map output/Figures/CS95_1.tiff",
     width = 4500, height = 3000, pointsize = 12, res = 600)
Z_score_region %>% 
    mutate(p_log = - log10(P)) %>% 
    mutate(cs_indi = case_when(rsid == sentinel ~ "Sentinel SNP",
                               rsid %in% CS_95[[1]]$V1 ~ "In CS",
                               TRUE ~ "Not in CS")) %>% 
    ggplot(aes(x = bp, y = p_log, color = cs_indi)) +
    geom_point() +
    geom_text(data = subset(Z_score_region %>% 
                                mutate(p_log = - log10(P)) %>% 
                                mutate(cs_indi = case_when(rsid == sentinel ~ "Sentinel SNP",
                                                           rsid %in% CS_95[[1]]$V1 ~ "In CS",
                                                           TRUE ~ "Not in CS")),
                          rsid %in% c(sentinel, CS_95[[1]]$V1)),
              aes(label = rsid), show.legend = F,
              hjust = 0, vjust = -0.5) +
    scale_color_manual("", values = c("#0099e5", "#caccd1", "#ff4c4c")) +
    xlab("Position") +
    ylab("- log10 (p-value)") +
    theme(axis.ticks = element_blank(), legend.background = element_blank(), 
          legend.key = element_blank(), panel.background = element_blank(), 
          panel.border = element_blank(), strip.background = element_blank(), 
          plot.background = element_blank(), complete = TRUE, legend.position = "bottom")
dev.off()

tiff("./Fine map output/Figures/CS95_2.tiff",
     width = 4500, height = 3000, pointsize = 12, res = 600)
Z_score_region %>% 
    mutate(p_log = - log10(P)) %>% 
    mutate(cs_indi = case_when(rsid == sentinel ~ "Sentinel SNP",
                               rsid %in% CS_95[[2]]$V1 ~ "In CS",
                               TRUE ~ "Not in CS")) %>% 
    ggplot(aes(x = bp, y = p_log, color = cs_indi)) +
    geom_point() +
    geom_text(data = subset(Z_score_region %>% 
                                mutate(p_log = - log10(P)) %>% 
                                mutate(cs_indi = case_when(rsid == sentinel ~ "Sentinel SNP",
                                                           rsid %in% CS_95[[2]]$V1 ~ "In CS",
                                                           TRUE ~ "Not in CS")),
                            rsid %in% c(sentinel, CS_95[[2]]$V1)),
              aes(label = rsid), show.legend = F,
              hjust = 0, vjust = -0.5) +
    scale_color_manual("", values = c("#0099e5", "#caccd1", "#ff4c4c")) +
    xlab("Position") +
    ylab("- log10 (p-value)") +
    theme(axis.ticks = element_blank(), legend.background = element_blank(), 
          legend.key = element_blank(), panel.background = element_blank(), 
          panel.border = element_blank(), strip.background = element_blank(), 
          plot.background = element_blank(), complete = TRUE, legend.position = "bottom")
dev.off()

tiff("./Fine map output/Figures/CS95_3.tiff",
     width = 4500, height = 3000, pointsize = 12, res = 600)
Z_score_region %>% 
    mutate(p_log = - log10(P)) %>% 
    mutate(cs_indi = case_when(rsid == sentinel ~ "Sentinel SNP",
                               rsid %in% CS_95[[3]]$V1 ~ "In CS",
                               TRUE ~ "Not in CS")) %>% 
    ggplot(aes(x = bp, y = p_log, color = cs_indi)) +
    geom_point() +
    geom_text(data = subset(Z_score_region %>% 
                                mutate(p_log = - log10(P)) %>% 
                                mutate(cs_indi = case_when(rsid == sentinel ~ "Sentinel SNP",
                                                           rsid %in% CS_95[[3]]$V1 ~ "In CS",
                                                           TRUE ~ "Not in CS")),
                            rsid %in% c(sentinel, CS_95[[3]]$V1)),
              aes(label = rsid), show.legend = F,
              hjust = 0, vjust = -0.5) +
    scale_color_manual("", values = c("#0099e5", "#caccd1", "#ff4c4c")) +
    xlab("Position") +
    ylab("- log10 (p-value)") +
    theme(axis.ticks = element_blank(), legend.background = element_blank(), 
          legend.key = element_blank(), panel.background = element_blank(), 
          panel.border = element_blank(), strip.background = element_blank(), 
          plot.background = element_blank(), complete = TRUE, legend.position = "bottom")
dev.off()

tiff("./Fine map output/Figures/CS99_1.tiff",
     width = 4500, height = 3000, pointsize = 12, res = 600)
Z_score_region %>% 
    mutate(p_log = - log10(P)) %>% 
    mutate(cs_indi = case_when(rsid == sentinel ~ "Sentinel SNP",
                               rsid %in% CS_99[[1]]$V1 ~ "In CS",
                               TRUE ~ "Not in CS")) %>% 
    ggplot(aes(x = bp, y = p_log, color = cs_indi)) +
    geom_point() +
    geom_text(data = subset(Z_score_region %>% 
                                mutate(p_log = - log10(P)) %>% 
                                mutate(cs_indi = case_when(rsid == sentinel ~ "Sentinel SNP",
                                                           rsid %in% CS_99[[1]]$V1 ~ "In CS",
                                                           TRUE ~ "Not in CS")),
                            rsid %in% c(sentinel, CS_99[[1]]$V1)),
              aes(label = rsid), show.legend = F,
              hjust = 0, vjust = -0.5) +
    scale_color_manual("", values = c("#0099e5", "#caccd1", "#ff4c4c")) +
    xlab("Position") +
    ylab("- log10 (p-value)") +
    theme(axis.ticks = element_blank(), legend.background = element_blank(), 
          legend.key = element_blank(), panel.background = element_blank(), 
          panel.border = element_blank(), strip.background = element_blank(), 
          plot.background = element_blank(), complete = TRUE, legend.position = "bottom")
dev.off()

tiff("./Fine map output/Figures/CS99_2.tiff",
     width = 4500, height = 3000, pointsize = 12, res = 600)
Z_score_region %>% 
    mutate(p_log = - log10(P)) %>% 
    mutate(cs_indi = case_when(rsid == sentinel ~ "Sentinel SNP",
                               rsid %in% CS_99[[2]]$V1 ~ "In CS",
                               TRUE ~ "Not in CS")) %>% 
    ggplot(aes(x = bp, y = p_log, color = cs_indi)) +
    geom_point() +
    geom_text(data = subset(Z_score_region %>% 
                                mutate(p_log = - log10(P)) %>% 
                                mutate(cs_indi = case_when(rsid == sentinel ~ "Sentinel SNP",
                                                           rsid %in% CS_99[[2]]$V1 ~ "In CS",
                                                           TRUE ~ "Not in CS")),
                            rsid %in% c(sentinel, CS_99[[2]]$V1)),
              aes(label = rsid), show.legend = F,
              hjust = 0, vjust = -0.5) +
    scale_color_manual("", values = c("#0099e5", "#caccd1", "#ff4c4c")) +
    xlab("Position") +
    ylab("- log10 (p-value)") +
    theme(axis.ticks = element_blank(), legend.background = element_blank(), 
          legend.key = element_blank(), panel.background = element_blank(), 
          panel.border = element_blank(), strip.background = element_blank(), 
          plot.background = element_blank(), complete = TRUE, legend.position = "bottom")
dev.off()

tiff("./Fine map output/Figures/CS99_3.tiff",
     width = 4500, height = 3000, pointsize = 12, res = 600)
Z_score_region %>% 
    mutate(p_log = - log10(P)) %>% 
    mutate(cs_indi = case_when(rsid == sentinel ~ "Sentinel SNP",
                               rsid %in% CS_99[[3]]$V1 ~ "In CS",
                               TRUE ~ "Not in CS")) %>% 
    ggplot(aes(x = bp, y = p_log, color = cs_indi)) +
    geom_point() +
    geom_text(data = subset(Z_score_region %>% 
                                mutate(p_log = - log10(P)) %>% 
                                mutate(cs_indi = case_when(rsid == sentinel ~ "Sentinel SNP",
                                                           rsid %in% CS_99[[3]]$V1 ~ "In CS",
                                                           TRUE ~ "Not in CS")),
                            rsid %in% c(sentinel, CS_99[[3]]$V1)),
              aes(label = rsid), show.legend = F,
              hjust = 0, vjust = -0.5) +
    scale_color_manual("", values = c("#0099e5", "#caccd1", "#ff4c4c")) +
    xlab("Position") +
    ylab("- log10 (p-value)") +
    theme(axis.ticks = element_blank(), legend.background = element_blank(), 
          legend.key = element_blank(), panel.background = element_blank(), 
          panel.border = element_blank(), strip.background = element_blank(), 
          plot.background = element_blank(), complete = TRUE, legend.position = "bottom")
dev.off()
