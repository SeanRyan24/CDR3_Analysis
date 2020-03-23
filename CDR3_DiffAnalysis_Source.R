rm(list = ls())
library(tidyverse)

#Load raw data and differential processed data
Vm_data_1 <- read.csv("VM1.csv")
Vm_data_2 <- read.csv("VM2.csv")
Vm_data_3 <- read.csv("VM3.csv")
Vm_data_4 <- read.csv("VM4.csv")
Vm_data_5 <- read.csv("VM5.csv")

N_data_1 <- read.csv("N1.csv")
N_data_2 <- read.csv("N2.csv")
N_data_3 <- read.csv("N3.csv")
N_data_4 <- read.csv("N4.csv")
N_data_5 <- read.csv("N5.csv")

diff_data <- read.csv("DiffAnalysis_Comb.csv")

#merging all Virtual memory TCR sequences into one file and all naive memory TCR sequences into another file
VM_data <- merge(Vm_data_1, Vm_data_2, all = TRUE)
VM_data <- merge(VM_data, Vm_data_3, all = TRUE)
VM_data <- merge(VM_data, Vm_data_4, all = TRUE)
VM_data <- merge(VM_data, Vm_data_5, all = TRUE)

N_data <- merge(N_data_1, N_data_2, all = TRUE)
N_data <- merge(N_data, N_data_3, all = TRUE)
N_data <- merge(N_data, N_data_4, all = TRUE)
N_data <- merge(N_data, N_data_5, all = TRUE)

##Adding a column that contains the P and N additions added in each CDR3 sequence
#first the Vm data:
VM_data<- VM_data %>% mutate(PN_Addn = JReadBegin - VReadEnd - 1)
#turn negative 1 into 0:
VM_data<- VM_data %>% mutate(PN_Addn = replace(PN_Addn, PN_Addn == -1, 0))
#Now the naive data:
N_data <- N_data %>% mutate(PN_Addn = JReadBegin - VReadEnd - 1)
#turn negative 1 into 0:
N_data <- N_data %>% mutate(PN_Addn = replace(PN_Addn, PN_Addn == -1, 0))
#rename copy
VM_data <- VM_data %>% rename(read_count = copy)
N_data <- N_data %>% rename(read_count = copy)
##cleaning the data for the non-weighted ONLY
#Deleting bad reads that are shown as an "*" in the data
Bad_Vm_reads <- subset(VM_data, CDR3.pep. %in% "*")
VM_data_clean <- anti_join(VM_data, Bad_Vm_reads)

#data set contains multiple entries for the same CDR3 regions these are removed by the duplicated function keepong only the one with the highest amount of reads
VM_data_clean <- VM_data_clean %>% arrange(desc(read_count)) 
VM_data_clean <- subset(VM_data_clean, !duplicated(CDR3.pep.))

#repeating these steps for the naive data set:
Bad_N_reads <- subset(N_data, CDR3.pep. %in% "*")
N_data_clean <- anti_join(N_data, Bad_N_reads)

N_data_clean <- N_data_clean %>% arrange(desc(read_count))
N_data_clean <- subset(N_data_clean, !duplicated(CDR3.pep.))

#seprarating raw Diff data into VM and N based on LogFC values
diff_data <- diff_data %>% rename(CDR = ï..CDR)
VM_diff_data <- diff_data %>% filter(logFC > 0)

N_diff_data <- diff_data %>% filter(logFC < 0)

#finding these differentially expressed TCRs in the raw data:
VM_data_clean <- VM_data_clean %>% rename(CDR = CDR3.pep.)
VM_diff_info <- inner_join(VM_data_clean, VM_diff_data, by = "CDR")

N_data_clean <- N_data_clean %>% rename(CDR = CDR3.pep.)
N_diff_info <- inner_join(N_data_clean, N_diff_data, by = "CDR")
#for unnormalized plots:
VM_Perc_indiv <- VM_diff_info %>% 
                  group_by(PN_Addn) %>% 
                  tally() %>% 
                  mutate(Percent_in_Repitore = (n/sum(n))*100)

N_Perc_indiv <- N_diff_info %>% 
                group_by(PN_Addn) %>% 
                tally() %>% 
                mutate(Percent_in_Repitore = (n/sum(n))*100)

####Doing weighting for all reads of these TCRs by using the raw uncleaned data set and allowing duplicates
VM_data <- VM_data %>% rename(CDR = CDR3.pep.)
VM_data_info_FC <- inner_join(VM_data, VM_diff_data, by = "CDR")
test_VM <- subset(VM_data_info_FC, !duplicated(CDR))

N_data <- N_data %>% rename(CDR = CDR3.pep.)
N_data_info_FC <- inner_join(N_data, N_diff_data, by = "CDR")
test_N <- subset(N_data_info_FC, !duplicated(CDR))

#For normalized plots:
VM_Perc_Norm <-  VM_data_info_FC%>% 
                  group_by(PN_Addn) %>% 
                  tally(wt = read_count) %>% 
                  mutate(Percent_in_Repitore = (n/sum(n))*100)
#now for the Naive data:
N_Perc_Norm <-  N_data_info_FC %>% 
                group_by(PN_Addn) %>% 
                tally(wt = read_count) %>% 
                mutate(Percent_in_Repitore = (n/sum(n))*100)
###Plotting
#Starting with unweighted data:
pl_VM1 <- VM_Perc_indiv %>% ggplot() + 
          aes(x = as.factor(PN_Addn), y = Percent_in_Repitore) +
          geom_histogram(stat = "identity", fill = "blue") +
          labs(y= "Percentage of TCRs", x = "Number of P and N Additions") +
          ggtitle("Virtual Memory Repitore P and N Addition Frequency")
pl_VM1

pl_N1 <- N_Perc_indiv %>% ggplot() + 
        aes(x = as.factor(PN_Addn), y = Percent_in_Repitore) +
        geom_histogram(stat = "identity", fill = "red") +
        labs(y= "Percentage of TCRs", x = "Number of P and N Additions") +
        ggtitle("Naive Repitore P and N Addition Frequency")
pl_N1

Pl_comb <- ggplot() + 
            geom_histogram(data = VM_Perc_indiv, stat = "identity", show.legend = TRUE,
                 aes(x = PN_Addn, y = Percent_in_Repitore, group = 1), color='blue', fill = "blue", alpha = 0.4) + 
            geom_histogram(data = N_Perc_indiv, stat = "identity", show.legend = TRUE,
                 aes(x = PN_Addn, y = Percent_in_Repitore, group = 1), color='red', fill = "red", alpha = 0.4) +
            scale_x_continuous(limits = c(-2, 16)) +
            labs(y= "Percentage of TCRs", x = "Number of P and N Additions") +
            ggtitle("Virtual Memory Repitore (Blue) vs Naive Repitore (Red) P and N Addition Frequency")
Pl_comb
##weighted count data:
#VM data:
pl_VM_norm <- VM_Perc_Norm %>% ggplot() + 
              aes(x = as.factor(PN_Addn), y = Percent_in_Repitore) +
              geom_histogram(stat = "identity", fill = "blue") +
              labs(y= "Percentage of Total TCR Reads", x = "Number of P and N Additions") +
              ggtitle("Virtual Memory Repitore P and N Addition Frequency per CDR3 Read Count") 
pl_VM_norm
#Plotting the Naive Data:
pl_N_norm <- N_Perc_Norm %>% ggplot() + 
              aes(x = as.factor(PN_Addn), y = Percent_in_Repitore) +
              geom_histogram(stat = "identity", fill = "red") +
              labs(y= "Percentage of Total TCR Reads", x = "Number of P and N Additions") +
              ggtitle("Naive Repitore P and N Addition Frequency per CDR3 Read Count")
pl_N_norm
#overlay of the per read data:
Pl_comb_norm <- ggplot() + 
                geom_histogram(data = VM_Perc_Norm, stat = "identity", show.legend = TRUE,
                 aes(x = PN_Addn, y = Percent_in_Repitore, group = 1), color='blue', fill = "blue", alpha = 0.4) + 
                geom_histogram(data = N_Perc_Norm, stat = "identity", show.legend = TRUE,
                 aes(x = PN_Addn, y = Percent_in_Repitore, group = 1), color='red', fill = "red", alpha = 0.4) +
                scale_x_continuous(limits = c(-2, 16)) +
                labs(y= "Percentage of Total TCR Reads", x = "Number of P and N Additions") +
                ggtitle("Virtual Memory Repitore (Blue) vs Naive Repitore (Red) P and N Addition Frequency per CDR3 Read Count")
Pl_comb_norm

#adding in statistics: 
#the paired Wilcoxon signed-rank test (two-tailed), which tests for distribution shift, 
#and the Kolmolgorov-Smirnov test, which tests for differences in distribution shape

#Adding zeros for paired test
N_indiv_stat <- N_Perc_indiv %>% add_row(PN_Addn = 13, n = 0, Percent_in_Repitore = 0) %>% arrange(PN_Addn)
VM_indiv_stat <- VM_Perc_indiv %>% 
                  add_row(PN_Addn = 12, n = 0, Percent_in_Repitore = 0) %>% 
                  add_row(PN_Addn = 14, n = 0, Percent_in_Repitore = 0) %>% 
                  arrange(PN_Addn)
VM_Norm_stat <- N_Perc_indiv %>% 
                add_row(PN_Addn = 15, n = 0, Percent_in_Repitore = 0) %>% 
                add_row(PN_Addn = 13, n = 0, Percent_in_Repitore = 0) %>% 
                arrange(PN_Addn)
#performing wilcoxon test
wilcox.test(VM_indiv_stat$Percent_in_Repitore, N_indiv_stat$Percent_in_Repitore, paired=TRUE, alternative = "less")
wilcox.test(VM_Norm_stat$Percent_in_Repitore, N_Perc_Norm$Percent_in_Repitore, paired=TRUE, alternative = "less")
#neither is significant
#preforming KS test
ks.test(VM_indiv_stat$Percent_in_Repitore, N_indiv_stat$Percent_in_Repitore, alternative = "less")
ks.test(VM_Norm_stat$Percent_in_Repitore, N_Perc_Norm$Percent_in_Repitore, alternative = "less")
#also not significant

