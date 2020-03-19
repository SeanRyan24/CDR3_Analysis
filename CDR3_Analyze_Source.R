rm(list = ls())
library(tidyverse)

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

#merging all Virtual memory TCR sequences into one file and all naive memory TCR sequences into another file
VM_data <- merge(Vm_data_1, Vm_data_2, all = TRUE)
VM_data <- merge(VM_data, Vm_data_3, all = TRUE)
VM_data <- merge(VM_data, Vm_data_4, all = TRUE)
VM_data <- merge(VM_data, Vm_data_5, all = TRUE)

N_data <- merge(N_data_1, N_data_2, all = TRUE)
N_data <- merge(N_data, N_data_3, all = TRUE)
N_data <- merge(N_data, N_data_4, all = TRUE)
N_data <- merge(N_data, N_data_5, all = TRUE)

##cleaning the data
#Deleting bad reads that are shown as an "*" in the data
Bad_Vm_reads <- subset(VM_data, CDR3.pep. %in% "*")
VM_data_clean <- anti_join(VM_data, Bad_Vm_reads)

#we are only looking for recurrent CDR3 regions so we will remove those with counts < 10, copy is a bad column name so we will rename it:
VM_data_clean <- VM_data_clean %>% rename(read_count = copy)
VM_data_clean <- subset(VM_data_clean, read_count > 10)

#data set contains multiple entries for the same CDR3 regions these are removed by the duplicated function keepong only the one with the highest amount of reads
VM_data_clean <- VM_data_clean %>% arrange(desc(read_count)) 
VM_data_clean <- subset(VM_data_clean, !duplicated(CDR3.pep.))

#repeating these steps for the naive data set:
Bad_N_reads <- subset(N_data, CDR3.pep. %in% "*")
N_data_clean <- anti_join(N_data, Bad_N_reads)

N_data_clean <- N_data_clean %>% rename(read_count = copy)
N_data_clean <- subset(N_data_clean, read_count > 10)

N_data_clean <- N_data_clean %>% arrange(desc(read_count))
N_data_clean <- subset(N_data_clean, !duplicated(CDR3.pep.))

##Adding a column that contains the P and N additions added in each CDR3 sequence
#first the Vm data:
VM_data_clean <- VM_data_clean %>% mutate(PN_Addn = JReadBegin - VReadEnd - 1)
#turn negative 1 into 0:
VM_data_clean <- VM_data_clean %>% mutate(PN_Addn = replace(PN_Addn, PN_Addn == -1, 0))
#Now the naive data:
N_data_clean <- N_data_clean %>% mutate(PN_Addn = JReadBegin - VReadEnd - 1)
#turn negative 1 into 0:
N_data_clean <- N_data_clean %>% mutate(PN_Addn = replace(PN_Addn, PN_Addn == -1, 0))

##We need to find the CDR3 sequences that exist only in VM or only in Naive
#first the VM data:
VM_not_in_N <- anti_join(VM_data_clean, N_data_clean, by = "CDR3.pep.")
#Now the Naive Data:
N_not_in_VM <- anti_join(N_data_clean, VM_data_clean, by = "CDR3.pep.")

##plot the frequency of P and N additions for each TCR repitore:
#First the VM data:
VM_Percentage <- VM_not_in_N %>% 
                  group_by(PN_Addn) %>% 
                  tally() %>% 
                  mutate(Percent_in_Repitore = (n/sum(n))*100)

pl_VM1 <- VM_Percentage %>% ggplot() + 
            aes(x = as.factor(PN_Addn), y = Percent_in_Repitore) +
            geom_histogram(stat = "identity", fill = "blue") +
            labs(y= "Percentage of TCRs", x = "Number of P and N Additions") +
            ggtitle("Virtual Memory Repitore P and N Addition Frequency")
pl_VM1
#Now thae Naive data:
N_Percentage <- N_not_in_VM %>% 
                group_by(PN_Addn) %>% 
                tally() %>% 
                mutate(Percent_in_Repitore = (n/sum(n))*100)

pl_N1 <- N_Percentage %>% ggplot() + 
          aes(x = as.factor(PN_Addn), y = Percent_in_Repitore) +
          geom_histogram(stat = "identity", fill = "red") +
          labs(y= "Percentage of TCRs", x = "Number of P and N Additions") +
          ggtitle("Naive Repitore P and N Addition Frequency")
pl_N1
#now creating an overlay of the two histograms:
Pl_comb <- ggplot() + 
            geom_histogram(data = VM_Percentage, stat = "identity", show.legend = TRUE,
                 aes(x = PN_Addn, y = Percent_in_Repitore, group = 1), color='blue', fill = "blue", alpha = 0.4) + 
            geom_histogram(data = N_Percentage, stat = "identity", show.legend = TRUE,
                 aes(x = PN_Addn, y = Percent_in_Repitore, group = 1), color='red', fill = "red", alpha = 0.4) +
            scale_x_continuous(limits = c(-2, 40)) +
            labs(y= "Percentage of TCRs", x = "Number of P and N Additions") +
            ggtitle("Virtual Memory Repitore (Blue) vs Naive Repitore (Red) P and N Addition Frequency")
Pl_comb
#Altering the data so that we are counting individual TCR reads in the population as opposed to just the CDR3 identities:
VM_Perc_Norm <-  VM_not_in_N %>% 
                group_by(PN_Addn) %>% 
                tally(wt = read_count) %>% 
                mutate(Percent_in_Repitore = (n/sum(n))*100)
#now for the Naive data:
N_Perc_Norm <-  N_not_in_VM %>% 
                group_by(PN_Addn) %>% 
                tally(wt = read_count) %>% 
                mutate(Percent_in_Repitore = (n/sum(n))*100)
#Plotting the VM data:
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
                scale_x_continuous(limits = c(-2, 40)) +
                labs(y= "Percentage of Total TCR Reads", x = "Number of P and N Additions") +
                ggtitle("Virtual Memory Repitore (Blue) vs Naive Repitore (Red) P and N Addition Frequency per CDR3 Read Count")
Pl_comb_norm
