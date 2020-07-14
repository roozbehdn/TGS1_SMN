library(data.table)
library(scales)
library(ggplot2)
library(ggrepel)
library(pheatmap)

library(scales)
library(ggrepel)

##### Input files ###############
new_SUPPA_WT_vs_SMN_dpsi_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\SUPPA_files\\SUPPA_DE_files\\union_scallop_grazia_diff_WT_vs_SMN.dpsi"
new_SUPPA_WT_vs_SMN_psivec_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\SUPPA_files\\SUPPA_DE_files\\union_scallop_grazia_diff_WT_vs_SMN.psivec"
new_SUPPA_WT_vs_TGS1_dpsi_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\SUPPA_files\\SUPPA_DE_files\\union_scallop_grazia_diff_WT_vs_TGS.dpsi"
new_SUPPA_WT_vs_TGS1_psivec_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\SUPPA_files\\SUPPA_DE_files\\union_scallop_grazia_diff_WT_vs_TGS.psivec"
SUPPA_WT_vs_SMN_dpsi_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\SUPPA_files\\SUPPA_DE_files\\AS_diff_WT_vs_SMN.dpsi"
SUPPA_WT_vs_SMN_psivec_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\SUPPA_files\\SUPPA_DE_files\\AS_diff_WT_vs_SMN.psivec"
SUPPA_WT_vs_TGS1_dpsi_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\SUPPA_files\\SUPPA_DE_files\\AS_diff_WT_vs_TGS.dpsi"
SUPPA_WT_vs_TGS1_psivec_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\SUPPA_files\\SUPPA_DE_files\\AS_diff_WT_vs_TGS.psivec"
new_SUPPA_AS_events_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\SUPPA_files\\SUPPA_AS_files\\union_scallop_grazia_uniq_sorted_SUPPA_events.ioe"
SUPPA_AS_events_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\SUPPA_files\\SUPPA_AS_files\\SUPPA_hg38_events_strict.ioe"
AS_diff_WT_vs_SMN_avglogtpm_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\SUPPA_files\\SUPPA_DE_files\\AS_diff_WT_vs_SMN_avglogtpm.tab"
AS_diff_WT_vs_TGS1_avglogtpm_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\SUPPA_files\\SUPPA_DE_files\\AS_diff_WT_vs_TGS_avglogtpm.tab"
DE_transcripts_TGS1_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\outputs\\transcripts_DE_analysis_WT_vs_TGS1.csv"
DE_transcripts_SMN_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\outputs\\transcripts_DE_analysis_WT_vs_SMN.csv"
hg38_name_id_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\utility_files\\hg38_gene_name_ids.txt"
sample_list_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\utility_files\\sample_list.txt"
#############################################


######### read in input files ###############
new_SUPPA_WT_vs_SMN_dpsi = fread(new_SUPPA_WT_vs_SMN_dpsi_file,sep="\t",header=TRUE)
new_SUPPA_WT_vs_SMN_psivec = fread(new_SUPPA_WT_vs_SMN_psivec_file,sep="\t",header=TRUE)
new_SUPPA_WT_vs_TGS1_dpsi = fread(new_SUPPA_WT_vs_TGS1_dpsi_file,sep="\t",header=TRUE)
new_SUPPA_WT_vs_TGS1_psivec = fread(new_SUPPA_WT_vs_TGS1_psivec_file,sep="\t",header=TRUE)
SUPPA_WT_vs_SMN_dpsi = fread(SUPPA_WT_vs_SMN_dpsi_file,sep="\t",header=TRUE)
SUPPA_WT_vs_SMN_psivec = fread(SUPPA_WT_vs_SMN_psivec_file,sep="\t",header=TRUE)
SUPPA_WT_vs_TGS1_dpsi = fread(SUPPA_WT_vs_TGS1_dpsi_file,sep="\t",header=TRUE)
SUPPA_WT_vs_TGS1_psivec = fread(SUPPA_WT_vs_TGS1_psivec_file,sep="\t",header=TRUE)
new_SUPPA_AS_events = fread(new_SUPPA_AS_events_file,sep="\t",header=TRUE)
SUPPA_AS_events = fread(SUPPA_AS_events_file,sep="\t",header=TRUE)
AS_diff_WT_vs_SMN_avglogtpm = fread(AS_diff_WT_vs_SMN_avglogtpm_file,sep="\t",header=FALSE)
AS_diff_WT_vs_TGS1_avglogtpm = fread(AS_diff_WT_vs_TGS1_avglogtpm_file,sep="\t",header=FALSE)
DE_transcripts_TGS1 = fread(DE_transcripts_TGS1_file,sep=",",header = TRUE,skip = 8)
DE_transcripts_SMN = fread(DE_transcripts_SMN_file,sep=",",header = TRUE,skip = 8)
hg38_name_id = fread(hg38_name_id_file,sep="\t",header=TRUE)
sample_list = fread(sample_list_file,sep="\t",header=TRUE)
#############################################

names(new_SUPPA_WT_vs_SMN_dpsi) = c("Event","dPSI","pval")
names(new_SUPPA_WT_vs_TGS1_dpsi) = c("Event","dPSI","pval")
names(SUPPA_WT_vs_SMN_dpsi) = c("Event","dPSI","pval")
names(SUPPA_WT_vs_TGS1_dpsi) = c("Event","dPSI","pval")
names(new_SUPPA_WT_vs_SMN_psivec) = c("Event","WT_1","WT_2","SMN_1","SMN_2")
names(new_SUPPA_WT_vs_TGS1_psivec) = c("Event","WT_1","WT_2","TGS1_1","TGS1_2")
names(SUPPA_WT_vs_SMN_psivec) = c("Event","WT_1","WT_2","SMN_1","SMN_2")
names(SUPPA_WT_vs_TGS1_psivec) = c("Event","WT_1","WT_2","TGS1_1","TGS1_2")


## find the position of events 
SUPPA_WT_vs_TGS1_psivec[,event_pos:=strsplit(Event,split=";")[[1]][2],by=Event]
SUPPA_WT_vs_SMN_psivec[,event_pos:=strsplit(Event,split=";")[[1]][2],by=Event]
SUPPA_WT_vs_TGS1_dpsi[,event_pos:=strsplit(Event,split=";")[[1]][2],by=Event]
SUPPA_WT_vs_SMN_dpsi[,event_pos:=strsplit(Event,split=";")[[1]][2],by=Event]
new_SUPPA_WT_vs_TGS1_psivec[,event_pos:=strsplit(Event,split=";")[[1]][2],by=Event]
new_SUPPA_WT_vs_SMN_psivec[,event_pos:=strsplit(Event,split=";")[[1]][2],by=Event]
new_SUPPA_WT_vs_TGS1_dpsi[,event_pos:=strsplit(Event,split=";")[[1]][2],by=Event]
new_SUPPA_WT_vs_SMN_dpsi[,event_pos:=strsplit(Event,split=";")[[1]][2],by=Event]
SUPPA_AS_events[,event_pos:=strsplit(event_id,split=";")[[1]][2],by=event_id]
new_SUPPA_AS_events[,event_pos:=strsplit(event_id,split=";")[[1]][2],by=event_id]


new_SUPPA_WT_vs_TGS1_psivec = new_SUPPA_WT_vs_TGS1_psivec[!event_pos%in%SUPPA_WT_vs_TGS1_psivec$event_pos]
new_SUPPA_WT_vs_SMN_psivec = new_SUPPA_WT_vs_SMN_psivec[!event_pos%in%SUPPA_WT_vs_SMN_psivec$event_pos]
new_SUPPA_WT_vs_TGS1_dpsi = new_SUPPA_WT_vs_TGS1_dpsi[!event_pos%in%SUPPA_WT_vs_TGS1_dpsi$event_pos]
new_SUPPA_WT_vs_SMN_dpsi = new_SUPPA_WT_vs_SMN_dpsi[!event_pos%in%SUPPA_WT_vs_SMN_dpsi$event_pos]
new_SUPPA_AS_events = new_SUPPA_AS_events[!event_pos%in%SUPPA_AS_events$event_pos]
new_SUPPA_WT_vs_TGS1_psivec[,novel:=1]
new_SUPPA_WT_vs_SMN_psivec[,novel:=1]
new_SUPPA_WT_vs_TGS1_dpsi[,novel:=1]
new_SUPPA_WT_vs_SMN_dpsi[,novel:=1]
SUPPA_WT_vs_TGS1_psivec[,novel:=0]
SUPPA_WT_vs_SMN_psivec[,novel:=0]
SUPPA_WT_vs_TGS1_dpsi[,novel:=0]
SUPPA_WT_vs_SMN_dpsi[,novel:=0]

SUPPA_WT_vs_TGS1_psivec = rbind(SUPPA_WT_vs_TGS1_psivec,new_SUPPA_WT_vs_TGS1_psivec)
SUPPA_WT_vs_SMN_psivec = rbind(SUPPA_WT_vs_SMN_psivec,new_SUPPA_WT_vs_SMN_psivec)
SUPPA_WT_vs_TGS1_dpsi = rbind(SUPPA_WT_vs_TGS1_dpsi,new_SUPPA_WT_vs_TGS1_dpsi)
SUPPA_WT_vs_SMN_dpsi = rbind(SUPPA_WT_vs_SMN_dpsi,new_SUPPA_WT_vs_SMN_dpsi)
SUPPA_AS_events = rbind(new_SUPPA_AS_events,SUPPA_AS_events)
SUPPA_WT_vs_TGS1_psivec = SUPPA_WT_vs_TGS1_psivec[!duplicated(event_pos)]
SUPPA_WT_vs_SMN_psivec = SUPPA_WT_vs_SMN_psivec[!duplicated(event_pos)]
SUPPA_WT_vs_TGS1_dpsi = SUPPA_WT_vs_TGS1_dpsi[!duplicated(event_pos)]
SUPPA_WT_vs_SMN_dpsi = SUPPA_WT_vs_SMN_dpsi[!duplicated(event_pos)]
SUPPA_AS_events = SUPPA_AS_events[!duplicated(paste(event_id))]

SUPPA_WT_vs_SMN_dpsi[,ENSG:=strsplit(strsplit(Event,split=";")[[1]][1],split = ".",fixed=TRUE)[[1]][1],by=Event]
SUPPA_WT_vs_SMN_dpsi = merge(SUPPA_WT_vs_SMN_dpsi,hg38_name_id,all.x = TRUE,all.y = FALSE,by.x="ENSG",by.y="gene_id")
names(SUPPA_WT_vs_SMN_dpsi) = c("ENSG","Event","dPSI","pval","event_pos","novel","gene_name")

SUPPA_WT_vs_TGS1_dpsi[,ENSG:=strsplit(strsplit(Event,split=";")[[1]][1],split = ".",fixed=TRUE)[[1]][1],by=Event]
SUPPA_WT_vs_TGS1_dpsi = merge(SUPPA_WT_vs_TGS1_dpsi,hg38_name_id,all.x = TRUE,all.y = FALSE,by.x="ENSG",by.y="gene_id")
names(SUPPA_WT_vs_TGS1_dpsi) = c("ENSG","Event","dPSI","pval","event_pos","novel","gene_name")

## dpsi>0 means that the event is higher expressed in the KO vs WT

##############################################
## first we get pie charts for WT vs SMN ######
A3_SMN_dpsi = SUPPA_WT_vs_SMN_dpsi[Event%like%"A3"]
A5_SMN_dpsi = SUPPA_WT_vs_SMN_dpsi[Event%like%"A5"]
AF_SMN_dpsi = SUPPA_WT_vs_SMN_dpsi[Event%like%"AF"]
AL_SMN_dpsi = SUPPA_WT_vs_SMN_dpsi[Event%like%"AL"]
MX_SMN_dpsi = SUPPA_WT_vs_SMN_dpsi[Event%like%"MX"]
RI_SMN_dpsi = SUPPA_WT_vs_SMN_dpsi[Event%like%"RI"]
SE_SMN_dpsi = SUPPA_WT_vs_SMN_dpsi[Event%like%"SE"]

par(mar=c(8,1,8,5))
slices <- c(nrow(A3_SMN_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(A3_SMN_dpsi[pval>=0.05]), nrow(A3_SMN_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Alternative 3'SS for SMN vs WT",col=c("cornflowerblue","blueviolet","deeppink"))
legend("topright", c("Significantly Increased", "Unchanged", "Significantly decreased"), cex = 0.8,fill = c("cornflowerblue","blueviolet","deeppink" ))

par(mar=c(8,1,8,5))
slices <- c(nrow(A5_SMN_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(A5_SMN_dpsi[pval>=0.05]), nrow(A5_SMN_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Alternative 5'SS for SMN vs WT",col=c("cornflowerblue","blueviolet","deeppink"))
legend("topright", c("Significantly Increased", "Unchanged", "Significantly decreased"), cex = 0.8,fill = c("cornflowerblue","blueviolet","deeppink" ))

par(mar=c(8,1,8,5))
slices <- c(nrow(AF_SMN_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(AF_SMN_dpsi[pval>=0.05]), nrow(AF_SMN_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Alternative first exon for SMN vs WT",col=c("cornflowerblue","blueviolet","deeppink"))
legend("topright", c("Significantly Increased", "Unchanged", "Significantly decreased"), cex = 0.8,fill = c("cornflowerblue","blueviolet","deeppink" ))

par(mar=c(8,1,8,5))
slices <- c(nrow(AL_SMN_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(AL_SMN_dpsi[pval>=0.05]), nrow(AL_SMN_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Alternative last exon for SMN vs WT",col=c("cornflowerblue","blueviolet","deeppink"))
legend("topright", c("Significantly Increased", "Unchanged", "Significantly decreased"), cex = 0.8,fill = c("cornflowerblue","blueviolet","deeppink" ))

par(mar=c(8,1,8,5))
slices <- c(nrow(MX_SMN_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(MX_SMN_dpsi[pval >=0.05]), nrow(MX_SMN_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Mutually exclusive exons for SMN vs WT",col=c("cornflowerblue","blueviolet","deeppink"))
legend("topright", c("Significantly Increased", "Unchanged", "Significantly decreased"), cex = 0.8,fill = c("cornflowerblue","blueviolet","deeppink" ))

par(mar=c(8,1,8,5))
slices <- c(nrow(RI_SMN_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(RI_SMN_dpsi[pval >=0.05]), nrow(RI_SMN_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Intron Retention for SMN vs WT",col=c("cornflowerblue","blueviolet","deeppink"))
legend("topright", c("Significantly Increased", "Unchanged", "Significantly decreased"), cex = 0.8,fill = c("cornflowerblue","blueviolet","deeppink" ))

par(mar=c(8,1,8,5))
slices <- c(nrow(SE_SMN_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(SE_SMN_dpsi[pval >=0.05]), nrow(SE_SMN_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Exon skipping for SMN vs WT",col=c("cornflowerblue","blueviolet","deeppink"))
legend("topright", c("Significantly Increased", "Unchanged", "Significantly decreased"), cex = 0.8,fill = c("cornflowerblue","blueviolet","deeppink" ))

######################################################
########## Now we do the same for TGS1  vs WT ########
A3_TGS1_dpsi = SUPPA_WT_vs_TGS1_dpsi[Event%like%"A3"]
A5_TGS1_dpsi = SUPPA_WT_vs_TGS1_dpsi[Event%like%"A5"]
AF_TGS1_dpsi = SUPPA_WT_vs_TGS1_dpsi[Event%like%"AF"]
AL_TGS1_dpsi = SUPPA_WT_vs_TGS1_dpsi[Event%like%"AL"]
MX_TGS1_dpsi = SUPPA_WT_vs_TGS1_dpsi[Event%like%"MX"]
RI_TGS1_dpsi = SUPPA_WT_vs_TGS1_dpsi[Event%like%"RI"]
SE_TGS1_dpsi = SUPPA_WT_vs_TGS1_dpsi[Event%like%"SE"]

par(mar=c(8,1,8,5))
slices <- c(nrow(A3_TGS1_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(A3_TGS1_dpsi[pval>=0.05]), nrow(A3_TGS1_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Alternative 3'SS for TGS1 vs WT",col=c("cornflowerblue","blueviolet","deeppink"))
legend("topright", c("Significantly Increased", "Unchanged", "Significantly decreased"), cex = 0.8,fill = c("cornflowerblue","blueviolet","deeppink" ))

par(mar=c(8,1,8,5))
slices <- c(nrow(A5_TGS1_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(A5_TGS1_dpsi[pval>=0.05]), nrow(A5_TGS1_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Alternative 5'SS for TGS1 vs WT",col=c("cornflowerblue","blueviolet","deeppink"))
legend("topright", c("Significantly Increased", "Unchanged", "Significantly decreased"), cex = 0.8,fill = c("cornflowerblue","blueviolet","deeppink" ))

par(mar=c(8,1,8,5))
slices <- c(nrow(AF_TGS1_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(AF_TGS1_dpsi[pval>=0.05]), nrow(AF_TGS1_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Alternative first exon for TGS1 vs WT",col=c("cornflowerblue","blueviolet","deeppink"))

par(mar=c(8,1,8,5))
slices <- c(nrow(AL_TGS1_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(AL_TGS1_dpsi[pval>=0.05]), nrow(AL_TGS1_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Alternative last exon for TGS1 vs WT",col=c("cornflowerblue","blueviolet","deeppink"))
legend("topright", c("Significantly Increased", "Unchanged", "Significantly decreased"), cex = 0.8,fill = c("cornflowerblue","blueviolet","deeppink" ))

par(mar=c(8,1,8,5))
slices <- c(nrow(MX_TGS1_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(MX_TGS1_dpsi[pval >=0.05]), nrow(MX_TGS1_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Mutually exclusive exons for TGS1 vs WT",col=c("cornflowerblue","blueviolet","deeppink"))
legend("topright", c("Significantly Increased", "Unchanged", "Significantly decreased"), cex = 0.8,fill = c("cornflowerblue","blueviolet","deeppink" ))

par(mar=c(8,1,8,5))
slices <- c(nrow(RI_TGS1_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(RI_TGS1_dpsi[pval >=0.05]), nrow(RI_TGS1_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Intron Retention for TGS1 vs WT",col=c("cornflowerblue","blueviolet","deeppink"))
legend("topright", c("Significantly Increased", "Unchanged", "Significantly decreased"), cex = 0.8,fill = c("cornflowerblue","blueviolet","deeppink" ))

par(mar=c(8,1,8,5))
slices <- c(nrow(SE_TGS1_dpsi[(pval < 0.05) & (dPSI > 0.1)]),nrow(SE_TGS1_dpsi[pval >=0.05]), nrow(SE_TGS1_dpsi[(pval < 0.05) & (dPSI < -0.1)]))
piepercent = round(100*slices/sum(slices), 1)
lbls <- c("Significantly Increased", "Unchanged", "Significantly decreased")
pie(slices, labels =  paste(piepercent,"%",sep=""), main="Exon skipping for TGS1 vs WT",col=c("cornflowerblue","blueviolet","deeppink"))
legend("topright", c("Significantly Increased", "Unchanged", "Significantly decreased"), cex = 0.8,fill = c("cornflowerblue","blueviolet","deeppink" ))


##########################################################################################
##### plot a 3way heatmap for all significant RI events in TGS1 and SMN samples  #########
##########################################################################################

RI_SMN_psivec = SUPPA_WT_vs_SMN_psivec[Event %like% "RI"]
RI_TGS1_psivec = SUPPA_WT_vs_TGS1_psivec[Event %like% "RI"]

six_way_matrix = merge(RI_SMN_psivec,RI_SMN_dpsi[,list(Event,dPSI,pval)],all.x=TRUE,all.y=TRUE,by.x="Event",by.y="Event")
six_way_matrix = six_way_matrix[!duplicated(Event)]
setnames(six_way_matrix,old=c("dPSI","pval"),new=c("dPSI_SMN","pval_SMN"))
six_way_matrix = merge(six_way_matrix,RI_TGS1_psivec[,list(Event,TGS1_1,TGS1_2)],all.x=TRUE,all.y=TRUE,by.x="Event",by.y="Event")
six_way_matrix = six_way_matrix[!duplicated(Event)]
six_way_matrix = merge(six_way_matrix ,RI_TGS1_dpsi[,list(Event,dPSI,pval)],all.x=TRUE,all.y=TRUE,by.x="Event",by.y="Event")
six_way_matrix = six_way_matrix[!duplicated(Event)]
setnames(six_way_matrix,old=c("dPSI","pval"),new=c("dPSI_TGS1","pval_TGS1"))
six_way_matrix_sig = six_way_matrix[ (pval_SMN < 0.05 & abs(dPSI_SMN)>0.1 ) | (pval_TGS1 < 0.05 & abs(dPSI_TGS1)>0.1)]
six_way_matrix_sig[,IR_status_SMN:="Not changed"]
six_way_matrix_sig[pval_SMN < 0.05 & dPSI_SMN > 0.1,IR_status_SMN:="Sig. increase"]
six_way_matrix_sig[pval_SMN < 0.05 & dPSI_SMN < -0.1,IR_status_SMN:="Sig. decrease"]
six_way_matrix_sig[,IR_status_TGS1:="Not changed"]
six_way_matrix_sig[pval_TGS1 < 0.05 & dPSI_TGS1 > 0.1,IR_status_TGS1:="Sig. increase"]
six_way_matrix_sig[pval_TGS1 < 0.05 & dPSI_TGS1 < -0.1,IR_status_TGS1:="Sig. decrease"]

annotation_row = data.table(
  IR_status_SMN = six_way_matrix_sig$IR_status_SMN,
  IR_status_TGS1 = six_way_matrix_sig$IR_status_TGS1
)
annotation_row$IR_status_SMN = as.factor(annotation_row$IR_status_SMN)
annotation_row$IR_status_TGS1 = as.factor(annotation_row$IR_status_TGS1)
annotation_row = data.frame(annotation_row)
rownames(annotation_row) = six_way_matrix_sig$Event
six_way_matrix_sig_mat = as.matrix(six_way_matrix_sig[,list(WT_1,WT_2,SMN_1,SMN_2,TGS1_1,TGS1_2)])
colnames(six_way_matrix_sig_mat)=c("1A-WT","1B-WT","4A-SMN KO6","5A-SMN KO26","2A-TGS1 KO4","3A-TGS1 KO50")
rownames(six_way_matrix_sig_mat) = six_way_matrix_sig$Event
ann_colors = list( IR_status_SMN  = c(`Not changed` = "lightyellow2", `Sig. increase` = "firebrick2", `Sig. decrease` = "darkgreen"),IR_status_TGS1  = c(`Not changed` = "lightyellow2", `Sig. increase` = "firebrick2", `Sig. decrease` = "darkgreen"))
pheatmap(six_way_matrix_sig_mat,show_rownames = FALSE,annotation_row = annotation_row, annotation_colors = ann_colors)
##########################################################
##########################################################

## the psi is computed for the black region (psi = black / (black + gray))

##################################################################
##### plot volcanot plots for dTPM vs dPSI for all IR events ######
##################################################################

## to get the TPM values for the transcripts involved in the intron retention


RI_TGS1_dpsi[,c("V5","V6","V7","V8","int_pos","int_pos1","int_pos2","chr","intron","gene_name","ENSG"):=NULL]
RI_TGS1_dpsi = merge(RI_TGS1_dpsi,SUPPA_AS_events[,list(event_id,alternative_transcripts,total_transcripts)],by.x= "Event",by.y="event_id",all.x=TRUE,all.y=FALSE)
RI_TGS1_dpsi_sig = RI_TGS1_dpsi[ pval < 0.05 & abs(dPSI)>0.1]
RI_TGS1_dpsi_sig[,major_isoform_down_reg:=0]
RI_TGS1_dpsi_sig[,major_isoform_up_reg:=0]
for (counter in 1:nrow(RI_TGS1_dpsi_sig)){
  alternative_transcripts = strsplit(RI_TGS1_dpsi_sig$alternative_transcripts[counter],split=",",fixed=TRUE)[[1]] # this transcript has the retained intron
  total_transcripts = strsplit(RI_TGS1_dpsi_sig$total_transcripts[counter],split=",",fixed=TRUE)[[1]] # all transcripts involved in the RI event 
  major_isoforms = setdiff(total_transcripts,alternative_transcripts)  # transcripts that do not have the intron
  if (major_isoforms %in% DE_transcripts_TGS1[b>0]$target_id){
    RI_TGS1_dpsi_sig[counter,major_isoform_down_reg:=1]
  }
  if (major_isoforms %in% DE_transcripts_TGS1[b<0]$target_id){
    RI_TGS1_dpsi_sig[counter,major_isoform_up_reg:=1]
  }
}
    
RI_SMN_dpsi[,c("V5","V6","V7","V8","int_pos","int_pos1","int_pos2","chr","intron","gene_name","ENSG"):=NULL]
RI_SMN_dpsi = merge(RI_SMN_dpsi,SUPPA_AS_events[,list(event_id,alternative_transcripts,total_transcripts)],by.x= "Event",by.y="event_id",all.x=TRUE,all.y=FALSE)
RI_SMN_dpsi_sig = RI_SMN_dpsi[ pval < 0.05 & abs(dPSI)>0.1]
RI_SMN_dpsi_sig[,major_isoform_down_reg:=0]
RI_SMN_dpsi_sig[,major_isoform_up_reg:=0]
for (counter in 1:nrow(RI_SMN_dpsi_sig)){
  alternative_transcripts = strsplit(RI_SMN_dpsi_sig$alternative_transcripts[counter],split=",",fixed=TRUE)[[1]] # this transcript has the retained intron
  total_transcripts = strsplit(RI_SMN_dpsi_sig$total_transcripts[counter],split=",",fixed=TRUE)[[1]] # all transcripts involved in the RI event 
  major_isoforms = setdiff(total_transcripts,alternative_transcripts)  # transcripts that do not have the intron
  if (major_isoforms %in% DE_transcripts_SMN[b>0]$target_id){
    RI_SMN_dpsi_sig[counter,major_isoform_down_reg:=1]
  }
  if (major_isoforms %in% DE_transcripts_SMN[b<0]$target_id){
    RI_SMN_dpsi_sig[counter,major_isoform_up_reg:=1]
  }
}

## now we want to add the TPM information from kallisto to the major transcripts for each RI event
directory = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\"
sample_id = sample_list$ID

kallisto_output_folders = list.files(directory, pattern = "kallisto_output_union_scallop_grazia_uniq", all.files = TRUE,recursive = FALSE)
kal_output_files = paste(directory,kallisto_output_folders,"\\abundance.tsv",sep="")

transcript_count = fread(kal_output_files[1],sep="\t",header=TRUE,select = c("target_id","tpm"))
consolidated_transcript_count  = transcript_count

for (counter in 2:length(kal_output_files)){
  transcript_count = fread(kal_output_files[counter],sep="\t",header=TRUE,select = c("tpm"))
  consolidated_transcript_count = cbind(consolidated_transcript_count,transcript_count)
}

consolidated_transcript_count[,target_id:=strsplit(as.character(target_id),split="|",fixed=TRUE)[[1]][1],by=1:nrow(consolidated_transcript_count)]
names(consolidated_transcript_count) = c("target_id",sample_id)


## plot the volcano plots for the intron retentions
RI_TGS1_dpsi[,transcript := strsplit(Event,split=";",fixed=TRUE)[[1]][1],by=Event]
RI_SMN_dpsi[,transcript := strsplit(Event,split=";",fixed=TRUE)[[1]][1],by=Event]
setnames(AS_diff_WT_vs_SMN_avglogtpm,"V2","average_logTPM")
setnames(AS_diff_WT_vs_TGS1_avglogtpm,"V2","average_logTPM")
RI_TGS1_dpsi = merge(RI_TGS1_dpsi,AS_diff_WT_vs_TGS1_avglogtpm[,list(V1,average_logTPM)],all.x=TRUE,all.y=FALSE,by.x="Event",by.y="V1")
RI_TGS1_dpsi[,sig:="not sig"]
RI_TGS1_dpsi[pval<0.05 & dPSI > 0.1,sig:="Up"]
RI_TGS1_dpsi[pval<0.05 & dPSI < -0.1,sig:="Down"]
RI_TGS1_dpsi[,log10pval:=-log10(pval),by=pval]
p <- ggplot(RI_TGS1_dpsi, aes(x=average_logTPM, y=dPSI, color=sig))
p + geom_point() + scale_color_manual(values=c("Up" = "red", "not sig" = "black", "nan" = "gray80","Down" = "blue"))
graphics.off()

RI_SMN_dpsi = merge(RI_SMN_dpsi,AS_diff_WT_vs_SMN_avglogtpm[,list(V1,average_logTPM)],all.x=TRUE,all.y=FALSE,by.x="Event",by.y="V1")
RI_SMN_dpsi[,sig:="not sig"]
RI_SMN_dpsi[pval<0.05 & dPSI > 0.1,sig:="Up"]
RI_SMN_dpsi[pval<0.05 & dPSI < -0.1,sig:="Down"]
RI_SMN_dpsi[,log10pval:=-log10(pval),by=pval]
p <- ggplot(RI_SMN_dpsi, aes(x=average_logTPM, y=dPSI, color=sig))
p + geom_point() + scale_color_manual(values=c("Up" = "red", "not sig" = "black", "nan" = "gray80","Down" = "blue"))



#########
### here we want to plot a volcano plot for the dTPM vs dPSI for intron retention events
#########


for (counter in 1:nrow(RI_TGS1_dpsi)){
  alternative_transcripts = strsplit(as.character(RI_TGS1_dpsi$alternative_transcripts[counter]),split=",",fixed=TRUE)[[1]] # this transcript has the retained intron
  total_transcripts = strsplit(as.character(RI_TGS1_dpsi$total_transcripts[counter]),split=",",fixed=TRUE)[[1]] # all transcripts involved in the RI event 
  major_isoforms = setdiff(total_transcripts,alternative_transcripts)  # transcripts that do not have the intron
  RI_TGS1_dpsi[counter,TPM_WT:=paste(as.vector(unlist(consolidated_transcript_count[target_id%in%major_isoforms,c(2,3)])),collapse=",")]
  RI_TGS1_dpsi[counter,TPM_TGS1:=paste(as.vector(unlist(consolidated_transcript_count[target_id%in%major_isoforms,c(4,5)])),collapse=",")]
}
RI_TGS1_dpsi_new = RI_TGS1_dpsi[TPM_TGS1!="" & TPM_WT!=""]
RI_TGS1_dpsi_new[,ave_TPM_TGS1:=mean(as.numeric(strsplit(TPM_TGS1,split = ",",fixed=TRUE)[[1]])),by=TPM_TGS1]
RI_TGS1_dpsi_new[,ave_TPM_WT:=mean(as.numeric(strsplit(TPM_WT,split = ",",fixed=TRUE)[[1]])),by=TPM_WT]
RI_TGS1_dpsi_new[,dTPM:=ave_TPM_TGS1 - ave_TPM_WT,by=1:nrow(RI_TGS1_dpsi_new)]
RI_TGS1_dpsi_new = setorder(RI_TGS1_dpsi_new,-sig)
p <- ggplot(RI_TGS1_dpsi_new[abs(dTPM)<25], aes(x=dTPM, y=dPSI, color=sig))
p + geom_point() + scale_color_manual(values=c("Up" = "chocolate3", "not sig" = "black", "nan" = "gray80","Down" = "cyan4"))  +  theme_bw()


## now the same for SMN
for (counter in 1:nrow(RI_SMN_dpsi)){
  alternative_transcripts = strsplit(RI_SMN_dpsi$alternative_transcripts[counter],split=",",fixed=TRUE)[[1]] # this transcript has the retained intron
  total_transcripts = strsplit(RI_SMN_dpsi$total_transcripts[counter],split=",",fixed=TRUE)[[1]] # all transcripts involved in the RI event 
  major_isoforms = setdiff(total_transcripts,alternative_transcripts)  # transcripts that do not have the intron
  RI_SMN_dpsi[counter,TPM_WT:=paste(as.vector(unlist(consolidated_transcript_count[target_id%in%major_isoforms,c(2,3)])),collapse=",")]
  RI_SMN_dpsi[counter,TPM_SMN:=paste(as.vector(unlist(consolidated_transcript_count[target_id%in%major_isoforms,c(4,5)])),collapse=",")]
}
RI_SMN_dpsi_new = RI_SMN_dpsi[TPM_SMN!="" & TPM_WT!=""]
RI_SMN_dpsi_new[,ave_TPM_SMN:=mean(as.numeric(strsplit(TPM_SMN,split = ",",fixed=TRUE)[[1]])),by=TPM_SMN]
RI_SMN_dpsi_new[,ave_TPM_WT:=mean(as.numeric(strsplit(TPM_WT,split = ",",fixed=TRUE)[[1]])),by=TPM_WT]
RI_SMN_dpsi_new[,dTPM:=ave_TPM_SMN - ave_TPM_WT,by=1:nrow(RI_SMN_dpsi_new)]
RI_SMN_dpsi_new = setorder(RI_SMN_dpsi_new,-sig)
p <- ggplot(RI_SMN_dpsi_new[abs(dTPM)<25], aes(x=dTPM, y=dPSI, color=sig))
p + geom_point() + scale_color_manual(values=c("Up" = "chocolate3", "not sig" = "black", "nan" = "gray80","Down" = "cyan4"))  +  theme_bw()



## build the final supplementary excel file for the significant AS events
SUPPA_TGS1_sig_AS_events = SUPPA_WT_vs_TGS1_dpsi[ pval < 0.05 & abs(dPSI)>0.1]
SUPPA_TGS1_sig_AS_events[,gene:=strsplit(Event,split=";",fixed=TRUE)[[1]][1],by=1:nrow(SUPPA_TGS1_sig_AS_events)]
SUPPA_TGS1_sig_AS_events[,AS_event_type:=strsplit(Event,split="[;:]")[[1]][2],by=1:nrow(SUPPA_TGS1_sig_AS_events)]
SUPPA_TGS1_sig_AS_events[is.na(gene_name),gene_name:=gene]
SUPPA_TGS1_sig_AS_events[,c("gene","ENSG"):=NULL]
SUPPA_TGS1_sig_AS_events = merge(SUPPA_TGS1_sig_AS_events,SUPPA_AS_events[,list(alternative_transcripts,total_transcripts,event_pos)],all.x = TRUE,all.y=FALSE,by.x="event_pos",by.y="event_pos")
SUPPA_TGS1_sig_AS_events = SUPPA_TGS1_sig_AS_events[!duplicated(Event)]
SUPPA_TGS1_sig_AS_events = merge(SUPPA_TGS1_sig_AS_events,SUPPA_WT_vs_TGS1_psivec[,list(event_pos,WT_1,WT_2,TGS1_1,TGS1_2)],all.x=TRUE,all.y=FALSE,by.x="event_pos",by.y="event_pos")
SUPPA_TGS1_sig_AS_events[,c("novel","event_pos"):=NULL]
setcolorder(SUPPA_TGS1_sig_AS_events,c("Event","gene_name","AS_event_type","dPSI","pval","WT_1","WT_2","TGS1_1","TGS1_2","alternative_transcripts","total_transcripts"))

SUPPA_SMN_sig_AS_events = SUPPA_WT_vs_SMN_dpsi[ pval < 0.05 & abs(dPSI)>0.1]
SUPPA_SMN_sig_AS_events[,gene:=strsplit(Event,split=";",fixed=TRUE)[[1]][1],by=1:nrow(SUPPA_SMN_sig_AS_events)]
SUPPA_SMN_sig_AS_events[,AS_event_type:=strsplit(Event,split="[;:]")[[1]][2],by=1:nrow(SUPPA_SMN_sig_AS_events)]
SUPPA_SMN_sig_AS_events[is.na(gene_name),gene_name:=gene]
SUPPA_SMN_sig_AS_events[,c("gene","ENSG"):=NULL]
SUPPA_SMN_sig_AS_events = merge(SUPPA_SMN_sig_AS_events,SUPPA_AS_events[,list(alternative_transcripts,total_transcripts,event_pos)],all.x = TRUE,all.y=FALSE,by.x="event_pos",by.y="event_pos")
SUPPA_SMN_sig_AS_events = SUPPA_SMN_sig_AS_events[!duplicated(Event)]
SUPPA_SMN_sig_AS_events = merge(SUPPA_SMN_sig_AS_events,SUPPA_WT_vs_SMN_psivec[,list(event_pos,WT_1,WT_2,SMN_1,SMN_2)],all.x=TRUE,all.y=FALSE,by.x="event_pos",by.y="event_pos")
SUPPA_SMN_sig_AS_events[,c("novel","event_pos"):=NULL]
setcolorder(SUPPA_SMN_sig_AS_events,c("Event","gene_name","AS_event_type","dPSI","pval","WT_1","WT_2","SMN_1","SMN_2","alternative_transcripts","total_transcripts"))

write.table(SUPPA_SMN_sig_AS_events,"G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\outputs\\AS_WT_vs_SMN.tsv",row.names=FALSE,quote=FALSE,sep="\t")
write.table(SUPPA_TGS1_sig_AS_events,"G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\outputs\\AS_WT_vs_TGS1.tsv",row.names=FALSE,quote=FALSE,sep="\t")