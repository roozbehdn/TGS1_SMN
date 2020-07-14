library(sleuth)
library(data.table)
library(dplyr)
library(pheatmap)




### input files #########
directory = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\"
sample_list_file = "G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\utility_files\\sample_list.txt"
load("/oak/stanford/groups/horence/Roozbeh/Lu_data/outputs/Transcript_scallop_DE_WT_vs_TGS1.RData")
scallop_transcripts_tmap_file = "/oak/stanford/groups/horence/Roozbeh/Lu_data/gffall.novel_transcripts_grazia.gtf.tmap"
transcripts_classification_file = "/oak/stanford/groups/horence/Roozbeh/Lu_data/excel_files/classification.csv"
transcripts_classification_Roozbeh_file = "/oak/stanford/groups/horence/Roozbeh/Lu_data/excel_files/classification_Roozbeh.tsv"
old_DE_transcripts_TGS1 = fread("/oak/stanford/groups/horence/Roozbeh/Lu_data/excel_files/TGS1_vs_CTRL.csv",sep=",",header=TRUE)
transcript_id_gene_name_file = "/oak/stanford/groups/horence/Roozbeh/Lu_data/transcript_id_gene_name.tsv" # I need this file to add gene names for annotated DE transcripts
############################

## read in input files #######
sample_list = fread(sample_list_file,sep="\t",header=TRUE)
transcripts_classification = fread(transcripts_classification_file,sep=",",header=TRUE)
transcripts_classification_Roozbeh = fread(transcripts_classification_Roozbeh_file,sep="\t",header=TRUE)
scallop_transcripts_tmap = fread(scallop_transcripts_tmap_file,sep="\t",header = TRUE)
transcript_id_gene_name = fread(transcript_id_gene_name_file,sep="\t",header=TRUE)
##############################


kallisto_output_folders = list.files(directory, pattern = "kallisto_output_union_scallop_grazia_uniq", all.files = TRUE,recursive = FALSE)[c(1,2,3,4)]
kal_dirs = paste(directory,kallisto_output_folders,sep="")
sample_id = sample_list$ID[c(1,2,3,4)]
s2c = data.table(sample_id,condition = c("WT","WT","TGS1","TGS1"))
colnames(s2c)=c("sample","condition")
s2c = cbind(s2c,path = kal_dirs)

so_wt_tgs <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, transformation_function = function(x) log2(x + 0.5))    # initializing the sleuth object
so_wt_tgs <- sleuth_fit(so_wt_tgs, ~condition, 'full')  # the full model with a parameters that represents the experimental condition
so_wt_tgs <- sleuth_fit(so_wt_tgs, ~1, 'reduced')   # the reduced model with the assumption that abudnances are equal
so_wt_tgs <- sleuth_lrt(so_wt_tgs, 'reduced', 'full')   # the actual test

#sleuth_table <- sleuth_results(so_wt_tgs, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
#head(sleuth_significant, 20)

## b>0 means that the transcripts is downregulated in the mutated sample (example: ENST00000373247.7|ENSG00000136877.15|OTTHUMG00000020716.6|OTTHUMT00000054251.3|FPGS-203|FPGS|2284|protein_coding|  )
## now we do the wald test
so_wt_tgs <- sleuth_wt(so_wt_tgs,'conditionWT','full')
id_dt=data.table(so_wt_tgs$obs_norm$target_id)
id_dt[,V2:=strsplit(V1,split="|",fixed=TRUE)[[1]][1],by=1:nrow(id_dt)]
so_wt_tgs$obs_norm$target_id=id_dt$V2
id_dt=data.table(so_wt_tgs$obs_raw$target_id)
id_dt[,V2:=strsplit(V1,split="|",fixed=TRUE)[[1]][1],by=1:nrow(id_dt)]
so_wt_tgs$obs_raw$target_id=id_dt$V2


results_table_wald_wt_tgs = sleuth_results(so_wt_tgs, 'conditionWT','full', test_type = 'wt')
results_table_wald_wt_tgs = data.table(results_table_wald_wt_tgs)
results_wald_significant_wt_tgs = results_table_wald_wt_tgs[qval <= 0.05 & abs(b) >= 2,list(target_id,mean_obs,b,pval,qval)]
results_wald_significant_wt_tgs[,transcript:=strsplit(target_id,split="|",fixed=TRUE)[[1]][1],by=target_id]
results_wald_significant_wt_tgs[,target_id:=NULL]
## plot_bootstrap(so, "ENST00000263734", units = "est_counts", color_by = "condition")

## finding upregulated novel transcripts with down regualted annotated transcripts

results_table_wald_wt_tgs[,target_id:=as.character(strsplit(target_id,split="|",fixed=TRUE)[[1]][1]),by=1:nrow(results_table_wald_wt_tgs)]

scallop_transcripts_tmap = merge(scallop_transcripts_tmap,results_table_wald_wt_tgs[,list(target_id,qval,b)],by.x="qry_id",by.y="target_id",all.x=TRUE,all.y=FALSE)
setnames(scallop_transcripts_tmap,old=c("qval","b"),new=c("qval_novel","b_novel"))
scallop_transcripts_tmap = merge(scallop_transcripts_tmap,results_table_wald_wt_tgs[,list(target_id,qval,b)],by.x="ref_id",by.y="target_id",all.x=TRUE,all.y=FALSE)
setnames(scallop_transcripts_tmap,old=c("qval","b"),new=c("qval_ref","b_ref"))
upreg_novel_with_downreg_ref_trans_wt_vt_tgs1 = scallop_transcripts_tmap[qval_novel<0.05 & b_novel<=-2 & qval_ref<0.05 & b_ref>=2]
write.table(upreg_novel_with_downreg_ref_trans_wt_vt_tgs1,"G:\\Shared drives\\Salzman Lab Team Drive\\Members\\Roozbeh\\Projects\\Lu_data\\outputs\\upreg_novel_with_downreg_ref_trans_wt_vt_tgs1.csv",quote=FALSE,row.names=FALSE,sep=",")


## compare with old DE transcripts and write output files for DE transcripts
old_DE_transcripts_TGS1 = old_DE_transcripts_TGS1[padj < 0.05 & abs(log2FoldChange)>2]
old_DE_transcripts_TGS1[,transcript_new:=transcript]
old_DE_transcripts_TGS1[transcript%like%"ENST",transcript_new:=strsplit(transcript,split=".",fixed=TRUE)[[1]][1],by=transcript]

out_dt= data.table(so_wt_tgs$tests$wt$full$conditionWT)
out_dt[,target_id:=strsplit(target_id,split="|",fixed=TRUE)[[1]][1],by=1:nrow(out_dt)]
out_dt_sig = out_dt[target_id%in%results_wald_significant_wt_tgs$transcript]
out_dt_sig[,target_id_new:=target_id]
out_dt_sig[target_id%like%"ENST",target_id_new:=strsplit(target_id_new,split=".",fixed=TRUE)[[1]][1],by=target_id]
out_dt_sig[,c("rss","iqr","failed_ise","smooth_sigma_sq","sigma_sq_pmax","smooth_sigma_sq_pmax","wald_stat","sigma_q_sq"):=NULL]
old_DE_transcripts_TGS1[,found_by_sleuth:=0]
old_DE_transcripts_TGS1[transcript_new %in% out_dt_sig$target_id_new, found_by_sleuth:=1]

out_dt_sig[,found_old_method:=0]
out_dt_sig[target_id_new %in% old_DE_transcripts_TGS1$transcript_new, found_old_method:=1]
out_dt_sig[,target_id_new:=NULL]
old_DE_transcripts_TGS1[,transcript_new:=NULL]

write.table(out_dt_sig,"/oak/stanford/groups/horence/Roozbeh/Lu_data/outputs/transcripts_DE_analysis_WT_vs_TGS1.csv",row.names=FALSE,sep=",",quote=FALSE)
write.table(old_DE_transcripts_TGS1,"/oak/stanford/groups/horence/Roozbeh/Lu_data/outputs/TGS1_vs_CTRL_Roozbeh.csv",row.names=FALSE,sep=",",quote=FALSE)


## plot heatmap with a side colorbar for the type of significant transcript
#first I want to combine the two classification files I have based on my results and grazia results
transcripts_classification_Roozbeh = transcripts_classification_Roozbeh[!novel%in%transcripts_classification$novel]
transcripts_classification = rbind(transcripts_classification[,list(novel,similar_to,name,utr3)],transcripts_classification_Roozbeh[,list(novel,similar_to,name,utr3)])
results_wald_significant_wt_tgs = merge(results_wald_significant_wt_tgs,transcripts_classification[!duplicated(novel)],all.x=TRUE,all.y=FALSE,by.x="transcript",by.y="novel")
results_wald_significant_wt_tgs[transcript%like%"ENS",transcript_type:="Annotated transcript"]
results_wald_significant_wt_tgs[transcript%like%"gene" & utr3=="longer" , transcript_type:="Novel with 3' extension"]
results_wald_significant_wt_tgs[transcript%like%"gene" & utr3=="shorter", transcript_type:="Novel with 3' shortening"]
results_wald_significant_wt_tgs[transcript%like%"gene" & utr3=="similar", transcript_type:="Novel with similar 3'"]
results_wald_significant_wt_tgs[is.na(similar_to)&!transcript%like%"ENS",transcript_type:="Novel intergenic"]
results_wald_significant_wt_tgs[!is.na(similar_to)&!transcript%like%"ENS" & is.na(utr3),transcript_type:="Novel opposite strand"]
results_wald_significant_wt_tgs[ transcript %like% "ENST",transcript := strsplit(transcript,split=".",fixed=TRUE)[[1]][1],by=transcript]
results_wald_significant_wt_tgs = merge(results_wald_significant_wt_tgs,transcript_id_gene_name,all.x=TRUE,all.y=FALSE,by.x="transcript",by.y="transcript_id")  # here we add gene names for all DE transcripts
results_wald_significant_wt_tgs[transcript %like% "ENST",name:=gene_name]
results_wald_significant_wt_tgs[,gene_name:=NULL]
write.table(results_wald_significant_wt_tgs,"/oak/stanford/groups/horence/Roozbeh/Lu_data/outputs/DE_transcripts_WT_vs_TGS1.csv",sep=",",row.names=FALSE,quote=FALSE)
results_wald_significant_wt_tgs[,heatmap_type:=transcript_type]
results_wald_significant_wt_tgs[transcript_type%like%"intergen" | transcript_type%like%"opposite",heatmap_type:="Other"]
results_wald_significant_wt_tgs$heatmap_type = as.factor(results_wald_significant_wt_tgs$heatmap_type)
results_wald_significant_wt_tgs[b>0,status:="down_regulated"]
results_wald_significant_wt_tgs[b<0,status:="up_regulated"]
results_wald_significant_wt_tgs = data.frame(results_wald_significant_wt_tgs)
rownames(results_wald_significant_wt_tgs)=results_wald_significant_wt_tgs$transcript
results_wald_significant_wt_tgs_dt = data.table(results_wald_significant_wt_tgs)

source("/oak/stanford/groups/horence/Roozbeh/Lu_data/scripts/plot_heat_map.R")
Breaks <- seq(0, 9, length = 100)
graphics.off()
pdf('/oak/stanford/groups/horence/Roozbeh/Lu_data/outputs/TGS1_transcripts_heatmp.pdf')
ann_colors = list(heatmap_type  = c(`Annotated transcript` = "tan2", `Novel with 3' extension` = "violetred2", `Novel with 3' shortening` = "cornflowerblue", `Novel with similar 3'` = "plum3", `Other` = "mediumpurple4" ),condition = c(TGS1="#1B9E77", WT = "tomato1"),status = c(down_regulated = "palegreen2",up_regulated = "orangered3"))
plot_heat_map(so_wt_tgs, transcripts = rownames(results_wald_significant_wt_tgs),show_rownames=FALSE,annotation_row = results_wald_significant_wt_tgs,annotation_colors = ann_colors,breaks =Breaks)
dev.off()
pdf('/oak/stanford/groups/horence/Roozbeh/Lu_data/outputs/TGS1_annotated_transcripts_heatmp.pdf')
transcripts = rownames(results_wald_significant_wt_tgs)[which(results_wald_significant_wt_tgs$heatmap_type== "Annotated transcript")]
plot_heat_map(so_wt_tgs, transcripts = transcripts,show_rownames=FALSE,annotation_row = data.frame(heatmap_type=c(rep("Annotated transcript",length(transcripts))),status = results_wald_significant_wt_tgs_dt[ transcript %in%transcripts]$status,row.names = transcripts),annotation_colors = ann_colors,breaks =Breaks)
dev.off()
pdf('/oak/stanford/groups/horence/Roozbeh/Lu_data/outputs/TGS1_transcripts_3_extended_heatmp.pdf')
transcripts = rownames(results_wald_significant_wt_tgs)[which(results_wald_significant_wt_tgs$heatmap_type== "Novel with 3' extension")]
plot_heat_map(so_wt_tgs, transcripts = transcripts,show_rownames=FALSE,annotation_row = data.frame(heatmap_type=c(rep("Novel with 3' extension",length(transcripts))),status = results_wald_significant_wt_tgs_dt[ transcript %in%transcripts]$status,row.names = transcripts),annotation_colors = ann_colors,breaks =Breaks)
dev.off()
pdf('/oak/stanford/groups/horence/Roozbeh/Lu_data/outputs/TGS1_transcripts_short_3_heatmp.pdf')
transcripts = rownames(results_wald_significant_wt_tgs)[which(results_wald_significant_wt_tgs$heatmap_type== "Novel with 3' shortening")]
plot_heat_map(so_wt_tgs, transcripts = transcripts,show_rownames=FALSE,annotation_row = data.frame(heatmap_type=c(rep("Novel with 3' shortening",length(transcripts))),status = results_wald_significant_wt_tgs_dt[ transcript %in%transcripts]$status,row.names = transcripts),annotation_colors = ann_colors,breaks =Breaks)
dev.off()
pdf('/oak/stanford/groups/horence/Roozbeh/Lu_data/outputs/TGS1_transcripts_similar_3_heatmp.pdf')
transcripts = rownames(results_wald_significant_wt_tgs)[which(results_wald_significant_wt_tgs$heatmap_type== "Novel with similar 3'")]
plot_heat_map(so_wt_tgs, transcripts = transcripts,show_rownames=FALSE,annotation_row = data.frame(heatmap_type=c(rep("Novel with similar 3'",length(transcripts))),status = results_wald_significant_wt_tgs_dt[ transcript %in%transcripts]$status,row.names = transcripts),annotation_colors = ann_colors,breaks =Breaks)
dev.off()
pdf('/oak/stanford/groups/horence/Roozbeh/Lu_data/outputs/TGS1_transcripts_other_heatmp.pdf')
transcripts = rownames(results_wald_significant_wt_tgs)[which(results_wald_significant_wt_tgs$heatmap_type== "Other")]
plot_heat_map(so_wt_tgs, transcripts = transcripts,show_rownames=FALSE,annotation_row = data.frame(heatmap_type=c(rep("Other",length(transcripts))),status = results_wald_significant_wt_tgs_dt[ transcript %in%transcripts]$status,row.names = transcripts),annotation_colors = ann_colors,breaks =Breaks)
dev.off()