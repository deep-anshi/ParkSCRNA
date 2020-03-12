########## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Creating different boxplot functions to be used in the later part of differentail epression
# gene analysis
# 
# Lab: Havrada 
# Author: Deepanshi Shokeen

########## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Boxplot for nlyzing sex in control and pd
box_plot_sex<- function(gene_name, data, out_dir, out_name, title){
  counts <- data@assays$SCT@counts[gene_name, rownames(data@meta.data)]
  counts[counts==0] <- NA
  cells <- rownames(pd.combined@meta.data)
  sex <- pd.combined@meta.data$sex
  sample <- pd.combined@meta.data$sample
  subject <- pd.combined@meta.data$subject
  df <- data.frame(counts, cells, sex, sample, subject)
  ppi=300
  png(paste0(out_dir, out_name), width=10*ppi, height=8.75*ppi, res=ppi)
  p<-ggplot(df, aes(x=sample, y=counts, color=sex)) + 
    geom_boxplot(width=0.3)+
    ggtitle(title)+
    labs(x="Sample group", y="Normalized Counts")
  print(p)
  dev.off()
}



# Boxplot for nlyzing different each sample in control and pd differenciated colrs for sex
box_plot_sex_sample<- function(gene_name, data, out_dir, out_name, title){
  counts <- data@assays$SCT@counts[gene_name, rownames(data@meta.data)]
  counts[counts==0] <- NA
  cells <- rownames(pd.combined@meta.data)
  sub_sample <- pd.combined@meta.data$subject2
  sample <- pd.combined@meta.data$sample
  sex <- pd.combined@meta.data$sex
  subject <- pd.combined@meta.data$subject
  df <- data.frame(counts, cells, sub_sample, sample, subject, sex)
  ppi=300
  png(paste0(out_dir, out_name), width=10*ppi, height=8.75*ppi, res=ppi)
  p<-ggplot(df, aes(x=sub_sample, y=counts, color=sex)) + 
    geom_boxplot(width=0.3)+
    ggtitle(title)+
    labs(x="Sample group", y="Normalized Counts")
  print(p)
  dev.off()
}


# Boxplot for nlyzing different each sample in control and pd differenciated colrs for sex
box_plot_sample<- function(gene_name, data, out_dir, out_name, title){
  counts <- data@assays$SCT@counts[gene_name, rownames(data@meta.data)]
  counts[counts==0] <- NA
  cells <- rownames(pd.combined@meta.data)
  sub_sample <- pd.combined@meta.data$subject2
  sample <- pd.combined@meta.data$sample
  subject <- pd.combined@meta.data$subject
  df <- data.frame(counts, cells, sub_sample, sample, subject)
  ppi=300
  png(paste0(out_dir, out_name), width=10*ppi, height=8.75*ppi, res=ppi)
  p<-ggplot(df, aes(x=sub_sample, y=counts, color=sample)) + 
    geom_boxplot(width=0.3)+
    ggtitle(title)+
    labs(x="Sample group", y="Normalized Counts")
  print(p)
  dev.off()
}

