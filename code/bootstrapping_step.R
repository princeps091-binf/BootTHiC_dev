library(GenomicRanges)
library(tidyverse)
library(valr)
library(furrr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(org.Hs.eg.db)

res_set <- c("1Mb", "500kb", "100kb", "50kb", "10kb", "5kb")
res_num <- c(1e6L, 5e5L, 1e5L, 5e4L, 1e4L, 5e3L)
names(res_num) <- res_set
###############################
#Snakemake
#clusters_folder <- snakemake@input
###############################
# Get Rda files in
get_obj_in_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
###############################
# Build GRange table
build_chr_cl_GRange_tbl<-function(res_folder,res_file,chromo){
  chr_spec_res<-get_obj_in_fn(paste0(res_folder,chromo,res_file))
  chr_cl_tbl<-tibble(chr=chromo,cl=names(chr_spec_res$cl_member),bins=lapply(chr_spec_res$cl_member,as.integer))
  chr_cl_tbl<-chr_cl_tbl %>% 
    mutate(res=str_split_fixed(cl,"_",2)[,1])
}
###############################
# Compute feature genome-annotation
build_rn_chr_content<-function(peakAnno,chr_feature_Grange,fn_file){
  rn_annotation<-sample(peakAnno@annoStat$Feature,size = length(chr_feature_Grange),prob = peakAnno@annoStat$Frequency/100,replace = T)
  #check number of peaks from that category
  n5<-length(grep("5'",as.character(rn_annotation)))
  n3<-length(grep("3'",as.character(rn_annotation)))
  nexon<-length(grep("Exon",as.character(rn_annotation)))
  nintron<-length(grep("Intron",as.character(rn_annotation)))
  n1kb<-length(grep("1kb",as.character(rn_annotation)))
  n2kb<-length(grep("2kb",as.character(rn_annotation)))
  n3kb<-length(grep("3kb",as.character(rn_annotation)))
  ndown<-length(grep("Down",as.character(rn_annotation)))
  ninter<-length(grep("Inter",as.character(rn_annotation)))
  
  n_vec<-c(n3,n5,ndown,nexon,ninter,nintron,n1kb,n2kb,n3kb)
  names(n_vec)<-fn_file
  
  rm(n5,n3,nexon,nintron,n1kb,n2kb,n3kb,ndown,ninter)
  n_vec<-n_vec[n_vec>0]
  return(n_vec)
  
}
###############################
# Generate random chromosome
rn_feature_GRange_build_fn<-function(n_vec,fn_bed_l,hg19_coord,tmp_cage_tbl,fn_file){
  
  
  rn_fn_coord_l<-vector('list',length(n_vec))
  names(rn_fn_coord_l)<-names(n_vec)
  for(f in names(n_vec)){
    # Consider repeating until valid shuffle? -> iterative try until ok
    tmp_n<-n_vec[f]
    if(grepl("inter_no.BED$",fn_file[f])){
      rn_fn_coord_l[[f]]<-structure(
        "message", 
        class = c("try-error", "character")
      )
      while("try-error" %in% class(rn_fn_coord_l[[f]])){
        rn_fn_coord_l[[f]]<-try(bed_shuffle(tmp_cage_tbl%>%sample_n(tmp_n),genome = hg19_coord,excl = fn_bed_l[[f]],within = T,max_tries=1e3),silent=T)
      }
      
    }
    if(!(grepl("inter_no.BED$",fn_file[f]))){
      rn_fn_coord_l[[f]]<-structure(
        "message", 
        class = c("try-error", "character")
      )
      while("try-error" %in% class(rn_fn_coord_l[[f]])){
        rn_fn_coord_l[[f]]<-try(valr::bed_shuffle(x = tmp_cage_tbl%>%sample_n(tmp_n),genome = hg19_coord,incl = fn_bed_l[[f]],within=T,max_tries=1e3),silent=T)
        
      }
      # Collect the successful shuffling by eliminating the shuffles producing try-error objects
    }
  }
  rn_fn_coord_tbl<-do.call(bind_rows,rn_fn_coord_l)
  rnp_Grange<-GRanges(seqnames=rn_fn_coord_tbl$chrom,
                      ranges = IRanges(start=rn_fn_coord_tbl$start,
                                       end=rn_fn_coord_tbl$end)
  )
  
  return(rnp_Grange)
}

###############################
# Filter input clusters to contain at least two feature-containing bins
filter_cluster_fn<-function(cl_GRange_tbl,chr_feature_GRange,chromo,nworker){
  tmp_bin<-cl_GRange_tbl %>% 
    group_by(res) %>% 
    summarise(bins=list(as.integer(unique(unlist(bins)))))
  tmp_bin<-tmp_bin %>% 
    mutate(cage.bin=pmap(list(res,bins),function(res,bins){
      bin_Grange<-GRanges(seqnames=chromo,
                          ranges = IRanges(start=bins,
                                           end=bins + res_num[res] -1)
      )
      return(countOverlaps(bin_Grange,chr_feature_GRange))
    }))
  tmp_bin<-tmp_bin %>% 
    mutate(cage.on=pmap(list(bins,cage.bin),function(bins,cage.bin){
      bins[which(cage.bin>0)]
    }))
  plan(multisession, workers=nworker)
  tmp_tbl<-cl_GRange_tbl %>% 
    mutate(ok.cl=future_pmap_lgl(list(res,bins),function(tmp_res,bins){
      res_bin<-tmp_bin %>% 
        filter(res==tmp_res) %>% 
        dplyr::select(cage.on) %>% 
        unnest(cols=c(cage.on)) %>% 
        unlist
      return(sum(as.integer(bins) %in% res_bin) > 1)
    }))
  plan(sequential)
  return(tmp_tbl)
}
###############################
# Load data
## load clusters
res_folder<-"~/Documents/multires_bhicect/data/GM12878/determinate_spec_res/"
res_file<-"_spec_res.Rda"

## load feature
feature_file<-"./data/CAGE_union_GM12878_Grange.Rda"
## Folder with functional bed files
fn_repo<-"~/Documents/multires_bhicect/data/epi_data/fn_BED/"
genome_file<-"~/Documents/multires_bhicect/data/hg19.genome"
###############################
# Genome-wide features
hg19_coord <- read_delim(genome_file, 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
names(hg19_coord)<-c("chrom","size")


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

feature_GRange<-get_obj_in_fn(feature_file)
###############################
chromo<-"chr22"
#cl_GRange_tbl<-get_obj_in_fn(paste0(cl_folder,chromo,cl_file))

chr_feature_GRange<-feature_GRange[seqnames(feature_GRange)==chromo]
peakAnno <- annotatePeak(chr_feature_GRange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)
rm(ChIPseekerEnv,envir = globalenv())

cl_res_tbl<-build_chr_cl_GRange_tbl(res_folder,res_file,chromo)

ok_cl_GRange_tbl<-filter_cluster_fn(cl_res_tbl,chr_feature_GRange,chromo,4) %>% 
  filter(ok.cl)

plan(cluster, workers = 3)
ok_cl_GRange_tbl<-ok_cl_GRange_tbl %>% mutate(GRange=future_pmap(list(chr,bins,res),function(chr,bins,res){
  return(GRanges(seqnames=chr,
                 ranges = IRanges(start=as.numeric(bins),
                                  end=as.numeric(bins)+res_num[res]-1
                 )))
  
  
}))
plan(sequential)

cl_list<-GRangesList(ok_cl_GRange_tbl$GRange)  

tmp_cage_tbl<-chr_feature_GRange %>% as_tibble %>% dplyr::select(seqnames,start,end)%>%dplyr::rename(chrom=seqnames)

# Compute feature genome-annotation
fn_folder<-paste0(fn_repo,chromo,"/")
fn_file<-grep('BED$',list.files(fn_folder),value = T)
fn_bed_l<-lapply(fn_file,function(f){
  read_bed(paste0(fn_folder,f),n_fields = 3)
})
names(fn_bed_l)<-fn_file

# Bootstrap
nboot<-1e4
obs_count<-countOverlaps(cl_list,chr_feature_GRange)
plan(cluster, workers=4)
boot_bool<-future_map(1:nboot,function(i){
  # Build random chromosome content
  n_vec<-build_rn_chr_content(peakAnno,chr_feature_GRange,fn_file)
  # Build random chromosome 
  rnp_Grange<-rn_feature_GRange_build_fn(n_vec,fn_bed_l,hg19_coord,tmp_cage_tbl,fn_file)
  # Compare with obs
  return(countOverlaps(cl_list,rnp_Grange) >= obs_count)
  
})
plan(sequential)


#Compute p-value
tmp_pval<-(apply(do.call(cbind,boot_bool),1,sum)+1)/(nboot+1)

ok_cl_GRange_tbl %>% 
  mutate(emp.pval=tmp_pval) %>% 
  ggplot(.,aes(-log10(emp.pval)))+
  geom_histogram()
