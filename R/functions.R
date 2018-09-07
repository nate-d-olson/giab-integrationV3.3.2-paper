## Workflow functions

## Utilities ###################################################################
get_file_list <- function(input_df, type){
    file_df <- dplyr::filter(input_df, file_type == type)
    as.list(file_df$file_source) %>%
        set_names(file_df$hgref) 
}

## File I/O ####################################################################
## Download file from url and saves to tmp file or user specified file name
get_from_ftp <- function(ftp, filename = NULL){
    if(grepl(x = ftp, pattern = "bed")){
        file_ext <- "bed"
    } else if( grepl(x = ftp, pattern = "vcf.gz")){
        file_ext <- "vcf.gz"
    } else {
        message("file extention not `bed` or `vcf.gz`")
        file_ext <- ""
    }
    
    if (is.null(filename)){
        filename <- tempfile(fileext = file_ext)
    }
    
    download.file(ftp, filename)
    
    filename
}

get_file_path <- function(file_source){
    if(grepl("ftp", file_source)){
        return(get_from_ftp(file_source))
    } else if(file.exists(file_source)){
        return(file_source)
    } else {
        stop(cat(file_source, " not ftp site or not a valid file path"))
    }
}

## Utilities ###################################################################
## Get non-N base genome and chromosome sizes, 
## requires BSgenome formatted object

get_chrom_sizes <- function(genome){
    
    get_alpha_freq <- function(i){
        genome[[i]] %>% 
            alphabetFrequency() %>% 
            as.workflow_data.frame()
    }
    
    alpha_freq_df <- as.list(1:22) %>% 
        map_dfc(get_alpha_freq)
    
    colnames(alpha_freq_df) <- paste0("chr",1:22)
    
    alpha_freq_df <- alpha_freq_df %>% 
        ## Removing bases not included in counts
        filter(chr1 != 0) %>% 
        add_column(base = c("A","C","G","T","N"))
    
    chromosome_lengths <- alpha_freq_df %>% 
        tidyr::gather(key = "chrom", value = "nbases", -base) %>% 
        group_by(chrom) %>% 
        mutate(base_type = if_else(base == "N", "N", "non_N")) %>% 
        group_by(chrom, base_type) %>% 
        summarise(n_bases = sum(nbases)) %>% 
        tidyr::spread(base_type, n_bases) %>% 
        mutate(len = N + non_N) %>% 
        dplyr::select(-N)
    
    ## workflow_data frame with total length and number of non-N bases
    workflow_data_frame(chrom = "genome", 
               non_N = sum(chromosome_lengths$non_N),
               len = sum(chromosome_lengths$len)) %>% 
        bind_rows(chromosome_lengths)
}



### High confidence coverage ###################################################
### Bed chromosome coverage lengths from bed
get_bed_cov_by_chrom <- function(bed_file){
    ## Read as tsv
    bed_df <- read_tsv(bed_file, 
                       col_names = c("chrom","start","end")) %>% 
        ## Compute region size
        mutate(region_size = end - start)
    
    ## Compute bases per chromosome
    chrom_cov <- bed_df %>% 
        group_by(chrom) %>% 
        summarise(nbases = sum(region_size))

    ## NBases workflow_data frame
    workflow_data_frame(chrom = "genome",
               nbases = sum(chrom_cov$nbases)) %>% 
        bind_rows(chrom_cov)
}

## Compute coverage of total size and non-N bases
get_hc_cov <- function(bed_source, chrom_lengths){
    ## get file path
    bed_file <- get_file_path(bed_source)
    
    cov_df <- get_chrom_cov_df(bed_file)
    
    cov_df %>% 
        ## add Chrom size and non-N base size
        left_join(chrom_lengths) %>% 
        ## Calculate coverage
        mutate(coverage = hc_size/len,
               coverage_nonN = hc_size/non_N)
}

### Generating VCF stats by chromosome ######################################### 
make_stats_df <- function(vcf_stats_dir){
    ## TODO -----------------------------
    ## list files
    ## read as workflow_data frame
    ## tidy
    ## return workflow_data frame
    workflow_data_frame()
}

get_vcf_stats <- function(vcf_source, vcf_type){
    vcf_file <- get_file_path(vcf_source)
    
    ## Defining and creating output directory
    hgref <- str_extract(vcf_source, "HG00.")
    vcf_stats_dir <- paste0("workflow_data/", hgref, "_", vcf_type, "_stats")
    if(!dir.exists(vcf_stats_dir)) dir.create(vcf_stats_dir)
    
    system(paste('bash bash/get_vcf_chrom_stats.sh ',
                 vcf_file,
                 vcf_stats_dir),
           intern = FALSE)
    
    ## Combine stats into a single workflow_data frame
    make_stats_df(vcf_stats_dir)
}

### High confidence variants in high confidence regions
get_hh_stats_df <- function(vcf_source, bed_source){
    ## VCF bed intersect
    bed_file <- get_file_path(bed_source)
    vcf_file <- get_file_path(vcf_source)
    
    # output directory
    hgref <- str_extract(vcf_source, "HG00.")
    hh_var_dir <- paste0(hgref, "_", vcf_type, "_stats")
    dir.create(paste0("workflow_data/", hh_var_dir))
    
    system(paste('bash bash/get_highhigh_vars.sh',
                 bed_file,
                 vcf_file, 
                 hh_var_dir),
           intern = FALSE)
    
    ## Get vcf stats by chromosome
    hh_vcf <- file.path("workflow_data", hh_var_dir, "highhigh.vcf")
    get_vcf_stats(h, hgref, vcf_type = "hh")
    
    ## Combine in to data frame
    ## TODO
    
    ## empty return for testing
    data_frame()
}

