## Workflow functions

## Utilities ###################################################################
 
## Generate a named list of files for downstream methods
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

get_chrom_sizes <- function(genome = "hs37d5"){
    if (genome == "hs37d5"){
        require(BSgenome.Hsapiens.1000genomes.hs37d5)
        genome <- BSgenome.Hsapiens.1000genomes.hs37d5
    } else {
        stop("Only coded for hs37d5")
    }
    
    get_alpha_freq <- function(i){
        genome[[i]] %>% 
            alphabetFrequency() %>% 
            data.frame()
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
    
    ## data frame with total length and number of non-N bases
    data_frame(chrom = "genome", 
               non_N = sum(chromosome_lengths$non_N),
               len = sum(chromosome_lengths$len)) %>% 
        bind_rows(chromosome_lengths)
}



### High confidence coverage ###################################################
### Bed chromosome coverage lengths from bed
get_bed_cov_by_chrom <- function(bed_file){
    ## Read as tsv
    bed_df <- read_tsv(bed_file, 
                       col_names = c("chrom","start","end"), col_types = "cii") %>% 
        ## Compute region size
        mutate(region_size = end - start)
    
    ## Compute bases per chromosome
    chrom_cov <- bed_df %>% 
        ## Changing to chromosome names to chr 
        mutate(chrom = paste0("chr", chrom)) %>% 
        group_by(chrom) %>% 
        summarise(nbases = sum(region_size))

    ## NBases data frame
    data_frame(chrom = "genome",
               nbases = sum(chrom_cov$nbases)) %>% 
        bind_rows(chrom_cov)
}

## Compute coverage of total size and non-N bases
get_hc_cov <- function(bed_source, chrom_lengths){
    ## get file path
    bed_file <- get_file_path(bed_source)
    
    cov_df <- get_bed_cov_by_chrom(bed_file)
    
    cov_df %>% 
        ## add Chrom size and non-N base size
        left_join(chrom_lengths) %>% 
        ## Calculate coverage
        mutate(coverage = nbases/len,
               coverage_nonN = nbases/non_N)
}

### Generating VCF stats by chromosome ######################################### 
make_stats_df <- function(vcf_stats_dir){
    ## list files
    stat_files <- list.files(vcf_stats_dir, full.names = TRUE) %>% 
        set_names(str_extract(., "(?<=_stats.).*(?=.stats)"))
    
    ## read as data frame
    stat_df <- stat_files %>% 
        map_dfr(read_delim, delim = ":", 
                trim_ws = TRUE, 
                col_names = c("Key","Value"),
                .id = "chrom")

    ## tidy
    ## TODO -----------------------------
    
    ## return data frame
    stat_df
}

get_vcf_stats <- function(vcf_source, vcf_type){
    vcf_file <- get_file_path(vcf_source)
    
    ## Defining and creating output directory
    hgref <- str_extract(vcf_source, "HG00.")
    vcf_stats_dir <- paste0("data_workflow/", hgref, "_", vcf_type, "_stats")
    if(!dir.exists(vcf_stats_dir)) dir.create(vcf_stats_dir)
    
    system(paste('bash bash/get_vcf_chrom_stats.sh ',
                 vcf_file,
                 vcf_stats_dir),
           intern = FALSE)
    
    ## Combine stats into a single data frame
    make_stats_df(vcf_stats_dir)
}

### High confidence variants in high confidence regions ########################
get_hh_stats_df <- function(vcf_source, bed_source){
    ## VCF bed intersect
    bed_file <- get_file_path(bed_source)
    vcf_file <- get_file_path(vcf_source)
    
    # output directory
    hgref <- str_extract(vcf_source, "HG00.")
    hh_var_dir <- paste0("data_workflow/", hgref, "_hh_vcf")
    if(!dir.exists(hh_var_dir)) dir.create(hh_var_dir)
    
    system(paste('bash bash/get_highhigh_vars.sh',
                 bed_file,
                 vcf_file, 
                 hh_var_dir),
           intern = FALSE)
    
    ## Get vcf stats by chromosome
    hh_vcf <- file.path(hh_var_dir, "highhigh.vcf.gz")
    get_vcf_stats(hh_vcf, vcf_type = "hh")
}


## High conf fraction of RefSeq coding sequence covered ########################
get_coding_seq_cov <- function(bed_source, refseq_coding_bed){
    
    # Get bed file
    bed_file <- get_file_path(bed_source) 
    
    # Intersect highconf regions with refseq coding
    tmp_bed <- tempfile(fileext = "bed")
    system2("bedtools", 
            args = c("intersect", "-a", bed_file, "-b", refseq_coding_bed),
            stdout = tmp_bed)
    
    # refseq  nbases
    refseq_cov <- get_bed_cov_by_chrom(refseq_coding_bed) %>% 
        dplyr::rename(ncoding_bases = nbases)
    
    # intersect nbases
    intersect_cov <- get_bed_cov_by_chrom(tmp_bed)
    
    ## Calculate fraction
    left_join(refseq_cov, intersect_cov) %>% 
        mutate(nbases = if_else(!is.na(nbases), as.numeric(nbases), 0)) %>% 
        mutate(frac_coding = nbases / ncoding_bases)
}


## Han Chinese Trio Mendelian Inconsistent #####################################
load_trio_vcf <- function(trio_vcffile){
    require(VariantAnnotation)
    require(BSgenome.Hsapiens.1000genomes.hs37d5)
    grch37 <- "BSgenome.Hsapiens.1000genomes.hs37d5" 
    
    readVcf(trio_vcffile,genome = grch37)
}

get_trio_inconsistent_df <- function(trioincon_vcf){
    ## Annotating Variant Type
    geno(trioincon_vcf)[["GT"]] %>% 
        data.frame(stringsAsFactors = FALSE) %>% 
        rownames_to_column(var = "vcf_rowname") %>% 
        add_column(indel = isIndel(trioincon_vcf, 
                                   singleAltOnly = FALSE)) %>% 
        add_column(snp = isSNV(trioincon_vcf, 
                               singleAltOnly = FALSE)) %>% 
        add_column(substitution = isSubstitution(trioincon_vcf, 
                                                 singleAltOnly = FALSE))
}

get_denovo_df <- function(trioincon_vcf, trio_incon_df){
    geno(trioincon_vcf)[["GT"]] %>% 
        data.frame(stringsAsFactors = FALSE) %>% 
        rownames_to_column(var = "vcf_rowname") %>% 
        left_join(trio_incon_df) %>% 
        filter(dad == "0/0", mom == "0/0", child != "1|1") 
}

## Bechmarking ##########################################################
load_benchmarking_results <- function(benchmark_dir){
    bench_dirs <- list.dirs(benchmark_dir, recursive = TRUE) %>%
        grep(pattern = "HG.*results",value = TRUE) %>% 
        paste0("/result_1")
    
    bench_dirs <- bench_dirs %>% set_names(str_extract(.,"HG00."))
    
    ## generating a list with benchmarking results
    bench_dirs %>% map(read_happy)
    
}

get_bench_summary_df <- function(benchmark_happy_dat){
    benchmark_happy_dat$extended %>% 
        filter(Type %in% c("INDEL","SNP"), 
               Subset %in% c("*", "notinalldifficultregions"),
               Subtype == "*",
               Filter == "PASS") %>% 
        dplyr::select(Type, Subset, METRIC.Recall, 
                      METRIC.Precision, METRIC.Frac_NA)
} 
