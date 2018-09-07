## Workflow plan
file.exists("inputs.csv")
file.exists("bash/get_vcf_chrom_stats.sh")
file.exists("bash/get_highhigh_vars.sh")

if(!dir.exists("workflow_data")) dir.create("workflow_data")

plan <- drake_plan(
    ## data frame with file sources
    input = read_csv("inputs.csv"),
    
    ## high confidence region coverage -----------------------------------------
    ## get non-N genome and chromosome sizes 
    # Error message - could not find function "as.workflow_data.frame"
    # chrom_size_df = get_chrom_sizes(BSgenome.Hsapiens.1000genomes.hs37d5),
    
    ## get coverage by chromosome
    hc_beds = get_file_list(input, type = "hc_bed"),
    
    hcr_cov_df = hc_beds %>% 
        map_dfr(get_hc_cov, bed_source = file_source, .id = "hgref"),
 

    ## high confidence variants per genome
    hc_vcfs = get_file_list(input, type = "hc_vcf"),
    
    hc_vcf_df = hc_vcfs %>% 
        map_dfr(get_vcf_stats, vcf_type = "hc", .id = "hgref"),
    
    ## high confidence variants in high confidence region per genome
    hh_vcf_df = map2_dfr(hc_vcfs, hc_beds, get_hh_stats_df, .id = "hgref"),
    strings_in_dots = "literals"
)
