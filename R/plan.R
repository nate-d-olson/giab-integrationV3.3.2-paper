## Workflow plan
library(drake)
file.exists("inputs.csv")
file.exists("bash/get_vcf_chrom_stats.sh")
file.exists("bash/get_highhigh_vars.sh")

if(!dir.exists("data_workflow")) dir.create("data_workflow")

plan <- drake_plan(
    ## data frame with file sources
    input = read_csv(file_in("inputs.csv")),
    
    ## high confidence region coverage -----------------------------------------
    ## get non-N genome and chromosome sizes 
    chrom_size_list = list(GRCh37 = get_chrom_sizes(genome = "hs37d5"),
                           GRCh38 = get_chrom_sizes(genome = "GRCh38")),
    
    ## get coverage by chromosome
    hc_beds = get_file_list(input, type = "hc_bed"),

    hcr_cov_df = hc_beds %>%
        map_dfr(get_hc_cov,
                chrom_lengths_list = chrom_size_list,
                .id = "hgref"),

    ## high confidence variants per genome
    hc_vcfs = get_file_list(input, type = "hc_vcf"),

    hc_vcf_df = hc_vcfs %>%
        map_dfr(get_vcf_stats, vcf_type = "hc", .id = "hgref"),

    ## high confidence variants in high confidence region per genome
    hh_vcf_df = map2_dfr(hc_vcfs, hc_beds, get_hh_stats_df, .id = "hgref"),

    ## Refseq coding coverage file downloaded from
    ## https://github.com/ga4gh/benchmarking-tools/blob/master/resources/stratification-bed-files/FunctionalRegions/refseq_union_cds.sort.coordsonly.bed.gz
    refseq_coding_bed_list = list(grch37 = file_in("data_ref/refseq_union_cds.sort.coordsonly.bed.gz"),
                             grch38 = file_in("data_ref/GRCh38_refseq_cds_merged.bed.gz")),
    
    hc_coding_df = hc_beds %>% 
        map_dfr(get_coding_seq_cov, refseq_coding_bed_list, .id = "hgref"),
    
    # ## Chinese trio mendelian analysis -----------------------------------------
    ## GRCh37
    hc_trioincon_37_vcf = file_in("data_hc/HG005_HG006_HG007_trioinconsistent.vcf.gz") %>%
        load_trio_vcf(),
    hc_trioincon_37_df = get_trio_inconsistent_df(hc_trioincon_37_vcf),
    hc_denovo_37_df = get_denovo_df(hc_trioincon_37_vcf, hc_trioincon_37_df),

    # GRCh38
    hc_trioincon_38_vcf = file_in("data_hc/HG005_HG006_HG007_GRCh38_trioinconsistent.vcf.gz") %>%
        load_trio_vcf(),
    hc_trioincon_38_df = get_trio_inconsistent_df(hc_trioincon_38_vcf),
    hc_denovo_38_df = get_denovo_df(hc_trioincon_38_vcf, hc_trioincon_38_df),
    
    ## Benchmarking Results
    benchmark_dat = load_benchmarking_results("data_benchmarking"),

    ## Summary table benchmarking data
    benchmark_summary_tbl = map_dfr(benchmark_dat,
                                    get_bench_summary_df,
                                    .id = "hgref") #,
    
    ## Additional Drake parameters
    # strings_in_dots = "literals"
)
