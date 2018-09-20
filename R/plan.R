## Workflow plan
file.exists("inputs.csv")
file.exists("bash/get_vcf_chrom_stats.sh")
file.exists("bash/get_highhigh_vars.sh")

if(!dir.exists("data_workflow")) dir.create("data_workflow")

plan <- drake_plan(
    ## data frame with file sources
    input = read_csv(file_in("inputs.csv")),
    
    ## high confidence region coverage -----------------------------------------
    ## get non-N genome and chromosome sizes 
    chrom_size_df = get_chrom_sizes(),
    
    ## get coverage by chromosome
    hc_beds = get_file_list(input, type = "hc_bed"),

    hcr_cov_df = hc_beds %>%
        map_dfr(get_hc_cov,
                chrom_lengths = chrom_size_df,
                .id = "hgref"),


    ## high confidence variants per genome
    hc_vcfs = get_file_list(input, type = "hc_vcf"),

    hc_vcf_df = hc_vcfs %>%
        map_dfr(get_vcf_stats, vcf_type = "hc", .id = "hgref"),

    ## high confidence variants in high confidence region per genome
    hh_vcf_df = map2_dfr(hc_vcfs, hc_beds, get_hh_stats_df, .id = "hgref"),

    ## Chinese trio mendelian analysis 
    hc_trioincon_vcf = load_trio_vcf("data_hc/HG005_HG006_HG007_trioinconsistent.vcf.gz"),
    hc_trioincon_df = get_trio_inconsistent_df(hc_trioincon_vcf),
    hc_denovo_df = get_denovo_df(hc_trioincon_vcf, hc_trioincon_df),
    
    ## Benchmarking Results
    benchmark_dat = load_benchmarking_results("data_benchmarking"),

    ## Summary table benchmarking data
    benchmark_summary_tbl = map_dfr(benchmark_dat, 
                                    get_bench_summary_df, 
                                    .id = "hgref"),
    
    ## Additional Drake parameters
    strings_in_dots = "literals"
)
