rule feature_counts:
    input:
        samples=get_featurecounts_input,  # list of sam or bam files
        annotation="resources/annotation.gff.gz",

    output:
        multiext(
            "results/count_tables/{experiment}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 2
    params:
        r_path="",  # implicitly sets the --Rpath flag
        extra="",
    log:
        "logs/feature_counts/{sample}.log",
    wrapper:
        "v1.3.2/bio/subread/featurecounts"