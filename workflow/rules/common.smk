import pandas as pd

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

# function to check config files for inclusion of optional workflow steps
def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))

def get_fq_merge(wildcards):
	unit = units.loc[wildcards.sample]
	if all(pd.isna(unit["fq1"])):
		accession = unit["sra"]
		if all(unit["read_format"] == "SE"):
			return expand(
				"data/sra/se/{accession}.fastq.gz", accession=accession)
		else:
			return expand(
				"data/sra/pe/{accession}_{read}.fastq.gz", accession=accession, read = wildcards.read[-1])
	if all(unit["read_format"] == "SE"):
		return units.loc[wildcards.sample, "fq1"].tolist()
	fq = "fq{}".format(wildcards.read[-1])
	return units.loc[wildcards.sample, fq].tolist()


def get_hisat2_input(wildcards):
	if not is_activated("mergeReads"):
		unit = units.loc[wildcards.sample]
		if all(pd.isna(unit["fq1"])):
			# SRA sample (always paired-end for now)
			accession = unit["sra"]
			if all(unit["read_format"] == "SE"):
				return expand("data/sra/se/{accession}.fastq.gz", accession=accession)
			else:
				return expand("data/sra/pe/{accession}_{read}.fastq.gz", accession=accession, read=[1,2])
		fastqs = units.loc[(wildcards.sample), ["fq1", "fq2"]]
		if len(fastqs) == 2:
			return [fastqs.fq1, fastqs.fq2]
		return fastqs.fq1
	unit = units.loc[wildcards.sample]
	if all(unit["read_format"] == "SE"):
		return ["data/merged/{sample}_single.fastq.gz"]
	return ["data/merged/{sample}_1.fastq.gz", "data/merged/{sample}_2.fastq.gz"]

def get_bam_merge(wildcards):
	unit =  units[units["sample_group"] == wildcards.sample_group]
	group = pd.unique(unit["sample_name"])
	return expand(
		"results/aligned_reads/filtered/{group}.bam", group=group)



def get_scaling_input(wildcards):
	stat_files = expand(
				["results/aligned_reads/stats/{sample}_unireads.idxstats"],
				sample = units["sample_name"]
			)
	return stat_files

def get_ind_spikeIn_input(wildcards):
	unit=units.loc[wildcards.sample]
	if all(unit["call_peaks"]):
		return "results/bigwigs/zscore_normalized/individual/{sample}.bw".format(sample = wildcards.sample)

def get_merged_spikeIn_input(wildcards):
	unit =  units[units["sample_group"] == wildcards.sample]
	if all(unit["call_peaks"]):
		return "results/bigwigs/zscore_normalized/merged/{sample}.bw".format(sample = wildcards.sample)


def get_featurecounts_input(wildcards):
	sample =  samples[samples["experiment"] == wildcards.experiment]
	in_samples = pd.unique(sample["sample_name"])
	return expand(
		"results/aligned_reads/filtered/{sample}.bam", sample=in_samples)

def get_contrast(wildcards):
	return config["diff_exp"]["experiments"][wildcards.experiment]["contrasts"][wildcards.contrast]

def get_final_output():
	final_output = []
	
		# z-score normalized bigwigs for individual replicates
	final_output.extend(expand(
					[
						"results/bigwigs/zscore_normalized/individual/{sample}.bw"
					],
					sample = units["sample_name"]
				)
			)

	
	# z-score normalized bigwigs for merged replicates
	final_output.extend(expand(
					[
						"results/bigwigs/zscore_normalized/merged/{sample}.bw"
					],
					sample = units["sample_group"]
				)
			)

	if config["use_spikeIn"]:
		# spikeIn-normalized bigwigs for individual replicates
		final_output.extend(expand(
						[
							"results/bigwigs/spikeIn_normalized/individual/{sample}.bw"
						],
						sample = units.loc[units["sample_name"]]
					)
				)

		# spikeIn-normalized bigwigs for merged replicates
		final_output.extend(expand(
						[
							"results/bigwigs/spikeIn_normalized/merged/{sample}.bw"
						],
						sample = units.loc[units["sample_group"]]
					)
				)
	# count_tables
	final_output.extend(expand(
							[
								"results/count_tables/{experiment}.featureCounts"
							],
							experiment = pd.unique(samples["experiment"])
						)
					)
	# RPKM tables
	final_output.extend(expand(
							[
								"results/count_tables/{experiment}_RPKM.tsv"
							],
							experiment = pd.unique(samples["experiment"])
						)
					)
	# DEseq results
	if config["run_diff_exp"]:
		experiments = pd.unique(samples["experiment"])
		for e in experiments:
			final_output.extend(expand(
							[
								"results/DEseq2/{experiment}_{contrast}_results.tsv"
							],
							experiment = e, contrast = config["diff_exp"]["experiments"][e]["contrasts"]
						)
					)
	return final_output
