import pandas as pd

patients = list(pd.unique(pd.read_table("../patients.csv")["patient"]))
samples = ["normal", "tumor"]

configfile: "config.yaml"

rule all:
     input: 
     	    expand("{patient}.{sample}.mosdepth.global.dist.txt", patient=patients, sample=samples),
     	    expand("{patient}.{sample}.mosdepth.summary.txt", patient=patients, sample=samples), 
	    expand("{patient}.{sample}.mosdepth.region.dist.txt", patient=patients, sample=samples), 
	    expand("{patient}.{sample}.per-base.bed.gz", patient=patients, sample=samples),
	    expand("{patient}.{sample}.regions.bed.gz", patient=patients, sample=samples),
	    expand("{patient}.{sample}.regions.bed.gz.csi", patient=patients, sample=samples),
	    expand("{patient}.{sample}.per-base.bed.gz.csi", patient=patients, sample=samples)

def get_patients(wildcards):
    return config["patients"][wildcards.patient][wildcards.sample]


rule mosD:
     input: 
     	bed="beds/35utr.merge.bed",
	bam=get_patients
     output:
        "{patient}.{sample}.mosdepth.global.dist.txt",
	"{patient}.{sample}.mosdepth.summary.txt", 
	"{patient}.{sample}.mosdepth.region.dist.txt",
	"{patient}.{sample}.per-base.bed.gz",
	"{patient}.{sample}.regions.bed.gz",
	"{patient}.{sample}.regions.bed.gz.csi",
	"{patient}.{sample}.per-base.bed.gz.csi"
     singularity: 
     	"/home/groups/carilee/software/containers/mosdepth_latest.sif"
     shell: 
     	"""
	mosdepth --by {input.bed} {wildcards.patient}.{wildcards.sample} {input.bam}
	"""