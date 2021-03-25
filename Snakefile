IDS, = glob_wildcards("/icgc/dkfzlsdf/analysis/B260/projects/HipSci/openAccess/endoderm_differentation/data_raw/run_21999/fastq/21999_{cell}_1_val_1.fq.gz")

print(IDS) 
print(expand("/dkfz/groups/OE0540/users/ruehle/rnavelocity/pre_kallisto/hg38_release101/velocity_files_coll/counts/run21999/21999_{cell}", cell=IDS))

rule all: 
 input: expand("/dkfz/groups/OE0540/users/ruehle/rnavelocity/pre_kallisto/hg38_release101/velocity_files_coll/counts/run21999/21999_{cell}", cell=IDS)

rule:
 input: 
  read1 = "/icgc/dkfzlsdf/analysis/B260/projects/HipSci/openAccess/endoderm_differentation/data_raw/run_21999/fastq/21999_{cell}_1_val_1.fq.gz",
  read2 = "/icgc/dkfzlsdf/analysis/B260/projects/HipSci/openAccess/endoderm_differentation/data_raw/run_21999/fastq/21999_{cell}_2_val_2.fq.gz"
 output: directory("/dkfz/groups/OE0540/users/ruehle/rnavelocity/pre_kallisto/hg38_release101/velocity_files_coll/counts/run21999/21999_{cell}")
 shell: "kallisto quant -i /dkfz/groups/OE0540/users/ruehle/rnavelocity/pre_kallisto/hg38_release101/velocity_files_coll/hs_cDNA_introns_101.idx -o {output} {input.read1} {input.read2}"
