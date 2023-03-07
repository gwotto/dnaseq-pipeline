# The GOSgene DNA-seq processing pipeline #


For a tutorial how to use this software, please consult the manual:
dnaseq-pipeline-manual.pdf

## Software dependencies

Needs the following software installed:

- BWA
- SAMtools
- Picard
- Sambamba
- samblaster
- Genome Analysis Toolkit (GATK)
- R with packages
    - ggplot2
	- gsalib
- Java >= 1.8
- python3

## Questions/Discussion points ##

* Potential conflicts can arise when there are multiple environment
  modules with the same name in different locations. To make sure to
  load right environment modules on the run node, set $MODULEPATH in
  the login node, and use qsub with the -V option, to pass the
  environment to the run nodes.

* The alignment part of the pipeline follows the standardisation
  effort of
  https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md

* picard MarkDuplicates is used instead of sambamba, because picard
  shows a slightly different behaviour (is this relevant to us?)
  (reference?)

* why are the bam files sorted twice (sambamba sort -> Picard
  MarkDuplicates -> sambamba sort)? first time is name sort, second
  time is query sort

* bwa use of the hidden -K option to seed seed for deterministic
  mapping

* The variant calling part of the pipeline follows GATK best
  practices:
  https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf.wdl

* output of HaplotypeCaller is GVCF format

* interval lists are optional input, used in base calibration, variant
  calling, combining and genotyping vcfs, but not used in variant
  recalibration.

* vcf calibration with chr19 test data gives error

	A USER ERROR has occurred: Bad input: Values for MQRankSum
    annotation not detected for ANY training variant in the input
    callset. VariantAnnotator may be used to add these annotations.

  Is this because the sample size is small and the samples are small?


## TODOs ##

* software in conda environment

* Include QC step, e.g. picard CollectMultipleMetrics

* optimise for scratch directory use, scratch directory use as the default

* Output files with discordant and split reads. This is needed by the
  lumpy SV detection software and can be achieved by samblaster in two
  ways. See also the [samblaster github
  page](https://github.com/GregoryFaust/samblaster).

  
  1. During the original alignment step
  
  ```sh
    bwa mem <idxbase> samp.r1.fq samp.r2.fq | samblaster -e -d samp.disc.sam -s samp.split.sam | samtools view -Sb - > samp.out.bam
	
  ```

  2. Pulling from a pre-existing BAM-file

  ```sh
    samtools view -h samp.bam | samblaster -a -e -d samp.disc.sam -s samp.split.sam -o /dev/null
  ```
    
* checking if programs used are in path

* configuration
  * bed file location
  * scratch directory True/False
  * modules True/False
  * paired end true false

* reference versioning
  to read cram files, the reference is needed, so it should be stored safely
	
* object serialization and persistence of config file

* x and y chromosome

* option --annotation-group for GenotypeGVCFs


## Bugs ##

* avoid race condition when executing makedirs in multiple threads,
  see
  https://stackoverflow.com/questions/42544885/error-when-mkdir-in-multi-threads-in-python
  and http://deepix.github.io/2017/02/02/eexists.html

* when a process is started repeatedly for the same sample and there
  is a lockfile, the error message overwrites the log file
  
* globbing of reference files in Ref.copy: only existing files are
  returned. If the pattern can not resolve to an existing file, the
  loop skips over it, without error messa

## Release Notes ##

* 0.1.5
   * picard quality metrics

* 0.1.4
  * pipeline running with python 3.6 now
	
* 0.1.3 (2020-06-18)
  * fixed bug causing crash when running without scratch dir. Now has
    separate temp directories for each sample, so processes dont
    interfere with same temp directory.

* 0.1.2
  * test on new cluster. interval list optional input for base
    calibration, variant calling, combining and genotyping vcfs.
  * version() function

* 0.1.1
  * writes date to output at the beginning of analysis steps
