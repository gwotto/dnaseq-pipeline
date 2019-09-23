# The GOSgene DNA-seq processing pipeline

(test master)

For a tutorial how to use this software, please consult the manual:
dnaseq-pipeline-manual.pdf

## TODOs

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
  
* Option to have genome reference files and vcf files for variant
  calibration in different directories
  
* checking if programs used are in path

* configuration
  * bed file location
  * scratch directory True/False
  * modules True/False
  * paired end true false

* reference
	
* object serialization and persistence of config file

* using cram instead of bam

* x and y chromosome

* option --annotation-group for GenotypeGVCFs

## Questions/Discussion points

The pipeline follows the standardisation effort of https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md

* picard MarkDuplicates is used instead of sambamba, because picard shows a slightly different behaviour (is this relevant to us?) (reference?)

* why are the bam files sorted twice (sambamba sort -> Picard MarkDuplicates -> sambamba sort)?

* bwa use of the hidden -K option to seed seed for deterministic mapping

## Bugs

* avoid race condition when executing makedirs in multiple threads, see
  https://stackoverflow.com/questions/42544885/error-when-mkdir-in-multi-threads-in-python
  and
  http://deepix.github.io/2017/02/02/eexists.html
