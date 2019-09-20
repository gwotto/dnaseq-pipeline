import sys
import os.path
import pickle

from pprint import pprint

## bindir = os.path.dirname(os.path.realpath(__file__))

bin_dir = os.path.dirname(os.path.realpath(__file__))

package_dir = os.path.normpath(os.path.join(bin_dir, '..'))

lib_dir = os.path.join(package_dir, "lib")

sys.path.append(lib_dir)

import dnaseq
import fastq
import bam

## get module environment working
## 
execfile('/usr/share/Modules/init/python.py')

print("testing version " + dnaseq.__version__())

print("removing loaded modules....")
module('purge')

module_list = ['fastqc/0.11.4', 'samtools/1.8', 'bwa/0.7.15', 'picard-tools/2.18.5', 'sambamba/0.5.9', 'gatk/4.0.6', 'samblaster/0.1.24']

print("now loading modules....")

for mod in module_list:
    module('load', mod)

user = os.environ['USER']

sample = 'sample-1'


fastq_dir =  os.path.join(package_dir,  'data/fastq/')

outdir = os.path.join(package_dir, 'test/test-outdir')

fastqc_outdir = os.path.join(outdir, 'fastqc')
mapper_outdir = os.path.join(outdir, 'bwa-raw')

scratchdir = os.path.join('/scratch0', user, sample)

reference_dir = os.path.join(package_dir, '/data/reference')
fasta = 'chr19.fa'

## indices are passed as a list, maybe better pass explicit list instead of globbing
sequence_index = ['chr19.fa*', 'chr19.dict']

fastq_dict = {}
fastq_dict['replicate_1'] = {}
fastq_dict['replicate_1']['fastq_1'] = 'sample-1_1_R1.fastq.gz'
fastq_dict['replicate_1']['fastq_2'] = 'sample-1_1_R2.fastq.gz'

fastq_dict['replicate_2'] = {}
fastq_dict['replicate_2']['fastq_1'] = 'sample-1_2_R1.fastq.gz'
fastq_dict['replicate_2']['fastq_2'] = 'sample-1_2_R2.fastq.gz'

fastq_dict['replicate_3'] = {}
fastq_dict['replicate_3']['fastq_1'] = 'sample-1_3_R1.fastq.gz'
fastq_dict['replicate_3']['fastq_2'] = 'sample-1_3_R2.fastq.gz'

print(fastq_dict)

fastq_obj = fastq.Fastq(sample = sample, fastq_dict = fastq_dict,
                        fastq_dir = fastq_dir)

print("sample 1: " + fastq_obj.sample)

print(fastq_obj.fastq_dict)

fastq_obj.runFastqc(fastqc_outdir = fastqc_outdir, scratchdir = scratchdir)

bwa_options = '-K 100000000 -Y -t 4'

b_obj = fastq_obj.runAlign(mapper_outdir = mapper_outdir,
                           reference_dir = reference_dir,
                           fasta = fasta, sequence_index = sequence_index,
                           mapper_options = bwa_options, scratchdir = scratchdir)

processed_outdir = os.path.join(outdir, 'bam-processed')

# sambamba_options = '-n'

b_obj = b_obj.runProcess(bam_outdir = processed_outdir, scratchdir = scratchdir,
                         reference_dir = reference_dir, fasta = fasta)

pprint(b_obj.bam_dict)

pprint(b_obj.index_dict)

#### calibration

calibrated_outdir = '/home/gwo/devel/dnaseq-pipeline/test/test-outdir/bam-calibrated'

known_sites = ['dbsnp_138.hg38.chr19.vcf.gz', 'Homo_sapiens_assembly38.known_indels.chr19.vcf.gz', 'Mills_and_1000G_gold_standard.indels.hg38.chr19.vcf.gz', 'hapmap_3.3.hg38.chr19.vcf.gz', '1000G_omni2.5.hg38.chr19.vcf.gz', '1000G_phase1.snps.high_confidence.hg38.chr19.vcf.gz']

known_sites_index = ['dbsnp_138.hg38.chr19.vcf.gz.tbi', 'Homo_sapiens_assembly38.known_indels.chr19.vcf.gz.tbi', 'Mills_and_1000G_gold_standard.indels.hg38.chr19.vcf.gz.tbi', 'hapmap_3.3.hg38.chr19.vcf.gz.tbi', '1000G_omni2.5.hg38.chr19.vcf.gz.tbi', '1000G_phase1.snps.high_confidence.hg38.chr19.vcf.gz.tbi']

b_out = b_obj.runCalibrate(outdir = calibrated_outdir,
                           reference_dir = reference_dir,
                           fasta = fasta, sequence_index = sequence_index,
                           intervals_file = '', interval_padding = '',
                           known_sites = known_sites,
                           known_sites_index = known_sites_index,
                           scratchdir = scratchdir)

vcf_outdir = os.path.join(outdir, 'gvcf-raw')

vcf_obj = b_out.runHaplotypeCaller(outdir = vcf_outdir,
                                   reference_dir = reference_dir,
                                   fasta = fasta, sequence_index = sequence_index,
                                   intervals_file = '', interval_padding = '',
                                   scratchdir = scratchdir)

pprint(vcf_obj.vcf_dict)
