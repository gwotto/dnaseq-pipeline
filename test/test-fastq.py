import sys
import os.path
import pickle

bindir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(bindir, "../lib/"))

import dnaseq
# import dnaseq.dnaseq
import dnaseq.fastq

## TODO create fastq object directly, not using pickle
## I did this before, but do not remember where
# fastq_obj_file = 'sample-1.pkl'
# fastq_obj = pickle.load(open(fastq_obj_file))

sample = 'sample-1'

fastq_dir = '/home/sejjgot/devel/dnaseq-pipeline/data/fastq-2/'

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

fastq_obj = dnaseq.fastq.Fastq(sample = sample, fastq_dict = fastq_dict,
                  fastq_dir = fastq_dir)

print("sample 1: " + fastq_obj.sample)

user = os.environ['USER']
sample = fastq_obj.sample


print(fastq_obj.fastq_dict)

outdir = '/home/sejjgot/devel/dnaseq-pipeline/test/test-outdir'

fastqc_outdir = os.path.join(outdir, 'fastqc')
mapper_outdir = os.path.join(outdir, 'bwa-raw')

scratchdir = os.path.join('/scratch0', user, sample)

fastq_obj.runFastqc(fastqc_outdir = fastqc_outdir, scratchdir = scratchdir)


reference_dir = '/home/sejjgot/devel/dnaseq-pipeline/data/reference'
fasta = 'chr19.fasta'
mapper_index_base = 'chr19.fasta*'

bwa_options = '-K 100000000 -Y -t 4'

bam_obj = fastq_obj.runAlign(mapper_outdir = mapper_outdir,
                             reference_dir = reference_dir,
                             fasta = fasta, mapper_index = mapper_index_base,
                             mapper_options = bwa_options, scratchdir = scratchdir)

print(bam_obj.bam_dict)
