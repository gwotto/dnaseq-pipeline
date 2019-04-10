import sys
import os.path
import pickle

bindir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(bindir, "../lib/"))

import dnaseq
# import dnaseq.dnaseq
import dnaseq.bam


## reference_dir = '/scratch0/sejjgot/sample-1/reference/'
reference_dir = '/home/gwo/devel/dnaseq-pipeline/data/reference'

fasta = 'chr19.fa'

## indices are passed as a list, maybe better pass explicit list instead of globbing
sequence_index = ['chr19.fa*', 'chr19.dict']

known_sites = ['dbsnp_138.hg38.chr19.vcf.gz', 'Homo_sapiens_assembly38.known_indels.chr19.vcf.gz', 'Mills_and_1000G_gold_standard.indels.hg38.chr19.vcf.gz', 'hapmap_3.3.hg38.chr19.vcf.gz', '1000G_omni2.5.hg38.chr19.vcf.gz', '1000G_phase1.snps.high_confidence.hg38.chr19.vcf.gz']

known_sites_index = ['dbsnp_138.hg38.chr19.vcf.gz.tbi', 'Homo_sapiens_assembly38.known_indels.chr19.vcf.gz.tbi', 'Mills_and_1000G_gold_standard.indels.hg38.chr19.vcf.gz.tbi', 'hapmap_3.3.hg38.chr19.vcf.gz.tbi', '1000G_omni2.5.hg38.chr19.vcf.gz.tbi', '1000G_phase1.snps.high_confidence.hg38.chr19.vcf.gz.tbi']

sample = 'sample-1'

bam_dir = '/home/gwo/devel/dnaseq-pipeline/test/test-outdir/bwa-raw/'

# bam_file = 'sample-1.bam'

bam_dict = {}

bam_dict['replicate_1'] = 'sample-1.bam'

bam_obj = dnaseq.bam.Bam(sample = sample, bam_dir = bam_dir,
                         bam_dict = bam_dict)

user = os.environ['USER']

outdir = '/home/gwo/devel/dnaseq-pipeline/test/test-outdir'

processed_outdir = os.path.join(outdir, 'bam-processed')
calibrated_outdir = os.path.join(outdir, 'bam-calibrated')

scratchdir = os.path.join('/scratch0', user, sample)

b_obj = bam_obj.runProcess(bam_outdir = processed_outdir, scratchdir = scratchdir)

from pprint import pprint

b_out = b_obj.runCalibrate(outdir = calibrated_outdir,
                           reference_dir = reference_dir,
                           fasta = fasta, sequence_index = sequence_index,
                           bed = '', interval_padding = '',
                           known_sites = known_sites,
                           known_sites_index = known_sites_index,
                           scratchdir = scratchdir)

pprint(b_obj.bam_dict)

pprint(b_obj.index_dict)
