import sys
import os

bindir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(bindir, "../lib/"))

import dnaseq
import dnaseq.dnaseq
import dnaseq.ref


reference_dir = '/home/gwo/devel/dnaseq-pipeline/data/reference'

user = os.environ['USER']

sample = "sample-1"

scratchdir = os.path.join('/scratch0', user, sample)

reference_scratchdir = os.path.join(scratchdir, 'reference')

fasta = 'chr19.fa'
sequence_index = 'chr19.fa*'

known_sites = ['dbsnp_138.hg38.chr19.vcf.gz', 'Homo_sapiens_assembly38.known_indels.chr19.vcf.gz', 'Mills_and_1000G_gold_standard.indels.hg38.chr19.vcf.gz', 'hapmap_3.3.hg38.chr19.vcf.gz', '1000G_omni2.5.hg38.chr19.vcf.gz', '1000G_phase1.snps.high_confidence.hg38.chr19.vcf.gz']

known_sites_index = ['dbsnp_138.hg38.chr19.vcf.gz.tbi', 'Homo_sapiens_assembly38.known_indels.chr19.vcf.gz.tbi', 'Mills_and_1000G_gold_standard.indels.hg38.chr19.vcf.gz.tbi', 'hapmap_3.3.hg38.chr19.vcf.gz.tbi', '1000G_omni2.5.hg38.chr19.vcf.gz.tbi', '1000G_phase1.snps.high_confidence.hg38.chr19.vcf.gz.tbi']

reference_files = [fasta, sequence_index] + known_sites + known_sites_index

ref1 = dnaseq.ref.Ref(reference_dir = reference_dir,
                     reference_files = reference_files)

ref1.copy(target_dir = reference_scratchdir)
