import sys
import os.path

from pprint import pprint

bindir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(bindir, "../lib/"))

import dnaseq
import vcf

print("testing version " + dnaseq.__version__())

## get module environment working
execfile('/usr/share/modules/init/python.py')

print("removing loaded modules....")
module('purge')

module_list = ['gatk/4.0.6', 'picard-tools/2.18.5']

print("now loading modules....")

for mod in module_list:
    module('load', mod)

user = os.environ['USER']

outdir = '/home/gwo/devel/dnaseq-pipeline/test/test-outdir'

genotype_outdir = os.path.join(outdir, 'vcf-genotype')
calibrated_outdir = os.path.join(outdir, 'vcf-calibrated')

project = 'cohort-1'

vcf_dir = '/home/gwo/devel/dnaseq-pipeline/test/test-outdir/gvcf-raw/'

vcf_dict = {}

vcf_dict['sample_1'] = 'sample-1.vcf.gz'

vcf_dict['sample_2'] = 'sample-2.vcf.gz'

vcf_obj = vcf.Vcf(project = project, vcf_dir = vcf_dir,
                  vcf_dict = vcf_dict)

##, tbi_dict = tbi_dict)

vcf_obj = vcf_obj.index()

scratchdir = os.path.join('/scratch0', user, project)

reference_dir = '/home/gwo/devel/dnaseq-pipeline/data/reference'

fasta = 'chr19.fa'

bed = 'chr19.bed'

interval_padding = 0

## indices are passed as a list, maybe better pass explicit list instead of globbing
sequence_index = ['chr19.fa.*', 'chr19.dict', 'chr19.bed']


## filter_string = ''
filter_string = '--filter-expression "ExcessHet > 54.69" --filter-name ExcessHet'

vcf_obj = vcf_obj.runGenotype(outdir = genotype_outdir,
                              reference_dir = reference_dir,
                              fasta = fasta, sequence_index= sequence_index,
                              intervals_file = bed,
                              interval_padding = interval_padding,
                              filter_string = filter_string,
                              scratchdir = scratchdir)

known_sites = ['dbsnp_138.hg38.chr19.vcf.gz', 'Homo_sapiens_assembly38.known_indels.chr19.vcf.gz', 'Mills_and_1000G_gold_standard.indels.hg38.chr19.vcf.gz', 'hapmap_3.3.hg38.chr19.vcf.gz', '1000G_omni2.5.hg38.chr19.vcf.gz', '1000G_phase1.snps.high_confidence.hg38.chr19.vcf.gz', 'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.chr19.vcf.gz']

known_sites_index = ['dbsnp_138.hg38.chr19.vcf.gz.tbi', 'Homo_sapiens_assembly38.known_indels.chr19.vcf.gz.tbi', 'Mills_and_1000G_gold_standard.indels.hg38.chr19.vcf.gz.tbi', 'hapmap_3.3.hg38.chr19.vcf.gz.tbi', '1000G_omni2.5.hg38.chr19.vcf.gz.tbi', '1000G_phase1.snps.high_confidence.hg38.chr19.vcf.gz.tbi', 'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.chr19.vcf.gz.tbi']

## {refdir} interpolates to the actual reference directory
snp_resource_string = ' --resource hapmap,known=false,training=true,truth=true,prior=15.0:{refdir}/hapmap_3.3.hg38.chr19.vcf.gz --resource omni,known=false,training=true,truth=false,prior=12.0:{refdir}/1000G_omni2.5.hg38.chr19.vcf.gz --resource 1000G,known=false,training=true,truth=false,prior=10.0:{refdir}/1000G_phase1.snps.high_confidence.hg38.chr19.vcf.gz --resource dbsnp,known=true,training=false,truth=false,prior=2.0:{refdir}/dbsnp_138.hg38.chr19.vcf.gz'

indel_resource_string = ' --resource mills,known=false,training=true,truth=true,prior=12:{refdir}/Mills_and_1000G_gold_standard.indels.hg38.chr19.vcf.gz --resource dbsnp,known=true,training=false,truth=false,prior=2.0:{refdir}/dbsnp_138.hg38.chr19.vcf.gz --resource axiomPoly,known=false,training=true,truth=false,prior=10:{refdir}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.chr19.vcf.gz'


vcf_obj = vcf_obj.runCalibrate(outdir = calibrated_outdir,
                               reference_dir = reference_dir,
                               fasta = fasta, sequence_index = sequence_index,
                               known_sites = known_sites,
                               known_sites_index = known_sites_index,
                               intervals_file = '', interval_padding = '',
                               snp_resource_string = snp_resource_string,
                               indel_resource_string = indel_resource_string,
                               scratchdir = scratchdir)


pprint(vcf_obj.vcf_dict)

pprint(vcf_obj.tbi_dict)

