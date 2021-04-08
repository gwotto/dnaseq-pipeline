"""
The vcf module
"""

import os.path
import shutil
import re

import dnaseq
import ref


class Vcf:
    """Vcf is a class and data structure that contains information
    about the vcf and index files of on or more samples of a project or
    cohort with the elements

    Attributes:
    -----------
    project (str):
       sample or cohort name.
    vcf_dict:
       dictionary containing the file names of vcf files (values) from
       the samples (keys, see details).
    tbi_dict:
       dictionary containing the file names of tbi index files (values)
       from the samples (keys, see details).
    vcf_dir(str)
       path to directory where vcf and index files are located
    -----------

    vcf_dict is a dictionary, with keys being the sample names of the
    cohort, and vcf file names as values. tbi_dict is a dictionary,
    with keys being the sample names of the cohort, and tbi index file
    names as values. vcf_dict and tbi_dict must have the same key
    names.

    |
    |-- sample_1 -> sample_1.vcf
    |
    |-- sample_2 -> sample_2.vcf

    """

    # The class "constructor" - It's actually an initializer 
    def __init__(self, project = None, vcf_dir = None, vcf_dict = None, tbi_dict = None):

        if project is None:
            self.project = ''
        else:
            self.project = project
            
        if vcf_dir is None:
            self.vcf_dir = ''
        else:
            self.vcf_dir = vcf_dir

        if vcf_dict is None:
            self.vcf_dict = {}
        else:
            self.vcf_dict = vcf_dict

        if tbi_dict is None:
            self.tbi_dict = {}
        else:
            self.tbi_dict = tbi_dict


    ## TODO checkFiles function
    def copy(self, target_dir):

        project = self.project
        vcf_dir = self.vcf_dir
        vcf_dict = self.vcf_dict
        tbi_dict = self.tbi_dict
        
        if not os.path.exists(target_dir):
            print("directory " + target_dir + " does not exist.")
            sys.exit()

        vcf_dict = dnaseq.dictRecurse(dictionary = vcf_dict,
                                      f = dnaseq.copyFile,
                                      source_dir = vcf_dir,
                                      target_dir = target_dir)

        tbi_dict = dnaseq.dictRecurse(dictionary = tbi_dict,
                                      f = dnaseq.copyFile,
                                      source_dir = vcf_dir,
                                      target_dir = target_dir)
        
        vcf_obj = Vcf(project = project, vcf_dir = target_dir,
                      vcf_dict = vcf_dict, tbi_dict = tbi_dict)

        return vcf_obj

    
    def index(self):

        sample = self.project
        vcf_dir = self.vcf_dir
        vcf_dict = self.vcf_dict

        tbi_dict = dnaseq.dictRecurse(dictionary = vcf_dict,
                                        f = indexVcfFile,
                                        vcf_dir = vcf_dir)

        self.tbi_dict = tbi_dict
        return self

    def vcfList(self):
        ## list of vcf file names
        vcf_list=[]

        vcf_dict = self.vcf_dict
        sample_keys = vcf_dict.keys()

        for sample in sample_keys:
            vcf_list.append(vcf_dict[sample])

        return vcf_list


    def runGenotype(self, outdir, reference_dir, fasta, sequence_index, scratchdir, intervals_file, interval_padding, filter_string):

        project = self.project
        vcf_dict = self.vcf_dict
        vcf_dir = self.vcf_dir

        print('\n***   combining vcf files   ***\n')

        vcf_outfile = project + ".vcf.gz"
        vcf_outpath = os.path.join(outdir, vcf_outfile)

        ## vcf_idx = vcf_outfile + ".tbi"
        
        if os.path.isfile(vcf_outpath):
            print("output vcf file " + vcf_outpath +
                  " exists, skipped genotyping vcf")
            
            vcf_dict =  {}
            ## TODO generate key programmatically
            vcf_dict['sample_1'] = vcf_outfile

            vcf_obj = Vcf(project = project, vcf_dict = vcf_dict, vcf_dir = outdir)

            vcf_obj = vcf_obj.index()
        
            return vcf_obj

        else:
            print("preparing to combine and genotype vcf files")
            
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            reference_files = [fasta, intervals_file] + sequence_index
            
            ref1 = ref.Ref(reference_dir = reference_dir,
                           reference_files = reference_files)

            if bool(scratchdir):

                print("Using temporary directory " + scratchdir)

                reference_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                               source_dir = reference_dir)

                ## better return a ref object that sits on scratch
                ref1.copy(target_dir = reference_scratchdir)
                
                vcf_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                         source_dir = vcf_dir)
                        
                vcf_obj = self.copy(target_dir = vcf_scratchdir)

                genotyped_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                                source_dir = outdir)
            
                vcf_obj = vcf_obj.genotype(outdir = genotyped_scratchdir,
                                           reference_dir = reference_scratchdir,
                                           fasta = fasta,
                                           intervals_file = intervals_file,
                                           interval_padding = interval_padding,
                                           filter_string = filter_string)

                print('deleting input scratch directory ' + vcf_scratchdir)
                shutil.rmtree(vcf_scratchdir)
        
                vcf_obj = vcf_obj.copy(target_dir = outdir)

                print('deleting output scratch directory ' + genotyped_scratchdir)
                shutil.rmtree(genotyped_scratchdir)

            else:

                vcf_obj = self.genotype(outdir = outdir,
                                        reference_dir = reference_dir,
                                        fasta = fasta,
                                        intervals_file = intervals_file,
                                        interval_padding = interval_padding,
                                        filter_string = filter_string)
                
        return vcf_obj



        
    def genotype(self, outdir, reference_dir, fasta,
                 intervals_file, interval_padding, filter_string):
        
        project = self.project
        vcf_dir = self.vcf_dir
        vcf_dict = self.vcf_dict
        tbi_dict = self.tbi_dict


        if not os.path.exists(outdir):
            print("directory " + outdir + " does not exist.")
            sys.exit()

        
        ## generate an input string for CombineGVCFs using the vcf file list
        vcf_list = self.vcfList()
        vcf_pathlist = [('-V ' + os.path.join(vcf_dir, i)) for i in vcf_list]
        input_string = ' '.join(vcf_pathlist)

        ## output file
        vcf_outfile = project + ".vcf.gz"
        vcf_outpath = os.path.join(outdir, vcf_outfile)


        ## the temporary file for the combined gvcs
        vcf_combined_tempdir = os.path.join(outdir, 'combined-temp')

        if not os.path.exists(vcf_combined_tempdir):
            os.makedirs(vcf_combined_tempdir)

        vcf_combined_temppath = os.path.join(vcf_combined_tempdir, vcf_outfile)

        
        ## the temporary file for the genotyped vcs
        vcf_geno_tempdir = os.path.join(outdir, 'geno-temp')

        if not os.path.exists(vcf_geno_tempdir):
            os.makedirs(vcf_geno_tempdir)

        vcf_geno_temppath = os.path.join(vcf_geno_tempdir, vcf_outfile)


        ## the temporary file for the filtered vcfs. not needed,
        ## because I skip MakeSitesOnlyVcf, so filtered vcfs go to
        ## vcf_outpath
        # vcf_filtered_tempdir = os.path.join(outdir, 'filtered-temp')

        # if not os.path.exists(vcf_filtered_tempdir):
        #     os.makedirs(vcf_filtered_tempdir)

        # vcf_filtered_temppath = os.path.join(vcf_filtered_tempdir, vcf_outfile)


        ## fasta = reference_files.pop(0)
        fasta_path = os.path.join(reference_dir, fasta)

        
        ## interval file not used in CombineGVCFs, or genotypeGVCFs but
        ## required for GenomicsDBImport
        interval_path = os.path.join(reference_dir, intervals_file)
      
        if os.path.isfile(interval_path):
            interval_option = ' -L ' + interval_path + ' -ip ' + str(interval_padding)
        else:
            interval_option = ''
            

        ## currently GenomicsDBImport runs only on one interval at a
        ## time, so I use combineGVCFs
        # db_dir = project + '_db'
        # ## TODO test: db_dirpath must not exist
        # db_dirpath = os.path.join(outdir, db_dir)
    
        # gatk_genomicsdbimport_c = 'gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport -verbosity=WARNING ' + input_string + ' --genomicsdb-workspace-path ' + db_dirpath + interval_option

        # print("running GATK GenomicsDBImport")
        # print(gatk_genomicsdbimport_c)

        # os.system(gatk_genomicsdbimport_c)
        gatk_combinegvcfs_c = 'gatk CombineGVCFs -R ' + fasta_path + ' ' + input_string + ' -O ' + vcf_combined_temppath
        
        print("running GATK CombineGVCFs")
        print(gatk_combinegvcfs_c)

        os.system(gatk_combinegvcfs_c)
        
        ## now run genotype
        gatk_genotypegvcf_c = 'gatk GenotypeGVCFs -verbosity=WARNING -R ' + fasta_path + ' -V ' + vcf_combined_temppath + ' -G StandardAnnotation --use-new-qual-calculator -O ' + vcf_geno_temppath
        
        print("running GATK GenotypeGVCFs")
        print(gatk_genotypegvcf_c)

        os.system(gatk_genotypegvcf_c)

        print('deleting temporary directory ' + vcf_combined_tempdir)
        shutil.rmtree(vcf_combined_tempdir)


        ## filter vcf

        gatk_filtervcf_c = 'gatk VariantFiltration ' + filter_string + ' -V ' + vcf_geno_temppath + ' -O ' + vcf_outpath

        print("running GATK VariantFiltration")
        print(gatk_filtervcf_c)

        os.system(gatk_filtervcf_c)

        print('deleting temporary directory ' + vcf_geno_tempdir)
        shutil.rmtree(vcf_geno_tempdir)


        ## sites only vcf (picard program) this is in the GATK best
        ## practices, but I do not use it, because I think we need the
        ## genotype information
        
        # picard_makesitesonlyvcf_c = 'MakeSitesOnlyVcf I=' + vcf_filtered_temppath + ' O=' + vcf_outpath

        # print("running picard MakeSitesOnlyVcf")
        # print(picard_makesitesonlyvcf_c)

        # os.system(picard_makesitesonlyvcf_c)

        # print('deleting temporary directory ' + vcf_filtered_tempdir)
        # shutil.rmtree(vcf_filtered_tempdir)

        vcf_out_dict = {}

        ## TODO generate key programmatically
        vcf_out_dict['sample_1'] = vcf_outfile
        
        vcf_obj = Vcf(project = project, vcf_dict = vcf_out_dict,
                      vcf_dir = outdir)

        vcf_obj = vcf_obj.index()
        
        return vcf_obj



    def runCalibrate(self, outdir, reference_dir, fasta, sequence_index,
                     known_sites, known_sites_index, snp_resource_string,
                     indel_resource_string, scratchdir):

        project = self.project
        vcf_dict = self.vcf_dict
        vcf_dir = self.vcf_dir

        print('\n***   variant calibration   ***\n')
        
        ## TODO
        ## only runs on first sample if there are multiple replicates
        vcf_key = sorted(vcf_dict)[0]

        ## TODO better construct the vcf filename from input vcf object
        vcf_outfile = project + ".vcf.gz"
        vcf_outpath = os.path.join(outdir, vcf_outfile)

        vcf_idx = vcf_outfile + ".tbi"

        
        if os.path.isfile(vcf_outpath):
            print("calibrated vcf file " + vcf_outpath +
                  " exists, skipped calibrating vcf")
            
            vcf_dict = {}
            vcf_dict[vcf_key] = vcf_outfile

            tbi_dict = {}
            tbi_dict[vcf_key] = vcf_idx

            vcf_obj = Vcf(project = project, vcf_dict = vcf_dict,
                            tbi_dict = tbi_dict, vcf_dir = outdir)

            return vcf_obj

        else:
            print("preparing to calibrate vcf file")
            
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            plot_outdir = os.path.join(outdir, 'plots')

            if not os.path.exists(plot_outdir):
                os.makedirs(plot_outdir)

            ## reference_files = [fasta, intervals_file] + sequence_index + known_sites + known_sites_index

            reference_files = [fasta] + sequence_index + known_sites + known_sites_index
            
            ref1 = ref.Ref(reference_dir = reference_dir,
                           reference_files = reference_files)

            if bool(scratchdir):

                print("Using temporary directory " + scratchdir)

                reference_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                               source_dir = reference_dir)

                ## better return a ref object that sits on scratch
                ref1.copy(target_dir = reference_scratchdir)
                
                vcf_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                         source_dir = vcf_dir)
                        
                vcf_obj = self.copy(target_dir = vcf_scratchdir)

                calibrated_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                                source_dir = outdir)

                plot_scratchdir = dnaseq.createScratchDir(scratchdir = calibrated_scratchdir,
                                                                source_dir = plot_outdir)
                
                vcf_obj = vcf_obj.calibrate(outdir = calibrated_scratchdir,
                                            plot_outdir = plot_scratchdir,
                                            reference_dir = reference_scratchdir,
                                            fasta = fasta, known_sites = known_sites,
                                            snp_resource_string = snp_resource_string,
                                            indel_resource_string = indel_resource_string)

                print('deleting input scratch directory ' + vcf_scratchdir)
                shutil.rmtree(vcf_scratchdir)
        
                vcf_obj = vcf_obj.copy(target_dir = outdir)

                ## copying the plot scratchdir
                dnaseq.copytree(plot_scratchdir, plot_outdir)
                
                print('deleting output scratch directory ' + calibrated_scratchdir)
                shutil.rmtree(calibrated_scratchdir)

            else:
                
                vcf_obj = self.calibrate(outdir = outdir,
                                         plot_outdir = plot_outdir,
                                         reference_dir = reference_dir,
                                         fasta = fasta, known_sites = known_sites,
                                         snp_resource_string = snp_resource_string,
                                         indel_resource_string = indel_resource_string)
                
        return vcf_obj

    


    def calibrate(self, outdir, plot_outdir, reference_dir, fasta, known_sites,
                  snp_resource_string, indel_resource_string):
 
        project = self.project
        vcf_dir = self.vcf_dir
        vcf_dict = self.vcf_dict
        tbi_dict = self.tbi_dict

        ## TODO
        ## only runs on first replicate if there are multiple replicates
        vcf_key = sorted(vcf_dict)[0]
        
        vcf_file = vcf_dict[vcf_key]
        vcf_path = os.path.join(vcf_dir, vcf_file)
            
        vcf_outfile = project + ".vcf.gz"
        vcf_outpath = os.path.join(outdir, vcf_outfile)
                
        if not os.path.exists(outdir):
            print("directory " + outdir + " does not exist.")
            sys.exit()

        ## fasta = reference_files.pop(0)
        fasta_path = os.path.join(reference_dir, fasta)


        ## interval not used in calibration step (see gatk documentation)
        # intervals_path = os.path.join(reference_dir, intervals_file)
        # if os.path.isfile(intervals_path):
        #     intervals_option = ' -L ' + intervals_path + ' -ip ' + str(interval_padding)
        # else:
        #     intervals_option = ''


        # ## directory for output plots
        # plot_outdir = os.path.join(outdir, 'plots')

        # if not os.path.exists(plot_outdir):
        #     os.makedirs(plot_outdir)
            
        ## temporary files and directories for snp recalibration
        snp_cal_tempdir = os.path.join(outdir, 'snp-cal-temp')

        if not os.path.exists(snp_cal_tempdir):
            os.makedirs(snp_cal_tempdir)

        snp_cal_temppath = os.path.join(snp_cal_tempdir, vcf_outfile)

        snp_tranches_file = vcf_file + '_snp-tranches'
        snp_tranches_outpath = os.path.join(plot_outdir, snp_tranches_file)

        snp_rscript_file = project + '_snp.output.plots.R'
        snp_rscript_outpath = os.path.join(plot_outdir, snp_rscript_file)

        ## interpolates {refdir} to the valye of ref_dir
        snp_resource_string = snp_resource_string.format(refdir = reference_dir)
        
        gatk_snprecalibrator_c = 'gatk VariantRecalibrator --verbosity=WARNING -R '  + fasta_path + ' -V ' + vcf_path + ' -O ' + snp_cal_temppath + ' ' + snp_resource_string + ' -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP --tranches-file ' + snp_tranches_outpath + ' --rscript-file ' +  snp_rscript_outpath
        
        print("running GATK VariantRecalibrator with SNPs")
        print(gatk_snprecalibrator_c)

        os.system(gatk_snprecalibrator_c)
            
        ## temporary files and directories for indel recalibration
        indel_cal_tempdir = os.path.join(outdir, 'indel-cal-temp')

        if not os.path.exists(indel_cal_tempdir):
            os.makedirs(indel_cal_tempdir)

        indel_cal_temppath = os.path.join(indel_cal_tempdir, vcf_outfile)

        indel_tranches_file = vcf_file + '_indel-tranches'
        indel_tranches_outpath = os.path.join(plot_outdir, indel_tranches_file)

        indel_rscript_file = project + '_indel.output.plots.R'
        indel_rscript_outpath = os.path.join(plot_outdir, indel_rscript_file)

        ## interpolates {refdir} to the valye of ref_dir
        indel_resource_string = indel_resource_string.format(refdir = reference_dir)
        
        gatk_indelrecalibrator_c = 'gatk VariantRecalibrator --verbosity=WARNING -R '  + fasta_path + ' -V ' + vcf_path + ' -O ' + indel_cal_temppath + ' ' + indel_resource_string + ' -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode INDEL --tranches-file ' + indel_tranches_outpath + ' --rscript-file ' +  indel_rscript_outpath
        
        print("running GATK VariantRecalibrator with Indels")
        print(gatk_indelrecalibrator_c)

        os.system(gatk_indelrecalibrator_c)


        ## temporary files and directories for snp recalibration
        snp_appl_tempdir = os.path.join(outdir, 'snp-appl-temp')

        if not os.path.exists(snp_appl_tempdir):
            os.makedirs(snp_appl_tempdir)

        snp_appl_temppath = os.path.join(snp_appl_tempdir, vcf_outfile)
        

        gatk_applyvqsr_snp_c = 'gatk ApplyVQSR -R ' + fasta_path + ' -V ' + vcf_path + ' -O ' + snp_appl_temppath + ' --truth-sensitivity-filter-level 99.0 --tranches-file ' + snp_tranches_outpath + ' --recal-file ' + snp_cal_temppath + ' -mode SNP'


        print("running GATK ApplyVQSR with SNPs")
        print(gatk_applyvqsr_snp_c)

        os.system(gatk_applyvqsr_snp_c)


        gatk_applyvqsr_indel_c = 'gatk ApplyVQSR -R ' + fasta_path + ' -V ' + snp_appl_temppath + ' -O ' + vcf_outpath + ' --truth-sensitivity-filter-level 99.0 --tranches-file ' + indel_tranches_outpath + ' --recal-file ' + indel_cal_temppath + ' -mode INDEL'


        print("running GATK ApplyVQSR with indels")
        print(gatk_applyvqsr_indel_c)

        os.system(gatk_applyvqsr_indel_c)
        
        print('deleting temporary directory ' + indel_cal_tempdir)
        shutil.rmtree(indel_cal_tempdir)

        print('deleting temporary directory ' + snp_cal_tempdir)
        shutil.rmtree(snp_cal_tempdir)

        vcf_obj = Vcf(project = project, vcf_dict = vcf_dict,
                      tbi_dict = tbi_dict,
                      vcf_dir = outdir)

        return vcf_obj
    

def indexVcfFile(vcf_file, vcf_dir):

    """ creates an tbi index for a vcf file

    Args:
        vcf_file (str): The name of the vcf file to be indexed
        vcf_dir (str): The path to the directory where the vcf file is
            located.

    Returns:
        The file name of the tbi index file

    """

    index_suffix = 'tbi'
    
    vcf_path = os.path.join(vcf_dir, vcf_file)
    
    index_file = vcf_file + '.' + index_suffix

    ## TODO create index if it does not exist
    index_path = os.path.join(vcf_dir, index_file)

    return index_file

