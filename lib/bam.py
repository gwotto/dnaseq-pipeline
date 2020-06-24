"""
The bam module
"""

import os.path
import shutil
import filecmp
import re
import errno

import dnaseq

import vcf
import ref


## TODO checkFiles function

class Bam:
    """Bam is a class and data structure that contains information
    about the bam and index files of a sample with the elements

    Attributes:
    -----------
    sample (str):
       sample name.
    bam_dict:
       dictionary containing the file names of bam files from the 
       sample (see details).
    index_dict:
       dictionary containing the file names of bam index files from
       the sample (see details).
    bam_dir(str)
       path to directory where bam files are located
    -----------

    bam_dict is a dictionary containing the names of bam files from
    the sample. There can be multiple bam files from one
    sample. index_dict is the dictionary containing the filenames of
    the bam index files.  bam_dict and index_dict have the same key
    names.


    |
    |-- bam_key_1
    |
    |-- bam_key_2

    """

    # The class "constructor" - It's actually an initializer 
    def __init__(self, sample = None, bam_dir = None, bam_dict = None, index_dict = None):

        if sample is None:
            self.sample = ''
        else:
            self.sample = sample
            
        if bam_dict is None:
            self.bam_dict = {}
        else:
            self.bam_dict = bam_dict

        if index_dict is None:
            self.index_dict = {}
        else:
            self.index_dict = index_dict

        if bam_dir is None:
            self.bam_dir = ''
        else:
            self.bam_dir = bam_dir


    ## not sure if this is needed    
    def bamList(self):
        ## list of bam file names
        bam_list=[]

        bam_dict = self.bam_dict
        rep_keys = bam_dict.keys()

        for rep in rep_keys:
            ## for fq in fastq_dict[rep].keys():
            bam_list.append(bam_dict[rep])##['bam-file'])

        return bam_list


    def merge(self, bam_outdir):
 
        sample = self.sample
        bam_dir = self.bam_dir
        bam_dict = self.bam_dict

        bam_outfile = sample + ".bam"
        bam_outpath = os.path.join(bam_outdir, bam_outfile)
        
        if not os.path.exists(bam_outdir):
            print("directory " + bam_outdir + " does not exist.")
            sys.exit()

        ## only rums on one (first) key of bam dictionary if there are
        ## multiple keys
        bam_key_list = bam_dict.keys()
        bam_key_list.sort()
        bam_key = bam_key_list[0]
           
        bam_list = self.bamList()

        bam_path_list = [bam_dir + "/" + bam for bam in bam_list]

        bam_string = ' '.join(bam_path_list)

        samtools_merge_c = "samtools merge " + bam_outpath + " " + bam_string

        print("running samtools merge")
        print(samtools_merge_c)

        os.system(samtools_merge_c)

        bam_dict = {}
        bam_dict[bam_key] = bam_outfile

        ## todo include bam indices
        b_obj = Bam(sample = sample, bam_dict = bam_dict,
                              bam_dir = bam_outdir)
        
        return b_obj

    
    def index(self):

        sample = self.sample
        bam_dir = self.bam_dir
        bam_dict = self.bam_dict

        index_dict = dnaseq.dictRecurse(dictionary = bam_dict,
                                        f = indexBamFile,
                                        bam_dir = bam_dir)

        self.index_dict = index_dict
        return self

    
    def convert2Cram(self, fasta):

        sample = self.sample
        bam_dir = self.bam_dir
        bam_dict = self.bam_dict

        bam_dict = dnaseq.dictRecurse(dictionary = bam_dict,
                                      f = bam2Cram,
                                      bam_dir = bam_dir,
                                      fasta = fasta)

        index_dict = {}

        c_obj = Bam(sample = sample, bam_dir = bam_dir,
                    bam_dict = bam_dict, index_dict = index_dict)

        return c_obj
    

    
    def copy(self, target_dir):

        sample = self.sample
        bam_dir = self.bam_dir
        bam_dict = self.bam_dict
        index_dict = self.index_dict
        
        if not os.path.exists(target_dir):
            print("directory " + target_dir + " does not exist.")
            sys.exit()

        bam_dict = dnaseq.dictRecurse(dictionary = bam_dict, f = dnaseq.copyFile,
                                      source_dir = bam_dir, target_dir = target_dir)

        index_dict = dnaseq.dictRecurse(dictionary = index_dict,
                                        f = dnaseq.copyFile,
                                        source_dir = bam_dir,
                                        target_dir = target_dir)
        
        b_obj = Bam(sample = sample, bam_dir = target_dir,
                    bam_dict = bam_dict, index_dict = index_dict)

        return b_obj    

    
    def delete(self):

        sample = self.sample
        bam_dir = self.bam_dir
        bam_dict = self.bam_dict
        index_dict = self.index_dict
        
        # if not os.path.exists(target_dir):
        #     print("directory " + target_dir + " does not exist.")
        #     sys.exit()

        bam_dict = dnaseq.dictRecurse(dictionary = bam_dict, f = dnaseq.deleteFile,
                                      dir = bam_dir)

        index_dict = dnaseq.dictRecurse(dictionary = index_dict,
                                        f = dnaseq.deleteFile,
                                        dir = bam_dir)

        ## empty the bama and index dictionaries
        self.bam_dict = {}

        self.index_dict = {}

        return self



    def runProcess(self, bam_outdir, reference_dir, fasta, scratchdir):


        ## TODO adjusting for object with multiple bam files
        
        ## 1. MC and MQ tags are added using samblaster (v0.1.24) with
        ## the parameters `-a --addMateTags`. (fastq.align())

        ## 2. Read group BAM files are merged together with `samtools
        ## merge` (v1.3.1-2). (fastq.align())

        ## 3. The resulting file is name-sorted with `sambamba sort -n`
        ## (v0.6.4).
        
        ## 4. Duplicates are marked using Picard MarkDuplicates (v2.4.1)
        ## with the parameter `ASSUME_SORT_ORDER=queryname`

        ## 5. then the results are coordinate sorted using `sambamba
        ## sort`.

        
        sample = self.sample
        bam_dir = self.bam_dir
        bam_dict = self.bam_dict

        print('\n***   bam file processing   ***\n')
        
        ## only rums on one (first) key of bam dictionary if there are
        ## multiple keys
        bam_key_list = bam_dict.keys()
        bam_key_list.sort()
        bam_key = bam_key_list[0]

        
        bam_outfile = sample + ".bam"
        bam_outpath = os.path.join(bam_outdir, bam_outfile)

        cram_outfile = sample + ".cram"
        cram_outpath = os.path.join(bam_outdir, cram_outfile)

        if os.path.isfile(cram_outpath):
            print("processed output file " + cram_outpath +
                  " exists, skipped processing bam")

            ## createing the output bam object with crm files
            bam_dict = {}
            bam_dict[bam_key] = cram_outfile

            ## create index_dict recursively, assuming the index files
            ## exist and have the form bam_file + ".bai"
            index_dict = dnaseq.dictRecurse(dictionary = bam_dict,
                                        f = lambda x: x + '.crai')
            
            b_obj = Bam(sample = sample, bam_dict = bam_dict,
                        index_dict = index_dict, bam_dir = bam_outdir)

        else:
            print("preparing to process bam file")
            
            if not os.path.exists(bam_outdir):
                os.makedirs(bam_outdir)

            reference_files = [fasta]
            
            ref1 = ref.Ref(reference_dir = reference_dir,
                           reference_files = reference_files)

            if bool(scratchdir):

                print("Using temporary directory " + scratchdir)
                                
                reference_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                               source_dir = reference_dir)
            

                ref1.copy(target_dir = reference_scratchdir)
                
                
                bam_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                         source_dir = bam_dir)

                b_obj = self.copy(target_dir = bam_scratchdir)

        
                processed_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                               source_dir = bam_outdir)

                b_obj = b_obj.process(bam_outdir = processed_scratchdir,
                                      reference_dir = reference_dir, fasta = fasta)

                ## currently, the input bam file is deleted by the process function
                print('deleting input scratch directory ' + bam_scratchdir)
                shutil.rmtree(bam_scratchdir)

                b_obj = b_obj.copy(target_dir = bam_outdir)

                print('deleting output scratch directory ' + processed_scratchdir)
                shutil.rmtree(processed_scratchdir)

            else:

                b_obj = self.process(bam_outdir = bam_outdir,
                                     reference_dir = reference_dir, fasta = fasta)
                
        return b_obj


    def process(self, bam_outdir, reference_dir, fasta):
 
        sample = self.sample
        bam_dir = self.bam_dir
        bam_dict = self.bam_dict

        # only processes the first bam file if there are multiple
        bam_key_list = bam_dict.keys()
        bam_key_list.sort()
        bam_key = bam_key_list[0]
        
        bam_file = bam_dict[bam_key]
        bam_path = os.path.join(bam_dir, bam_file)
        
        bam_outfile = sample + '.bam'
        bam_outpath = os.path.join(bam_outdir, bam_outfile)
        
        bam_tempdir = os.path.join(bam_outdir, (sample + '-temp'))
        ## tmp dir for MarkDuplicates
        tmp_dir = os.path.join(bam_tempdir, 'tmp')
        
        fasta_path = os.path.join(reference_dir, fasta)

        if not os.path.exists(bam_outdir):
            print("directory " + bam_outdir + " does not exist.")
            sys.exit()
            
        if not os.path.exists(bam_tempdir):
            os.makedirs(bam_tempdir)

        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)


        bam_temppath = os.path.join(bam_tempdir, sample + "_sambamba_1.bam")

        ## increasing the memory limit seems to prevent the 'Too many
        ## open files error'. The number of open temp files seems to
        ## be smaller than with the default
        ## TODO configure memory and other java variables on pipeline start
        sambamba_sort1_c = 'sambamba sort --memory-limit 8GB -n -t 8 --tmpdir=' + bam_tempdir + ' -o ' + bam_temppath + ' ' + bam_path
        print("running sambamba")
        print(sambamba_sort1_c)

        os.system(sambamba_sort1_c)

        ## TODO remove input file if on scratch
        # print("deleting input file: " + bam_path)
        # os.remove(bam_path)

        ## Duplicates are marked using Picard MarkDuplicates (v2.4.1)
        ## with the parameter `ASSUME_SORT_ORDER=queryname

        bam_path = bam_temppath        
        bam_temppath = os.path.join(bam_tempdir, sample + "_picard.bam")

        metrics_temppath = os.path.join(bam_tempdir, sample + "_metrics.bam")
    
        picard_mark_duplicates_c = "MarkDuplicates I=" + bam_path + " O=" + bam_temppath + " M=" + metrics_temppath + " TMP_DIR=" + tmp_dir + " ASSUME_SORT_ORDER=queryname VERBOSITY=WARNING" 

        print("running MarkDuplicates")
        print(picard_mark_duplicates_c)

        os.system(picard_mark_duplicates_c)

        ## remove input file
        print("removing input file: " + bam_path)
        os.remove(bam_path)

        ## sort
        
        bam_path = bam_temppath
        bam_temppath = os.path.join(bam_tempdir, sample + "_sambamba_2.bam")

        ## increasing the memory limit seems to prevent the 'Too many
        ## open files error'. The number of open temp files seems to
        ## be smaller than with the default
        ## TODO configure memory and other java variables on pipeline start
        sambamba_sort2_c = 'sambamba sort --memory-limit 8GB -t 8 --tmpdir=' + bam_tempdir + ' -o '   + bam_temppath + ' ' + bam_path
        print("running sambamba")
        print(sambamba_sort2_c)

        os.system(sambamba_sort2_c)    
    
        print("copying " + bam_temppath + " to " + bam_outpath)
        shutil.copy(bam_temppath, bam_outpath)

        print('deleting temporary directory ' + bam_tempdir)
        shutil.rmtree(bam_tempdir)
        
        bam_dict = {}
        bam_dict[bam_key] = bam_outfile

        index_dict = {}
        
        ## todo include bam indices
        b_obj = Bam(sample = sample, bam_dict = bam_dict,
                        index_dict = index_dict,
                        bam_dir = bam_outdir)
        
        ## index
        b_obj =  b_obj.index()
        
        ## converting to cram
        c_obj =  b_obj.convert2Cram(fasta = fasta_path)

        b_obj = b_obj.delete()
        
        ## index
        c_obj =  c_obj.index()

        return c_obj


    def runCalibrate(self, outdir, reference_dir, fasta, sequence_index,
                     intervals_file, interval_padding, known_sites, known_sites_index,
                     scratchdir):

        sample = self.sample
        bam_dir = self.bam_dir
        bam_dict = self.bam_dict

        print('\n***   base recalibration   ***\n')

        ## only runs on first replicate if there are multiple replicates
        bam_key_list = bam_dict.keys()
        bam_key_list.sort()
        bam_key = bam_key_list[0]
        
        bam_outfile = sample + ".bam"
        bam_outpath = os.path.join(outdir, bam_outfile)

        cram_outfile = sample + ".cram"
        cram_outpath = os.path.join(outdir, cram_outfile)

        cram_idx = cram_outfile + ".crai"
        
        if os.path.isfile(cram_outpath):
            print("calibrated output file " + cram_outpath +
                  " exists, skipped processing bam")

            ## maybe just return self?
            bam_dict = {}
            bam_dict[bam_key] = cram_outfile

            index_dict = {}
            index_dict[bam_key] = cram_idx

            b_obj = Bam(sample = sample, bam_dict = bam_dict,
                            index_dict = index_dict, bam_dir = outdir)

            return b_obj

        else:
            print("preparing to calibrate bam file")


            ## avoiding race condition in creating directory when running
            ## multiple threads. python 3.2 has os.makedirs(folder,
            ## exist_ok=True)
            try:
                os.makedirs(outdir)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    ## Reraise if failed for reasons other than existing already
                    raise

            reference_files = [fasta, intervals_file] + sequence_index + known_sites + known_sites_index
            
            ref1 = ref.Ref(reference_dir = reference_dir,
                           reference_files = reference_files)

            
            if bool(scratchdir):

                print("Using temporary directory " + scratchdir)
                
                reference_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                               source_dir = reference_dir)
                
                ## better return a ref object that sits on scratch
                ref1.copy(target_dir = reference_scratchdir)

                bam_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                         source_dir = bam_dir)

                b_obj = self.copy(target_dir = bam_scratchdir)

                calibrated_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                                source_dir = outdir)

                b_obj = b_obj.calibrate(outdir = calibrated_scratchdir,
                                        reference_dir = reference_scratchdir,
                                        fasta = fasta, intervals_file = intervals_file,
                                        interval_padding = interval_padding,
                                        known_sites = known_sites)

                print('deleting input scratch directory ' + bam_scratchdir)
                shutil.rmtree(bam_scratchdir)
                
                b_obj = b_obj.copy(target_dir = outdir)

                print('deleting output scratch directory ' + calibrated_scratchdir)
                shutil.rmtree(calibrated_scratchdir)

            else:

                b_obj = self.calibrate(outdir = outdir,
                                       reference_dir = reference_dir,
                                       fasta = fasta, intervals_file = intervals_file,
                                       interval_padding = interval_padding,
                                       known_sites = known_sites)

        return b_obj



    def calibrate(self, outdir, reference_dir, fasta, intervals_file, interval_padding, known_sites):
 
        sample = self.sample
        bam_dir = self.bam_dir
        bam_dict = self.bam_dict
        index_dict = self.index_dict

        ## only runs on first replicate if there are multiple replicates
        bam_key_list = bam_dict.keys()
        bam_key_list.sort()
        bam_key = bam_key_list[0]
        
        bam_file = bam_dict[bam_key]
        bam_path = os.path.join(bam_dir, bam_file)

        if not os.path.exists(outdir):
            print("directory " + outdir + " does not exist.")
            sys.exit()

        recal_outfile = sample + "_recal.txt"
        recal_outpath = os.path.join(outdir, recal_outfile)

        ## fasta = reference_files.pop(0)
        fasta_path = os.path.join(reference_dir, fasta)

        ## interval files are optional
        intervals_path = os.path.join(reference_dir, intervals_file)
        
        if os.path.isfile(intervals_path):
            intervals_option = ' -L ' + intervals_path + ' -ip ' + str(interval_padding)
        else:
            intervals_option = ''
        
        bam_outfile = sample + ".bam"
        bam_outpath = os.path.join(outdir, bam_outfile)

        # cram_outfile = sample + ".cram"
        # cram_outpath = os.path.join(outdir, cram_outfile)

        ## know_sites_arg
        known_sites_string = ''
    
        for f1 in known_sites:
            f1_path = os.path.join(reference_dir, f1)
            f1_string = '--known-sites ' + f1_path
            known_sites_string = known_sites_string + f1_string + ' '


        gatk_recalibrate_c = "gatk BaseRecalibrator --verbosity=WARNING --preserve-qscores-less-than 6  -R " + fasta_path + " -I " + bam_path + " -O " + recal_outpath + intervals_option + ' ' + known_sites_string

        
        print("running GATK BaseRecalibrator")
        print(gatk_recalibrate_c)

        os.system(gatk_recalibrate_c)

        gatk_applybqsr_c = "gatk ApplyBQSR -verbosity=WARNING --read-filter MappingQualityNotZeroReadFilter --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -R " + fasta_path + " -I " + bam_path + " --bqsr-recal-file " + recal_outpath + ' -O ' + bam_outpath + intervals_option

        
        print("running GATK ApplyBQSR")
        print(gatk_applybqsr_c)

        os.system(gatk_applybqsr_c)

        ## remove input bam, currently this is done by the calling function
        # print('deleting scratch directory ' + bam_dir)
        # shutil.rmtree(bam_dir)

        ## output bam
        
        bam_dict = {}
        bam_dict[bam_key] = bam_outfile

        # apparently gatk calibration generates index files, I
        # integrate them into the bam object to be able to delete them

        index_dict = dnaseq.dictRecurse(dictionary = bam_dict,
                                        f = lambda x: re.sub('.bam$', '.bai', x))

        b_obj = Bam(sample = sample, bam_dict = bam_dict,
                    index_dict = index_dict,
                    bam_dir = outdir)
        
        ## converting to cram

        c_obj =  b_obj.convert2Cram(fasta = fasta_path)

        b_obj = b_obj.delete()
        
        ## index
        c_obj =  c_obj.index()
        
        return c_obj
    


    def runHaplotypeCaller(self, outdir, reference_dir, fasta, intervals_file, interval_padding, sequence_index, scratchdir):
    
        sample = self.sample
        bam_dict = self.bam_dict
        bam_dir = self.bam_dir

        print('\n***   variant calling   ***\n')

        vcf_outfile = sample + ".vcf.gz"
        vcf_outpath = os.path.join(outdir, vcf_outfile)

        vcf_idx = vcf_outfile + ".tbi"

        bam_key_list = bam_dict.keys()
        bam_key_list.sort()
        bam_key = bam_key_list[0]
        
        if os.path.isfile(vcf_outpath):
            print("vcf output file " + vcf_outpath +
                  " exists, skipped HaplotypeCaller")

            vcf_dict = {}
            vcf_dict[bam_key] = vcf_outfile

            tbi_dict = {}
            tbi_dict[bam_key] = vcf_idx

            
            vcf_obj = vcf.Vcf(project = sample, vcf_dict = vcf_dict, 
                              tbi_dict = tbi_dict, vcf_dir = outdir)

            return vcf_obj

        else:
            print("preparing to run HaplotypeCaller")
            
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

                bam_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                       source_dir = bam_dir)

                b_obj = self.copy(target_dir = bam_scratchdir)

                vcf_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                         source_dir = outdir)

                vcf_obj = b_obj.haplotypeCaller(outdir = vcf_scratchdir,
                                                reference_dir = reference_scratchdir,
                                                fasta = fasta, intervals_file = intervals_file,
                                                interval_padding = interval_padding)

                print('deleting input scratch directory ' + bam_scratchdir)
                shutil.rmtree(bam_scratchdir)
                
                ## copy output from scratch to output directory
                vcf_obj = vcf_obj.copy(target_dir = outdir)

                print('deleting output scratch directory ' + vcf_scratchdir)
                shutil.rmtree(vcf_scratchdir)

            else:

                vcf_obj = self.haplotypeCaller(outdir = outdir,
                                               reference_dir = reference_dir,
                                               fasta = fasta, intervals_file = intervals_file,
                                               interval_padding = interval_padding)

        return vcf_obj



    def haplotypeCaller(self, outdir, reference_dir, fasta, intervals_file, interval_padding):
 
        sample = self.sample
        bam_dir = self.bam_dir
        bam_dict = self.bam_dict

        ## only runs on first replicate if there are multiple replicates
        bam_key_list = bam_dict.keys()
        bam_key_list.sort()
        bam_key = bam_key_list[0]
        
        bam_file = bam_dict[bam_key]
        bam_path = os.path.join(bam_dir, bam_file)

        fasta_path = os.path.join(reference_dir, fasta)

        ## intervals option is optional
        intervals_path = os.path.join(reference_dir, intervals_file)
        if os.path.isfile(intervals_path):
            intervals_option = ' -L ' + intervals_path + ' -ip ' + str(interval_padding)
        else:
            intervals_option = ''

        vcf_outfile = sample + ".vcf.gz"
        vcf_outpath = os.path.join(outdir, vcf_outfile)

        vcf_idx = vcf_outfile + ".tbi"
        
        if not os.path.exists(outdir):
            print("directory " + outdir + " does not exist.")
            sys.exit()
            
        ## --java-options "-Xmx4g"
            
        hc_c = 'gatk HaplotypeCaller --verbosity=ERROR -R ' + fasta_path + ' -I ' + bam_path + ' -O ' + vcf_outpath + ' -ERC GVCF' + intervals_option

        print("running HaplotypeCaller")
        print(hc_c)

        os.system(hc_c)
                
        ## TODO: currently this runs through the whole dictionary, not
        ## only the first key as above
        vcf_dict = dnaseq.dictRecurse(dictionary = bam_dict,
                                    f = lambda x: re.sub('.cram', '.vcf.gz', x))

        ## index file names are vcf file name with '.tbi' suffix
        tbi_dict = dnaseq.dictRecurse(dictionary = vcf_dict,
                                        f = lambda x: x + '.tbi')
        
        v_obj = vcf.Vcf(project = sample, vcf_dict = vcf_dict,
                        tbi_dict = tbi_dict,
                        vcf_dir = outdir)
        
        return v_obj




def indexBamFile(bam_file, bam_dir):

    """ creates an index for a bam file

    It deals with bam or cram files as input

    Args:
        bam_file (str): The name of the bam or cram file to be indexed
        bam_dir (str): The path to the directory where the bam or cram
            file is located.

    Returns:
        The file name of the index file file

    """
    
    file_suffix = bam_file.split('.')[-1]
    
    if (file_suffix != 'bam') and (file_suffix != 'cram'):
        print("file for indexing: " + bam_file + " does not have the correct suffix (\'bam\' or \'cram\')")
        sys.exit()
        
    else:

        ## TODO check if index file exists
        
        index_suffix = re.sub('m$', 'i', file_suffix)
        
        bam_outpath = os.path.join(bam_dir, bam_file)
        samtools_index_c = "samtools index " + bam_outpath
            
        print("running samtools index")
        print(samtools_index_c)

        os.system(samtools_index_c)
    
        index_file = bam_file + '.' + index_suffix

        return index_file



def bam2Cram(bam_file, bam_dir, fasta):

    """ converts a bam file to cram format

    Args:
        bam_file (str): The name of the bam file to be converted
        bam_dir (str): The path to the directory where the bam file is
            located.
        fasta (str): The path to the genome reference file (fasta format).
    Returns:
        The file name of the index file file

    """
    
    bam_path = os.path.join(bam_dir, bam_file)

    cram_file = re.sub('.bam$', '.cram', bam_file)

    cram_outpath = os.path.join(bam_dir, cram_file)

    samtools_cram_c = "samtools view -C -T " + fasta + " -o " + cram_outpath + " " + bam_path
            
    print("converting bam to cram")
    print(samtools_cram_c)
    os.system(samtools_cram_c)

    return cram_file
