"""
The fastq module
"""

import os
import os.path
import shutil
import filecmp
import glob
import re
import errno

import ref
import bam
import fastqc
import dnaseq


class Fastq:
    """Fastq is a class and data structure that contains information
    about the fastq files of a sample with the elements

    Attributes:
    -----------
    sample (str):
       sample name.
    fastq_dict:
       dictionary modeling the relationships of the fastq files (see details).
    fastq_dir(str)
       path to directory where fastq files are located
    -----------

    fastq_dict is a dictionary of dictionaries that models the
    relationships of fastq files. Fastq files from a sample can have
    multiple replicates, (e.g. multiple lanes). Each replicate has a
    pair of fastq files (forward and reverse).


    |
    |-- fastq_key_1
    |             |-fastq_1
    |             |-fastq_2
    |
    |-- fastq_key_2
                  |-fastq_1
                  |-fastq_2

    """
    
    # The class "constructor" - It's actually an initializer 
    def __init__(self, sample, fastq_dict, fastq_dir):
        self.sample = sample
        self.fastq_dict = fastq_dict
        self.fastq_dir = fastq_dir

    def fastqList(self):

        """
        write fastq file names in a fastq object to a list
        """
        ## list of fastq file names
        fastq_list=[]

        fastq_dict = self.fastq_dict
        rep_keys = fastq_dict.keys()

        for rep in fastq_dict.keys():
            for fq in fastq_dict[rep].keys():
                fastq_list.append(fastq_dict[rep][fq])

        return fastq_list

    def fastqRepList(self):

        """
        transposes the fasq_dict dictionary, returning a list of lists:
        1. level: fastq_keys, usually forward and reverse
        2. level: rep_keys, replicates for the same read direction of the same sample
        """
        
        fastq_dict = self.fastq_dict

        rep_keys = fastq_dict.keys()
        list.sort(rep_keys)
        
        fastq_keys = fastq_dict[rep_keys[0]].keys()
        list.sort(fastq_keys)
        
        fastq_rep_list = []
        
        for i in range(len(fastq_keys)):
            fastq_rep_list.append([])
                
        for i in range(len(rep_keys)):
            for j in range(len(fastq_keys)):
                fastq_rep_list[j].append(fastq_dict[rep_keys[i]][fastq_keys[j]])

        return(fastq_rep_list)
        
    
    def concat(self, concat_dir):
        
        ## concatenate fastq files and generate an array of the file names for bwa
        ## returns a list of files with concatenated input files
        sample = self.sample
        fastq_dir = self.fastq_dir
        
        if not os.path.exists(concat_dir):
            os.makedirs(concat_dir)

        fastq_rep_list = self.fastqRepList()

        ## use list comprehension for list of lists
        fastq_rep_list = [[ fastq_dir + "/" + rep for rep in reps] for reps in fastq_rep_list]
        
        concat_file_list = []
        
        for i in range(len(fastq_rep_list)):
            file_count = str(i + 1)
            concat_file = concat_dir + "/" + sample + ".tmp." + str(i + 1) + ".fastq.gz";

            concat_file_list.append(concat_file)
            
            concat_string = ' '.join(fastq_rep_list[i])

            print('concatenating fastq reads')
            
            concat_c = 'zcat ' + concat_string + ' | gzip > ' + concat_file

            print(concat_c)
            
            os.system(concat_c)

        return(concat_file_list)


    def runAlign(self, mapper_outdir, reference_dir, fasta,
                 sequence_index, mapper_options, scratchdir):

        ## move fastq and reference to scratch, run bwa, and copy files back
        sample = self.sample
        fastq_dir = self.fastq_dir
        fastq_dict = self.fastq_dict

        print('\n***   bwa   ***\n')
        
        mapper_outfile = sample + ".bam"
        mapper_outpath = os.path.join(mapper_outdir, mapper_outfile) 
        
        ## if the path exists, just construct the bam object
        if os.path.isfile(mapper_outpath):
            print("bwa output file " + mapper_outpath + " exists, skipped running bwa")

            rep_key_list = fastq_dict.keys()
            rep_key_list.sort()
            rep_key = rep_key_list[0]

            bam_dict = {}
            bam_dict[rep_key] = mapper_outfile

            bam_obj = bam.Bam(sample = sample, bam_dict = bam_dict,
                                  bam_dir = mapper_outdir)

            
        else:
            print("preparing to run bwa")
            
            if not os.path.exists(mapper_outdir):
                os.makedirs(mapper_outdir)

            reference_files = [fasta] + sequence_index

            ref1 = ref.Ref(reference_dir = reference_dir,
                     reference_files = reference_files)    

            if bool(scratchdir):

                print("Using temporary directory " + scratchdir)


                fastq_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                           source_dir = fastq_dir)

                fastq_obj = self.copy(target_dir = fastq_scratchdir)

                mapper_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                            source_dir = mapper_outdir)
                
                reference_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                               source_dir = reference_dir)
            
                ## better return a ref object that sits on scratch
                ref1.copy(target_dir = reference_scratchdir)

                ## return bam object
                bam_obj = fastq_obj.align(mapper_outdir = mapper_scratchdir,
                                          reference_dir = reference_scratchdir,
                                          fasta = fasta,
                                          mapper_options = mapper_options)

                ## the fastq files currently are removed by the align function 
                print('deleting input scratch directory ' + fastq_scratchdir)
                shutil.rmtree(fastq_scratchdir)
                
                bam_obj = bam_obj.copy(target_dir = mapper_outdir)

                print('deleting output scratch directory ' + mapper_scratchdir)
                shutil.rmtree(mapper_scratchdir)

            else:
                
                bam_obj = self.align(mapper_outdir = mapper_outdir,
                                     reference_dir = reference_dir,
                                     fasta = fasta,
                                     mapper_options = mapper_options)

        return bam_obj


        
    def align(self, mapper_outdir, reference_dir, fasta, mapper_options):

        sample = self.sample
        fastq_dir = self.fastq_dir
        fastq_dict = self.fastq_dict
        
        fasta_path = os.path.join(reference_dir, fasta)

        mapper_tempdir = os.path.join(mapper_outdir, 'temp')
                    
        if not os.path.exists(mapper_tempdir):
            os.makedirs(mapper_tempdir)
            
        bam_dict = {}

        rep_key_list = fastq_dict.keys()
        list.sort(rep_key_list)
        
        for rep in rep_key_list:

            fastq_list = []
            
            for fq in fastq_dict[rep].keys():
                fastq_list.append(fastq_dict[rep][fq])

            ## DEBUG
            for fastq in fastq_list:
                print("fastq_dir:" + fastq_dir)
                print("fastq:" + fastq)
                
            fastq_rep_list = [fastq_dir + "/" + fastq for fastq in fastq_list]

            fastq_string = ' '.join(fastq_rep_list)
            
            fastq_file = fastq_list[0]
            
            fastq_base = fastq_file.split('.')[0]

            bwa_outfile = fastq_base + ".bam"
            bwa_outpath = os.path.join(mapper_tempdir, bwa_outfile)

            ## read groups, see
            ## https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
            ## library is missing, e.g. "\tLB:" + sample
            readgroup_string = "@RG\tID:" + fastq_base + "\tSM:" + sample + "\tPL:Illumina"

            ## MC and MQ tags are added using samblaster (v0.1.24)
            ## with the parameters `-a --addMateTags`

            # using verbosity level 2
            
            bwa_c = "bwa mem " + mapper_options + " -v 2 -R '" + readgroup_string + "' " + fasta_path + " " + fastq_string + " | samblaster -a --addMateTags | samtools view -b - > " + bwa_outpath

            print("running bwa")
            print(bwa_c)

            os.system(bwa_c)
            bam_dict[rep] = bwa_outfile

            ## TODO remove input fastq when on scratch
            # for fq1 in fastq_rep_list:
            #     print("removing " + fq1)
            #     os.remove(fq1)


        b_obj = bam.Bam(sample = sample, bam_dict = bam_dict,
                        bam_dir = mapper_tempdir)
        
        if (len(b_obj.bamList()) > 1):
            
            ## merge bam here (maybe as an option),
            ## TODO check if there are multiple bam
            print('merging bam files...')
            b_obj = b_obj.merge(bam_outdir = mapper_outdir)
            
        elif (len(b_obj.bamList()) == 1):

            k = b_obj.bam_dict.keys()[0]
            bam_path = os.path.join(b_obj.bam_dir, b_obj.bam_dict[k])
            bam_outfile = sample + ".bam"
            bam_outpath = os.path.join(mapper_outdir, bam_outfile)

            print("moving " + bam_path + " to " + bam_outpath)
            shutil.move(bam_path, bam_outpath)

            b_obj.bam_dir = mapper_outdir

            b_obj.bam_dict[k] = bam_outfile

            
        print('deleting temporary directory ' + mapper_tempdir)

        ## leave tempdir for debugging
        shutil.rmtree(mapper_tempdir)
        
        return b_obj

    
    def copy(self, target_dir):

        sample = self.sample
        fastq_dir = self.fastq_dir
        fastq_dict = self.fastq_dict

        ## in general, this should be generated by the calling function
        # if not os.path.exists(target_dir):
        #     os.makedirs(target_dir)

        fastq_dict = dnaseq.dictRecurse(dictionary = fastq_dict,
                                        f = dnaseq.copyFile,
                                        source_dir = fastq_dir,
                                        target_dir = target_dir)
    
        fastq_obj = Fastq(sample = sample, fastq_dict = fastq_dict,
                          fastq_dir = target_dir)

        return fastq_obj

    

    def runFastqc(self, fastqc_outdir, scratchdir):

        ## TODO run this iteratively through fastq_dict
        sample = self.sample
        fastq_dir = self.fastq_dir
        # first check of files on fastq_outdir exist
        fastq_list = self.fastqList()


        print('\n***   fastqc   ***\n')


        ## avoiding race condition in creating directory when running
        ## multiple threads. python 3.2 has os.makedirs(folder,
        ## exist_ok=True)
        try:
            os.makedirs(fastqc_outdir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                ## Reraise if failed for reasons other than existing already
                raise
        
        # if not os.path.exists(fastqc_outdir):
        #     os.makedirs(fastqc_outdir)
        
        fastqc_file_list = []
        # maybe create a fastqc object
        
        for fq in fastq_list:

            ## get rid of double extension fastq.gz, os.path.splitext does not work here
            fastqc_file = fq.split('.')[0] + "_fastqc.html"
            # fastqc_outpath = os.path.join(fastqc_outdir, fastqc_file)

            fastqc_file_list.append(fastqc_file)
            
        fastqc_obj = fastqc.Fastqc(sample = sample,
                                   fastqc_file_list = fastqc_file_list,
                                   fastqc_dir = fastqc_outdir)
            

        if (fastqc_obj.checkFiles()):

            print("fastqc files of sample " + sample
                  + " exist in directory " + fastqc_outdir
                  + ", skipped running fastqc")

        else:
            print("preparing to run fastqc")

            
            if bool(scratchdir):

                print("Using temporary directory " + scratchdir)

                fastq_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                           source_dir = fastq_dir)

                fastq_obj = self.copy(target_dir = fastq_scratchdir)

                ## maybe check for success
                fastqc_scratchdir = dnaseq.createScratchDir(scratchdir = scratchdir,
                                                            source_dir = fastqc_outdir)
            
                fastqc_obj = fastq_obj.fastqc(fastqc_scratchdir)

                print('deleting input scratch directory ' + fastq_scratchdir)
                shutil.rmtree(fastq_scratchdir)

                fastqc_obj = fastqc_obj.copy(target_dir = fastqc_outdir)

                print('deleting output scratch directory ' + fastqc_scratchdir)

                shutil.rmtree(fastqc_scratchdir)

            else:

                fastqc_obj = self.fastqc(fastqc_outdir)

            return fastqc_obj

        
    def fastqc(self, fastqc_outdir):

        sample = self.sample
        fastq_dir = self.fastq_dir

        fastq_list = self.fastqList()
                    
        fastqc_file_list = []
            
        for f1 in fastq_list:

            ## get rid of double extension fastq.gz, os.path.splitext does not work here
            fastq_base = f1.split('.')[0]
            
            fastq_path = os.path.join(fastq_dir, f1)
            
            ## fastqc_path = fastqc_outdir + fastq_base + '_fastqc.zip'
            ## print(fastqc_path)
            
            fastqc_c = "fastqc -q -o " + fastqc_outdir + " -f fastq " + fastq_path

            print("running fastqc")
            print(fastqc_c)

            os.system(fastqc_c)

            fastqc_file = f1.split('.')[0] + "_fastqc.html"
            fastqc_file_list.append(fastqc_file)
            

        fastqc_obj = fastqc.Fastqc(sample = sample,
                                   fastqc_file_list = fastqc_file_list,
                                   fastqc_dir = fastqc_outdir)
        
        return fastqc_obj
        
            
    def checkFiles(self):
        
        sample = self.sample
        fastq_dir = self.fastq_dir
        fastq_dict = self.fastq_dict

        if not os.path.isdir(fastq_dir):
            print("directory " + fastq_dir + " does not exist.")
            sys.exit()
                  
        fastq_list = self.fastqList()
        fastq_path_list = [os.path.join(fastq_dir, fq) for fq in fastq_list]

        test_list = [os.path.isfile(fp) for fp in fastq_path_list]
        
        return all(test_list)

