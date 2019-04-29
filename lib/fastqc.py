import os.path
import shutil
import glob
import sys

class Fastqc:


    """Fastqc is a class and data structure that contains information
    about the fastqc files of a sample with the elements

    Attributes:
    -----------
    sample (str):
       sample name.
    fastqc_file_list:
       list of fastqc output files for the sample.
    fastqc_dir(str)
       path to directory where fastqc files are located

    """


    # The class "constructor" - It's actually an initializer 
    def __init__(self, sample, fastqc_file_list, fastqc_dir):
        self.sample = sample
        ## problem: how to construct this list
	self.fastqc_file_list = fastqc_file_list
        self.fastqc_dir = fastqc_dir

    ## TODO dictionary  instead of file list
    def checkFiles(self):
        
        sample = self.sample
        fastqc_dir = self.fastqc_dir
        fastqc_file_list = self.fastqc_file_list

        if not os.path.isdir(fastqc_dir):
            print("directory " + fastqc_dir + " does not exist.")
            sys.exit()
                  
        # fastq_list = self.fastqList()
        fastqc_path_list = [os.path.join(fastqc_dir, fq) for fq in fastqc_file_list]

        test_list = [os.path.isfile(fp) for fp in fastqc_path_list]
        
        return all(test_list)


    ## TODO recursive copying from dictionary
    def copy(self, target_dir):
        
        sample = self.sample
        fastqc_dir = self.fastqc_dir
        fastqc_file_list = self.fastqc_file_list
        
        for fqc_file in fastqc_file_list:

            ## get rid of double extension fastq.gz, os.path.splitext does not work hee
            fastqc_base = os.path.splitext(fqc_file)[0]

            fastqc_base = fastqc_dir + "/" + fastqc_base + '*'

            for fqc_path in glob.glob(fastqc_base):
                
                fqc_base = os.path.basename(fqc_path)
                
                fqc_outpath = os.path.join(target_dir, fqc_base)
                
                print('copying: ' + fqc_path + ' to ' + fqc_outpath)

                shutil.copy(fqc_path, fqc_outpath)

