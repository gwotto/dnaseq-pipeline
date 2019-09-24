import os.path
import shutil
import glob

import dnaseq

class Ref:
    # The class "constructor" - It's actually an initializer 
    def __init__(self, reference_dir, reference_files):
        self.reference_dir = reference_dir
        self.reference_files = reference_files

    ## NB pythonic way is not to use getters and setters


    def copy (self, target_dir):

        reference_dir = self.reference_dir
        reference_files = self.reference_files
        target_dir = target_dir
        
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
            
        print("reference files to copy: " + ', '.join(reference_files))
        # loop through the input filenames
        ## TODO what to do with directories?
        for r1 in reference_files:
            
            reference_path = os.path.join(reference_dir, r1)
            
            ## needs a '*'
            ## problem with globbing: only existing files
            ## are returned. If the pattern can not resolve to an
            ## existing file, the loop skips over it
            for r2 in glob.glob(reference_path):
            
                r2_base = os.path.basename(r2)
                
                reference_scratchpath = os.path.join(target_dir, r2_base)
            
                ## TODO check if it is the same file                
                dnaseq.copyFile(file = r2_base, source_dir = reference_dir, target_dir = target_dir)
