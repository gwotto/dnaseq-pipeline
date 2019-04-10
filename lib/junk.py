## junk yard of functions currently not needed



   def _runBwaOld(self, mapper_outdir, reference_dir, fasta,
                   mapper_index, mapper_options):

        sample = self.sample
        fastq_dir = self.fastq_dir

        fasta_path = os.path.join(reference_dir, fasta)
        bwa_outfile = sample + ".bam"
        bwa_outpath = os.path.join(mapper_outdir, bwa_outfile)

        concat_dir = mapper_outdir + '/' + 'fastq-concat'

        concat_file_list = self.fastqConcat(concat_dir=concat_dir)

        concat_string = ' '.join(concat_file_list)

        # Popen(["bwa", "mem", fasta_path, fastq_path], stdout = out_fh)

        readgroup_string = "@RG\tID:" + sample + "\tSM:" + sample + "\tLB:" + sample + "\tPL:Illumina"
        
        bwa_c = "bwa mem " + mapper_options + " -R '" + readgroup_string + "' " + fasta_path + " " + concat_string + " | samtools view -b - > " + bwa_outpath

        print("running bwa")
        print(bwa_c)

        os.system(bwa_c)

        ## remove concatenated fastq files
        shutil.rmtree(concat_dir) 

        ## TODO do this in a function programatically, in case there
        ## are multible bam files
        bam_dict = {}
        bam_dict['replicate_1'] = bwa_outfile

        bam_obj = bam.bam_obj(sample = sample, bam_dict = bam_dict,
                              bam_dir = mapper_outdir)

        return bam_obj

    def mergeBamOld(self, bam_outdir, scratchdir):

        sample = self.sample
        
        bam_outfile = sample + ".bam"
        bam_outpath = os.path.join(bam_outdir, bam_outfile)

        if os.path.isfile(bam_outpath):
            print("processed output file " + bam_outpath +
                  " exists, skipped processing bam")


        else:
            print("preparing to merge bam files")
            
            if not os.path.exists(bam_outdir):
                os.makedirs(bam_outdir)

            bam_scratchdir = os.path.join(scratchdir, 'bwa-raw')

            bam_scratchobj = self.copyBam(target_dir = bam_scratchdir)

            merged_scratchdir = os.path.join(scratchdir, 'bam-merged')

            bam_scratchobj = bam_scratchobj._mergeBam(bam_outdir = merged_scratchdir)

            bam_obj = bam_scratchobj.copyBam(target_dir = bam_outdir)
            
            return bam_obj



    def _mergeBam(self, bam_outdir):


        ## 1. MC and MQ tags are added using samblaster (v0.1.24) with
        ## the parameters `-a --addMateTags`.

        ## 2. Read group BAM files are merged together with `samtools
        ## merge` (v1.3.1-2).

        ## 3. The resulting file is name-sorted with `sambamba sort -n`
        ## (v0.6.4).
        
        ## 4. Duplicates are marked using Picard MarkDuplicates (v2.4.1)
        ## with the parameter `ASSUME_SORT_ORDER=queryname`

        ## 5. then the results are coordinate sorted using `sambamba
        ## sort`.
 
        sample = self.sample
        bam_dir = self.bam_dir
        bam_dict = self.bam_dict

        bam_outfile = sample + ".bam"
        bam_outpath = os.path.join(bam_outdir, bam_outfile)
        
        if not os.path.exists(bam_outdir):
            os.makedirs(bam_outdir)
            
        bam_list = self.bamList()

        bam_path_list = [bam_dir + "/" + bam for bam in bam_list]

        bam_string = ' '.join(bam_path_list)

        samtools_merge_c = "samtools merge " + bam_outpath + " " + bam_string

        print("running samtools merge")
        print(samtools_merge_c)

        os.system(samtools_merge_c)

        bam_dict = {}
        bam_dict['replicate_1'] = bam_outfile

        b_obj = bam_obj(sample = sample, bam_dict = bam_dict,
                              bam_dir = bam_outdir)

        
        ## todo include bam indices
        # b_obj = bam_obj(sample = sample, bam_dict = bam_dict,
        #                 index_dict = index_dict,
        #                 bam_dir = bam_outdir)

        return b_obj


    
    def _processBamOld(self, bam_outdir):


        ## 1. MC and MQ tags are added using samblaster (v0.1.24) with
        ## the parameters `-a --addMateTags`.

        ## 2. Read group BAM files are merged together with `samtools
        ## merge` (v1.3.1-2).

        ## 3. The resulting file is name-sorted with `sambamba sort -n`
        ## (v0.6.4).
        
        ## 4. Duplicates are marked using Picard MarkDuplicates (v2.4.1)
        ## with the parameter `ASSUME_SORT_ORDER=queryname`

        ## 5. then the results are coordinate sorted using `sambamba
        ## sort`.
 
        sample = self.sample
        bam_dir = self.bam_dir
        bam_dict = self.bam_dict

        bam_tempdir = os.path.join(bam_outdir, 'temp')
        
        if not os.path.exists(bam_outdir):
            os.makedirs(bam_outdir)
            
        if not os.path.exists(bam_tempdir):
            os.makedirs(bam_tempdir)
        
        bam_dict = dnaseq.dictRecurse(dictionary = bam_dict, f = processBamFile,
                                        bam_dir = bam_dir, bam_outdir = bam_outdir,
                                        bam_tempdir = bam_tempdir)

        index_dict = dnaseq.dictRecurse(dictionary = bam_dict, f = indexBamFile,
                                        bam_dir = bam_outdir)
        

        ## todo include bam indices
        b_obj = bam_obj(sample = sample, bam_dict = bam_dict,
                        index_dict = index_dict,
                        bam_dir = bam_outdir)

        return b_obj

    
    
def processBamFile(bam_file, bam_dir, bam_outdir, bam_tempdir):

    # bam_outpath = os.path.join(bam_outdir, bam_outfile)

    bam_base = bam_file.split('.')[0]

    bam_path = os.path.join(bam_dir, bam_file)
    bam_temppath = os.path.join(bam_tempdir, bam_base + "_sambamba_1.bam")

    # TOD maybe give processed file a 'processed' suffix
    bam_outfile = bam_base + '.bam'
    bam_outpath = os.path.join(bam_outdir, bam_file)
    
    sambamba_sort1_c = 'sambamba sort -n -o ' + bam_temppath + ' ' + bam_path
    print("running sambamba")
    print(sambamba_sort1_c)

    os.system(sambamba_sort1_c)    

    ## Duplicates are marked using Picard MarkDuplicates (v2.4.1)
    ## with the parameter `ASSUME_SORT_ORDER=queryname

    ## TODO  location of marked_dup_metrics.txt file
    bam_inpath = bam_temppath

    bam_temppath = os.path.join(bam_tempdir, bam_base + "_picard.bam")

    metrics_temppath = os.path.join(bam_tempdir, bam_base + "_metrics.bam")
    
    picard_mark_duplicates_c = "MarkDuplicates I=" + bam_inpath + " O=" + bam_temppath + " M=" + metrics_temppath + " ASSUME_SORT_ORDER=queryname" 

    print("running MarkDuplicates")
    print(picard_mark_duplicates_c)

    os.system(picard_mark_duplicates_c)

    ## sort and index
        
    bam_inpath = bam_temppath
    bam_temppath = os.path.join(bam_tempdir, bam_base + "_sambamba_2.bam")

    sambamba_sort2_c = 'sambamba sort -o ' + bam_temppath + ' ' + bam_inpath
    print("running sambamba")
    print(sambamba_sort2_c)

    os.system(sambamba_sort2_c)    
    
    print("copying " + bam_temppath + " to " + bam_outpath)
            
    shutil.copy(bam_temppath, bam_outpath)

    print('deleting temporary directory ' + bam_tempdir)
    shutil.rmtree(bam_tempdir)

    return bam_outfile


def indexBamFile(bam_file, bam_dir):

    ## maybe this can be replaced by a lambda function
    bam_inpath = os.path.join(bam_dir, bam_file)
    
    bam_index_file = bam_file + ".bai"
    bam_index_outpath = os.path.join(bam_dir, bam_index_file)

    samtools_index_c = "samtools index " + bam_inpath

    print("running samtools index")
    print(samtools_index_c)

    os.system(samtools_index_c)

    return bam_index_file


def calibrateBamFile(bam_file, sample, bam_dir, reference_dir, fasta, known_sites, outdir):

    bam_inpath = os.path.join(bam_dir, bam_file)

    recal_outfile = sample + "_recal.txt"
    recal_outpath = os.path.join(outdir, recal_outfile)

    ## fasta = reference_files.pop(0)
    fasta_path = os.path.join(reference_dir, fasta)

    bam_outfile = sample + ".bam"
    bam_outpath = os.path.join(outdir, bam_outfile)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ## know_sites_arg
    known_sites_string = ''
    
    for f1 in known_sites:
        f1_path = os.path.join(reference_dir, f1)
        f1_string = '-knownSites ' + f1_path
        known_sites_string = known_sites_string + f1_string + ' '
    
    gatk_recalibrate_c = "BaseRecalibrator --preserve_qscores_less_than 6 -dfrac .1 -nct 4  -R " + fasta_path + " -I " + bam_inpath + " -o " + recal_outpath + ' ' + known_sites_string

    print("running GATK BaseRecalibrator")
    print(gatk_recalibrate_c)

    os.system(gatk_recalibrate_c)

    gatk_printreads_c = "PrintReads --read_filter MappingQualityZero -preserveQ 6  -SQQ 10 -SQQ 20 -SQQ 30 --disable_indel_quals -R " + fasta_path + " -I " + bam_inpath + " -BQSR " + recal_outpath + ' -o ' + bam_outpath 

    print("running GATK PrintReads")
    print(gatk_printreads_c)

    os.system(gatk_printreads_c)
    
    return bam_outfile
