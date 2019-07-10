import sys
import os.path
import yaml
import argparse
import pickle
import subprocess
import shutil
import socket

bindir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(bindir, "../lib/"))

import dnaseq
import fastq

program = os.path.basename(sys.argv[0])
version = dnaseq.__version__()

date = subprocess.check_output('date')

## == parse import arguments ==
parser = argparse.ArgumentParser(description='running dnaseq pipeline',
                                 prog=program,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-r', '--run-mode', choices=['cluster', 'server', 'test'],
                    dest = 'run_mode', default='test', metavar = '',
                    help='mode in which to run, can be cluster, server, or test.')

parser.add_argument('-v', '--version',
                    ## metavar = '',
                    action = 'version', version='%(prog)s ' + version,
                    help='prints out the version of the program')

## this is a trick to list the following arguments as required arguments instead of optional arguments, see https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument('-c', '--config-file', dest = 'config_file', required = True,
                    help = 'path to yaml configuration file',
                    type=lambda x: dnaseq.is_valid_file(parser, x))

requiredNamed.add_argument('-f', '--fastq-object-file', dest = 'fastq_obj_file',
                    required = True,
                    help = 'path to serialized fastq object', metavar = 'FILE',
                    type = lambda x: dnaseq.is_valid_file(parser, x))

args = parser.parse_args()
   

yml_file = args.config_file
fastq_obj_file = args.fastq_obj_file
run_mode = args.run_mode

print('running ' + program + ' version ' + version)
print('\nstarting at: ' + date)

## == configurations from yaml file ==
yml_fh = open(yml_file, 'r')  # return an open file handle
cfg = yaml.load(yml_fh)

outdir = cfg['outdir']
mapper = cfg['mapper']

reference_dir = cfg['reference-dir']
scratchdir = cfg['scratchdir']
fasta = cfg['fasta']
sequence_index = cfg['sequence-index']

intervals_file = cfg['intervals-file']
interval_padding = cfg['interval-padding']

known_sites = cfg['known-sites']
known_sites_index = cfg['known-sites-index']


## == running in test (standalone) mode

## if running in test mode (without the initialize_alignment.py
## wrapper) and modules are used, load them here

modules = cfg['modules']

if run_mode == 'test' and modules: ## or run_mode == 'server':

   print("using environment modules")
   
   module_lib = cfg['module-lib']
   module_init = cfg['module-init']
   module_list = cfg['module-list']

   ## to get module environment working
   exec(open(module_init).read())

   ## setting the modules library
   module('use', '-a', module_lib)
   
   print('removing loaded modules....')
   module('purge')

   ## necessary to remove all white space here
   module_list = module_list.replace(" ", "")
   module_list = module_list.split(',')

   print("now loading modules....")

   for mod in module_list:
      module('load', mod)
    
   ## write to log file
   module('list')

      
## == fastq object ==

## load fastq object
fastq_obj = pickle.load(open(fastq_obj_file))

## remove fastq object file, unless we are in test mode
if not run_mode == 'test':
   os.remove(fastq_obj_file)

sample = fastq_obj.sample

## == output directory and lock file ==
if not os.path.exists(outdir):
    os.makedirs(outdir)

lfile = sample + '_alignment.lock'

lfile_path = os.path.join(outdir, lfile)

## problem: if there is running process, this overwrites its log file
if os.path.isfile(lfile_path):
    print('There is a file ' + lfile_path + ' that locks the sequence analysis job for ' +
          sample + '. Is there is another instance of this pipeline running for this sample?')
    sys.exit(0)
else:
   print("creating lock file " + lfile_path)
   fh = open(lfile_path, 'w+')
   fh.close
    
## == node environment
user = os.environ['USER']

if bool(scratchdir):
   scratchdir = os.path.join(scratchdir, user, sample)

   if not os.path.exists(scratchdir):
      os.makedirs(scratchdir)

if run_mode == 'cluster':
   print('User ' + user + ' running job ' + os.environ['JOB_NAME'] + ' with ID ' + os.environ['JOB_ID'] + ' on ' + os.environ['HOSTNAME'])
else:
   print('User ' + user + ' running sample ' + sample + ' on ' + os.uname()[1])

   
## == fastqc ==

fastqc_outdir = os.path.join(outdir, 'fastqc')

print('\ninitializing fastqc quality control...')

fastq_obj.runFastqc(fastqc_outdir = fastqc_outdir, scratchdir = scratchdir)


##  == alignment ==

bam_outdir = os.path.join(outdir, 'mapped-raw')

if mapper == "bwa":

    print('\ninitializing read alignment...')

    ## TODO try to cacth errors
    mapper_options=cfg['mapper-options']

    bam_obj = fastq_obj.runAlign(mapper_outdir = bam_outdir,
                                 reference_dir = reference_dir,
                                 fasta = fasta, mapper_options = mapper_options,
                                 sequence_index = sequence_index,
                                 scratchdir = scratchdir)

## == processing ==

print('\ninitializing bam processing...')

processed_outdir = os.path.join(outdir, 'mapped-processed')

bam_obj = bam_obj.runProcess(bam_outdir = processed_outdir,
                             scratchdir = scratchdir,
                             reference_dir = reference_dir,
                             fasta = fasta)

## == base calibration ==

print('\ninitializing base calibration...')

calibrated_outdir = os.path.join(outdir, 'mapped-calibrated')

bam_obj = bam_obj.runCalibrate(outdir = calibrated_outdir,
                               reference_dir = reference_dir,
                               fasta = fasta, sequence_index = sequence_index,
                               intervals_file = intervals_file,
                               interval_padding = interval_padding,
                               known_sites = known_sites,
                               known_sites_index = known_sites_index,
                               scratchdir = scratchdir)


## == variant calling ==

print('\ninitializing variant calling...')

vcf_outdir = os.path.join(outdir, 'variants-raw')

vcf_obj = bam_obj.runHaplotypeCaller(outdir = vcf_outdir,
                                     reference_dir = reference_dir,
                                     fasta = fasta, intervals_file = intervals_file,
                                     interval_padding = interval_padding,
                                     sequence_index = sequence_index,
                                     scratchdir = scratchdir)


## == clean up ==

print('\ncleaning up...')

## remove scratch directory, unless we are in test mode
if (not run_mode == 'test') and bool(scratchdir):
   print('deleting schratch directory ' + scratchdir)
   shutil.rmtree(scratchdir)

if os.path.isfile(lfile_path):
   print('deleting lock file ' + lfile_path)
   os.remove(lfile_path)
