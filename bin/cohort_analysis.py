import sys
import os.path
import yaml
import argparse
import pickle
import subprocess
import shutil
import socket
from datetime import datetime

bindir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(bindir, "../lib/"))

import dnaseq
import vcf

program = os.path.basename(sys.argv[0])
version = dnaseq.__version__()
host = socket.gethostname()

## == parse import arguments ==
parser = argparse.ArgumentParser(description='running cohort variant pipeline',
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

requiredNamed.add_argument('-f', '--vcf-object-file', dest = 'vcf_obj_file',
                    required = True,
                    help = 'path to serialized vcf object', metavar = 'FILE',
                    type = lambda x: dnaseq.is_valid_file(parser, x))

args = parser.parse_args()
   

yml_file = args.config_file
vcf_obj_file = args.vcf_obj_file
run_mode = args.run_mode

print('\nProgram: ' + program)
print('\nVersion: ' + version)
print('\nHost: ' + host)
print('\nStart time: ' + str(datetime.now()))

## == configurations from yaml file ==
yml_fh = open(yml_file, 'r')  # return an open file handle
cfg = yaml.load(yml_fh, Loader=yaml.FullLoader)

outdir = cfg['outdir']

reference_dir = cfg['reference-dir']
scratchdir = cfg['scratchdir']
fasta = cfg['fasta']
sequence_index = cfg['sequence-index']

intervals_file = cfg['intervals-file']
interval_padding = cfg['interval-padding']

known_sites = cfg['known-sites']
known_sites_index = cfg['known-sites-index']
snp_resource_string = cfg['snp-resource-string']
indel_resource_string = cfg['indel-resource-string']

filter_string = cfg['filter-string']

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

      
## == vcf object ==

## load vcf object
vcf_obj = pickle.load(open(vcf_obj_file, 'rb'))

## remove vcf object file, unless we are in test mode
if not run_mode == 'test':
   os.remove(vcf_obj_file)

project = vcf_obj.project

## == output directory and lock file ==
if not os.path.exists(outdir):
    os.makedirs(outdir)

lfile = project + '_variant_calling.lock'

lfile_path = os.path.join(outdir, lfile)

## problem: if there is running process, this overwrites its log file
if os.path.isfile(lfile_path):
    print('There is a file ' + lfile_path + ' that locks the vcf analysis job for ' +
          project + '. Is there is another instance of this pipeline running for this sample?')
    sys.exit(0)
else:
   print("creating lock file " + lfile_path)
   fh = open(lfile_path, 'w+')
   fh.close
    
## == node environment
user = os.environ['USER']

if bool(scratchdir):
   scratchdir = os.path.join(scratchdir, user, project)

   if not os.path.exists(scratchdir):
      os.makedirs(scratchdir)

if run_mode == 'cluster':
   print('User ' + user + ' running job ' + os.environ['JOB_NAME'] + ' with ID ' + os.environ['JOB_ID'] + ' on ' + host)
else:
   print('User ' + user + ' running cohort ' + project + ' on ' + os.uname()[1])

   
## == combining and genotyping vcf files ==

geno_outdir = os.path.join(outdir, 'variants-genotyped')

from pprint import pprint
pprint(vcf_obj.vcf_dict)

print('\ninitializing combining and genotyping vcf files...')
print('\nstarting at: ' + str(datetime.now()))

vcf_obj = vcf_obj.runGenotype(outdir = geno_outdir,
                              reference_dir = reference_dir,
                              fasta = fasta, sequence_index= sequence_index,
                              intervals_file = intervals_file,
                              interval_padding = interval_padding,
                              filter_string = filter_string,
                              scratchdir = scratchdir)


##  == variant calibration ==

calibrated_outdir = os.path.join(outdir, 'variants-calibrated')

print('\ninitializing variant calibration...')
print('\nstarting at: ' + str(datetime.now()))

vcf_obj = vcf_obj.runCalibrate(outdir = calibrated_outdir,
                               reference_dir = reference_dir,
                               fasta = fasta, sequence_index = sequence_index,
                               known_sites = known_sites,
                               known_sites_index = known_sites_index,
                               snp_resource_string = snp_resource_string,
                               indel_resource_string = indel_resource_string,
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

print('\nFinish time: ' + str(datetime.now()))
