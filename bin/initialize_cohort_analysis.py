import os.path
import subprocess
import yaml
import re
import pickle
import argparse
import sys
import uuid

## path to library files
bindir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(bindir, "../lib/"))

import dnaseq
import vcf

program = os.path.basename(sys.argv[0])
version = dnaseq.__version__()


## parse import arguments
parser = argparse.ArgumentParser(description = "initializing dnaseq cohort analysis pipeline")

## this is a trick to list the following arguments as required arguments instead of optional arguments, see https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument("-c", "--config-file", dest = "config_file", required = True,
                    help = "path to yaml configuration file",
                    type = lambda x: dnaseq.is_valid_file(parser, x))

parser.add_argument("-r", "--run-mode", choices=['cluster', 'server', 'test'],
                    dest = "run_mode", default='test', metavar = '',
                    help="mode in which to run, can be cluster, server, test")

parser.add_argument('-v', '--version',
                    ## metavar = '',
                    action = 'version', version='%(prog)s ' + version,
                    help='prints out the version of the program')

args = parser.parse_args()


print('running ' + program + ' version ' + version)
print('\nstarting at: ' + subprocess.check_output('date'))

## get configurations from yaml file
yml_file = args.config_file
## returns an open file handle
yml_fh = open(yml_file, 'r')

cfg = yaml.load(yml_fh)

## run mode
run_mode = args.run_mode

cohorts_file = cfg['cohorts-file']
outdir = cfg['outdir']
reference_dir = cfg['reference-dir']
vcf_dir = cfg['vcf-dir']
grep_column = cfg['grep-column']
grep_term = cfg['grep-term']
pe = cfg['paired-end']


## == output and log directory ==

if not os.path.exists(outdir):
   os.makedirs(outdir)

log_dir = os.path.join(outdir, 'logs')
if not os.path.exists(log_dir):
   os.makedirs(log_dir)
   
## == environment modules ==

modules = cfg['modules']

if(modules):
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


## == cohort dictionary ==

cohort_dict = dnaseq.Vividict()

## read in data
with open(cohorts_file) as f:
    ## skip first line
   next(f)
    
   for line in f:

      ## skip lines that start with non-word characters,
      ## e.g. commented-out lines
      if re.match("\W+",line):
         continue

      in_array = line.split()
      
      ## use grep to select lines
      if not re.search(grep_term, in_array[grep_column]):
         continue
            
      cohort = in_array[0]

      if not cohort in cohort_dict:            
         cohort_dict[cohort]["replicate_count"] = 0
            
      cohort_dict[cohort]["replicate_count"] += 1
            
      replicate_key = 'sample_' + str(cohort_dict[cohort]["replicate_count"]) 

      vcf_file = in_array[1]
        
      cohort_dict[cohort]['cohort-dict'][replicate_key] = vcf_file
                                        
print(cohort_dict)



## == vcf object ==

for cohort in cohort_dict:
    
   vcf = vcf.Vcf(project = cohort,
                 vcf_dict = cohort_dict[cohort]['cohort-dict'],
                 tbi_dict = {},
                 vcf_dir = vcf_dir)

   vcf = vcf.index()
   
   ## to make the file name unique
   tmp = uuid.uuid4().hex
   binary_file_name = cohort + '_' + tmp + '.pkl'

   binary_file = open(binary_file_name, mode = 'wb')
   my_pickled_vcf = pickle.dump(vcf, binary_file)
   binary_file.close()
    
   print("\nInitialising job for cohort: " + cohort)

   pipeline_command = 'python ' + bindir + '/cohort_analysis.py --config-file ' + yml_file + ' --vcf-object-file ' + binary_file_name + ' --run-mode ' + run_mode

   ## qsub variables

   if run_mode == 'cluster' or run_mode == 'test':
      mem = cfg['mem']
      tscratch = cfg['tscratch']
      time = cfg['time']

      qsub_options = '-S /bin/bash -o ' + log_dir + ' -e ' + log_dir + ' -cwd -l tmem=' + mem + ',tscratch=' + tscratch + ',h_rt=' + time + ' -V -N ' + cohort + '_variant_calling'

   
   if run_mode == 'cluster':
      print "\nrunning the pipeline on the sge queue"

      print('\nvariant calling command: ' + pipeline_command)
      
      print('\nqsub options: ' + qsub_options)

      qsub_command = 'qsub ' + qsub_options + '<<EOF\n' +  pipeline_command + '\nEOF'

      print('\nqsub command: ' + qsub_command)

      os.system(qsub_command)


   if run_mode == 'test':
      print "\ndry run for tests"

      print('\nvariant calling command: ' + pipeline_command)
      
      print('\nqsub options: ' + qsub_options)

      qsub_command = 'qsub ' + qsub_options + '<<EOF\n' +  pipeline_command + '\nEOF'

      print('\nqsub command: ' + qsub_command)

      
   if run_mode == 'server':
      print "\nrunning the pipeline as a subprocess on the server"

      print('\nvariant calling command: ' + pipeline_command)

      stdout_log = os.path.join(log_dir, cohort + '_stdout.txt')
      fh_stdout = open(stdout_log, 'w')

      stderr_log = os.path.join(log_dir, cohort + '_stderr.txt')
      fh_stderr = open(stderr_log, 'w')

      
      # subprocess needs the command as a list
      p = subprocess.Popen(pipeline_command.split(' '), stdout=fh_stdout,
                           stderr=fh_stderr)
