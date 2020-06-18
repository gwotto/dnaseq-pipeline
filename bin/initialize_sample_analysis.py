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
import fastq

program = os.path.basename(sys.argv[0])
version = dnaseq.__version__()

## parse import arguments
parser = argparse.ArgumentParser(description = "initializing dnaseq sample analysis pipeline")

## this is a trick to list the following arguments as required arguments instead of optional arguments, see https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument("-c", "--config-file", dest = "config_file", required = True,
                    help = "path to yaml configuration file",
                    type = lambda x: dnaseq.is_valid_file(parser, x))

parser.add_argument("-r", "--run-mode", choices=['cluster', 'server', 'test'],
                    dest = "run_mode", default='test', metavar = '',
                    help="mode in which to run, can be cluster, server, test, default is test")

parser.add_argument('-v', '--version',
                    ## metavar = '',
                    action = 'version', version='%(prog)s ' + version,
                    help='prints out the version of the program')

args = parser.parse_args()

## run mode
run_mode = args.run_mode

print('running ' + program + ' version ' + version)
print('\nstarting at: ' + subprocess.check_output('date'))

## get configurations from yaml file
yml_file = args.config_file
## returns an open file handle
yml_fh = open(yml_file, 'r')

## TODO get warning:
## YAMLLoadWarning: calling yaml.load() without Loader=... is deprecated, as the default Loader is unsafe. Please read https://msg.pyyaml.org/load for full details.
cfg = yaml.load(yml_fh)

samples_file = cfg['samples-file']
outdir = cfg['outdir']
reference_dir = cfg['reference-dir']
fastq_dir = cfg['fastq-dir']
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

## == sample dictionary ==

## dictionary
id_dict = dnaseq.Vividict()

## read in data
with open(samples_file) as f:
    ## skip first line
   next(f)
    
   for line in f:

      ## skip lines that start with non-word characters,
      ## e.g. commented-out lines
      if re.match("\W+",line):
         continue

      in_array=line.split()
        
      ## use grep to select lines
      if not re.search(grep_term,in_array[grep_column]):
         continue
            
      sample=in_array[0]

      ## initialize sample
      if not sample in id_dict:            
         id_dict[sample]["replicate_count"] = 0
            
      id_dict[sample]["replicate_count"] += 1
            
      replicate_key = 'replicate_' + str(id_dict[sample]["replicate_count"]) 
      fastq_1_id = in_array[1]
        
      id_dict[sample]["fastq_dict"][replicate_key]["fastq_1"] = fastq_1_id
                
      if pe:
         fastq_2_id = in_array[2]
         id_dict[sample]["fastq_dict"][replicate_key]["fastq_2"] = fastq_2_id
        
print('id and sample dictionary:\n')
print(id_dict)
print('\n')


## == fastq object and alignment command ==

for sample in id_dict:
    
   ## TODO replicate count in Fastq object
   fastq_obj = fastq.Fastq(sample = sample,
                           fastq_dict = id_dict[sample]["fastq_dict"],
                           fastq_dir = fastq_dir)

   ## to make the file name unique
   tmp = uuid.uuid4().hex
   binary_file_name = sample + '_' + tmp + '.pkl'

   binary_file = open(binary_file_name, mode = 'wb')
   my_pickled_fastq = pickle.dump(fastq_obj, binary_file)
   binary_file.close()
    
   print("\nInitialising job for sample: " + sample)

   pipeline_command = 'python ' + bindir + '/sample_analysis.py --config-file ' + yml_file + ' --fastq-object-file ' + binary_file_name + ' --run-mode ' + run_mode

   ## qsub variables

   if run_mode == 'cluster' or run_mode == 'test':
      mem = cfg['mem']
      tscratch = cfg['tscratch']
      time = cfg['time']

      qsub_options = '-S /bin/bash -o ' + log_dir + ' -e ' + log_dir + ' -cwd -l tmem=' + mem + ',tscratch=' + tscratch + ',h_rt=' + time + ' -V -N ' + sample + '_alignment'

   
   if run_mode == 'cluster':
      print('\nrunning the pipeline on the sge queue')

      print('\nalignment command: ' + pipeline_command)
      
      print('\nqsub options: ' + qsub_options)

      qsub_command = 'qsub ' + qsub_options + '<<EOF\n' +  pipeline_command + '\nEOF'

      print('\nqsub command: ' + qsub_command)

      os.system(qsub_command)


   if run_mode == 'test':
      print('\ndry run for tests')

      print('\nalignment command: ' + pipeline_command)
      
      print('\nqsub options: ' + qsub_options)

      qsub_command = 'qsub ' + qsub_options + '<<EOF\n' +  pipeline_command + '\nEOF'

      print('\nqsub command: ' + qsub_command)

      
   if run_mode == 'server':
      
      print('\nrunning the pipeline as a subprocess on the server')

      print('\nalignment command: ' + pipeline_command)

      stdout_log = os.path.join(log_dir, sample + '_stdout.txt')
      fh_stdout = open(stdout_log, 'w')

      stderr_log = os.path.join(log_dir, sample + '_stderr.txt')
      fh_stderr = open(stderr_log, 'w')

      
      # subprocess needs the command as a list
      p = subprocess.Popen(pipeline_command.split(' '), stdout=fh_stdout,
                           stderr=fh_stderr)
      
