import sys
import argparse
import os
import os.path
import shutil
import filecmp
import re

def __version__():
    version = 'v0.1.4'
    return(version)

class Vividict(dict):
    """ enables nested dictionaries with path creation on the fly
    """
    
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value


## use argparse to parse input options
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg
        ## return open(arg, 'r')  # return an open file handle


def dictRecurse(dictionary, f, **args):

    """ runs a function recursively on elements nested dictionaries

    Args:
        dictionary: The name of the (nested) dictionary
        f: The function to be applied to the dictionary elements.
        args: additional arguments for the function.

    Example:
    
    """
    
    new = {}
    for k, v in dictionary.items():
        if isinstance(v, dict):
            new[k] = dictRecurse(dictionary = v, f = f, **args)
        elif not (v == ''):
            new[k] = f(v, **args)
    return new

def copyFile(file, source_dir, target_dir):
    
    """ copies files from a source directory to a target directory

    If the original directory 'dir' is /path/to/out/dir, and the
    scratch directory is /path/to/scratch/ a directory
    /path/to/scratch/dir is created. If there is already a file with
    the same name in the target directory, filecmp.cmp is used to test
    if the file is identical to the file in the source directory. If
    not the file is copied from source to target.

    Args:
        file (str): The name of the file to be copied
        source_dir (str): The path to the original directory (full
            path).
        target_dir (str): The target directory where the file will be
            copied to (full path).

    Returns:
        The file name of the copied file

    """

    source_path = os.path.join(source_dir, file)
    target_path = os.path.join(target_dir, file)
    
    if not os.path.isfile(source_path):
        print('file ' + source_path + ' does not exist,or is not a file skipped.')

        ## TODO get the right return value here
        ## return is equivalent to return None
        ## sys.exit()
            
    else:
        
    ## the file is copied if it does not exist or is different from source file
        if not os.path.exists(target_path) or not filecmp.cmp(source_path,
                                                            target_path):
            print('copying: ' + source_path + ' to ' + target_path)
            shutil.copy(source_path, target_path)
        else:
            print('file ' + target_path + ' already exists, skipped copying.')

    ## returns the file name of the copied file, maybe better return
    ## the path, but to do this I have to modify calling functions
    return file


def copytree(src, dst, symlinks=False, ignore=None):

    """copies a directory tree

    copytree implementation by Mital Vora, because with
    shutils.copytree, the target directory must not exist

    Args:
        src: The directory to be copied

        dst: the target directory, which is created if it does not
            exist.

        symlinks: 
        ignore:

    Returns:

    """

    if not os.path.exists(dst):
        os.makedirs(dst)
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            copytree(s, d, symlinks, ignore)
        else:
            if not os.path.exists(d) or os.stat(s).st_mtime - os.stat(d).st_mtime > 1:
                shutil.copy2(s, d)


def deleteFile(file, dir):

    """ deletes files specified by file name and directory

    Helper function to delete files

    Args:
        file (str): The name of the file to be copied
        dir (str): The path to the directory (full path).
    Returns:
        The file name of the deleted file

    """
    
    file_path = os.path.join(dir, file)
        
    if not os.path.isfile(file_path):
        print('file ' + file_path + ' does not exist or is not a file, skipped.')

        ## TODO get the right return value here
        ## return is equivalent to return None
        ## sys.exit()
            
    else:
        
        print('deleting: ' + file_path)
        os.remove(file_path)

        ## returns the file name of the copied file, maybe better return
        ## the path, but to do this I have to modify calling functions
        return file


def createScratchDir(source_dir, scratchdir):

    """ creates a directory on scratch, mirroring the output directory
    structure.

    If the original directory 'dir' is /path/to/out/dir, and the
    scratch directory is /path/to/scratch/ a directory
    /path/to/scratch/dir is created.

    Args:
        source_dir (str): The original directory (full path).
        scratchdir (str): The target (scratch) directory in which the
            new directory will be created. 
    
    Returns:
        The path to the newly created directory

    """

    ## needs to get rid of trailing '/' in order for basename working
    ## correctly
    source_dir = re.sub("/+$", "", source_dir)
    
    dir = os.path.join(scratchdir,
                       os.path.basename(source_dir))

    if not os.path.exists(dir):
        print("creating directory " + dir)
        os.makedirs(dir)

    return(dir)
