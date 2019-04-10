import sys
import shutil
import os
import unittest
from unittest.mock import patch

bindir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(bindir, "../lib/"))

## check how to avoid having to specify the full path
import dnaseq.dnaseq


## unit tests, see
## https://docs.python.org/3/library/unittest.html#module-unittest

## TODO use unittest.mock to avoid creating files, see
## https://stackoverflow.com/questions/37159714/python-creating-a-mock-or-fake-directory-with-files-for-unittesting/37161270

## https://www.toptal.com/python/an-introduction-to-mocking-in-python

## https://stackoverflow.com/questions/32748212/how-to-test-a-function-that-creates-a-directory

class createScratchTestCase(unittest.TestCase):
    
    def setUp(self):
        self.source_dir = os.path.join('unittest-out', 'source-dir')
        self.scratchdir = os.path.join('unittest-out', 'scratchdir')

    
    ## TODO use mock patches
    # @patch('dnaseq.dnaseq.os.path')
    # @patch('dnaseq.dnaseq.createScratchDir')
    
    def test_create_scratch_dir(self):        

        ## check what the mock_exists does
        ## mock_exists.return_value = True

        dnaseq.dnaseq.createScratchDir(source_dir = self.source_dir,
                                       scratchdir = self.scratchdir)

        self.assertTrue(os.path.exists(self.scratchdir), "scratch directory does not exist")

    def tearDown(self):
        shutil.rmtree(self.scratchdir)

        

class CopyFileTestCase(unittest.TestCase):


    def setUp(self):

        self.file = 'file.txt'
        self.source_dir = os.path.join('unittest-out', 'source_dir')
        self.target_dir = os.path.join('unittest-out', 'target_dir')

        self.source_path = os.path.join(self.source_dir, self.file)
        self.target_path = os.path.join(self.target_dir, self.file)


        if not os.path.exists(self.source_dir):
            os.makedirs(self.source_dir)

        if not os.path.exists(self.target_dir):
            os.makedirs(self.target_dir)
                   
        f = open(self.source_path,"w+")
        f.write("Hello, World!")
        f.close()

        ## TODO use mock patches
        ## @patch('dnaseq.dnaseq.os.path')
        ## @patch('dnaseq.dnaseq.copyFile')
        # @patch('dnaseq.dnaseq.shutil')

    def test_copy_file(self):        


        dnaseq.dnaseq.copyFile(file = self.file, source_dir = self.source_dir,
                  target_dir = self.target_dir)

        self.assertTrue(os.path.exists(self.target_dir), "target directory does not exist")
        
        self.assertTrue(os.path.isfile(self.target_path), "target file does not exist")

    def test_copy_file_return(self):        

        # check return value, should be file name
        rv = dnaseq.dnaseq.copyFile(file = self.file, source_dir = self.source_dir,
                                    target_dir = self.target_dir)

        self.assertEqual(self.file, rv)


    def tearDown(self):
        shutil.rmtree(self.source_dir)
        shutil.rmtree(self.target_dir)

        
if __name__ == '__main__':
    unittest.main()
