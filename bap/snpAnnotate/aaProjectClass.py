import click
import os
import os.path
import sys
import shutil
import random
import string
import itertools
import time
import platform
from ruamel import yaml
from ..bapHelp import *

class aaProject():
	def __init__(self, script_dir, input, output, name, ncores,
		bwa_path, bwa_index):
		
				
		#----------------------------------
		# Assign straightforward attributes
		#----------------------------------
		self.script_dir = script_dir
		self.output = output
		self.name = name			
		self.bamfile = input			

		# Handle bwa coordination
		self.bwa = get_software_path('bwa', bwa_path)
		if(os.path.isfile(bwa_index)):
			self.bwa_index = bwa_index
		else: 
			sys.exit("Could not find the BWA index: %s" % bwa_index)
		
		if(name == "default"):
			filename, file_extension = os.path.splitext(self.bamfile)
			self.name = os.path.basename(filename)

		
	#--------------------------------------------------------------------------------
	# Define a method to dump the object as a .yaml/dictionary for use in other files
	#--------------------------------------------------------------------------------
	def __iter__(self):
		
		yield 'script_dir', self.script_dir
		yield 'mode', self.mode
		yield 'output', self.output
		yield 'bamfile', self.bamfile
		yield 'name', self.name

		
		yield 'bwa', self.bwa
		yield 'bwa_index', self.bwa_index
		
	