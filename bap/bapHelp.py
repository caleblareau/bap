import time
import shutil
import re
import os
import sys
import csv
import gzip
import glob
import subprocess
from functools import partial
from importlib.metadata import version, PackageNotFoundError


def bap_version():
	'''
	Return the installed bap-atac version, or "unknown" if the package
	metadata cannot be located (e.g. running from an uninstalled checkout).
	'''
	try:
		return version("bap-atac")
	except PackageNotFoundError:
		return "unknown"


def run_cmd(args, log_file=None, log_stderr_only=False, shell=False, check=True):
	'''
	Run an external command from a list of arguments (no shell by default),
	raising a clear error if it fails. This replaces the historical
	os.system(<string>) calls, avoiding shell-injection and quoting issues.

	args:            list of command arguments (or a string when shell=True)
	log_file:        optional path; command output is written here. By default
	                 both stdout and stderr are captured (equivalent to `&>`).
	log_stderr_only: when True, only stderr is captured to log_file and stdout
	                 is inherited by the terminal (equivalent to `2>`).
	shell:           only set True for commands that need a shell pipeline
	check:           when True, exit with a clear message on a non-zero return
	                 code; set False to let the caller inspect the result
	'''
	if log_file is not None:
		with open(log_file, 'w') as fh:
			if log_stderr_only:
				completed = subprocess.run(args, stderr=fh, shell=shell)
			else:
				completed = subprocess.run(args, stdout=fh, stderr=subprocess.STDOUT, shell=shell)
	else:
		completed = subprocess.run(args, shell=shell)
	if check and completed.returncode != 0:
		printed = args if isinstance(args, str) else " ".join(str(a) for a in args)
		sys.exit(gettime() + "ERROR: command failed (exit %d): %s" % (completed.returncode, printed))
	return completed


def string_hamming_distance(str1, str2):
	'''
	Fast hamming distance over 2 strings known to be of same length.
	In information theory, the Hamming distance between two strings of equal
	length is the number of positions at which the corresponding symbols
	are different.
	eg "karolin" and "kathrin" is 3.
	'''
	return sum(c1 != c2 for c1, c2 in zip(str1, str2))

def intersection(lst1, lst2):
	lst3 = [value for value in lst1 if value in lst2]
	return lst3

def make_folder(folder):
	"""
	Function to only make a given folder if it does not already exist
	"""
	if not os.path.exists(folder):
		os.makedirs(folder)

def file_len(fname):
	"""
	Function to count number of lines in a file
	"""
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

def verify_file(filename):
	"""
	Ensure that file can both be read and exists
	"""
	try:
		extension = filename.split('.')[-1]
		if extension == "gz":
			fp = gzip.open(filename, 'rt')
		else:
			fp = open(filename)
		one = fp.readline()
		fp.close()
	except IOError as err:
		sys.exit(gettime() + "Error reading the file {0}: {1}".format(filename, err))
	return(filename)


def get_software_path(tool, abs_path):
	'''
	Function takes a tool name and a possible absolute path
	specified by the user input and returns the absolute path
	to the tool
	'''
	tool_path = shutil.which(tool)
	if(str(tool_path) == "None"):
		if(abs_path == ""):
			sys.exit("ERROR: cannot find "+tool+" in environment; add it to user PATH environment or specify executable using a flag.")
	if(abs_path != ""):
		if(os.path.isfile(abs_path)):
			tool_path = abs_path
	return(tool_path)

def rev_comp(seq):
	'''
	Fast Reverse Compliment
	'''  
	tbl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
	return ''.join(tbl[s] for s in seq[::-1])

def gettime(): 
	'''
	Matches `date` in unix
	'''
	return(time.strftime("%a ") + time.strftime("%b ") + time.strftime("%d ") + time.strftime("%X ") + 
		time.strftime("%Z ") + time.strftime("%Y")+ ": ")
		

def findIdx(list1, list2):
	'''
	Return the indices of list1 in list2
	'''
	return [i for i, x in enumerate(list1) if x in list2]

def check_R_packages(required_packages, R_path):
	'''
	Determines whether or not R packages are properly installed
	'''
	completed = subprocess.run(
		[R_path, "--vanilla", "--slave", "-e", "cat(rownames(installed.packages()), sep='\\n')"],
		stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
	installed_packages = completed.stdout.strip().split("\n")
	if(not set(required_packages) < set(installed_packages)):
		sys.exit("ERROR: cannot find the following R package: " + str(set(required_packages) - set(installed_packages)) + "\n" + 
			"Install it in your R console and then try rerunning proatac (but there may be other missing dependencies).")
		

def inferSampleVectors(input):
	'''
	Given input file / directory, creates 3 vectors of samplename, fastq1, fastq2
	'''
	def list_duplicates_of(seq,item):
		start_at = -1
		locs = []
		while True:
			try:
				loc = seq.index(item,start_at+1)
			except ValueError:
				break
			else:
				locs.append(loc)
				start_at = loc
		return locs

	samplenames = []
	fastq1 = []
	fastq2 = []

	# Import sample table if that's what is provided
	if(os.path.isfile(input)):
		sniff_range = 4096
		delimiters = ';\t, '
		sniffer = csv.Sniffer()

		with open(input, 'r') as infile:
			# Determine dialect
			dialect = sniffer.sniff(
				infile.read(sniff_range), delimiters=delimiters
			)
			infile.seek(0)

			# Sniff for header
			has_header = sniffer.has_header(infile.read(sniff_range))
			infile.seek(0)
			reader = csv.reader(infile, dialect)
			for row in reader:
				if(len(row) > 0):
					while '' in row:
						row.remove('')
					samplenames.append(row[0])
					fastq1.append(verify_file(row[1]))
					fastq2.append(verify_file(row[2]))
	
	# Otherwise figure out .fastq files from directory, do the merging, and infer sample names
	else:
		files = sorted(glob.glob(input.rstrip("/") + "/*.fastq.gz"))
		if(len(files) < 2):
			sys.exit("ERROR: input determined to be a directory but no paired .fastq.gz files found; QUITTING")
		files1 = [filename.replace("_R1", "_1").replace("_R2", "_1").replace("_2", "_1") for filename in files]
		samples = [os.path.basename(sample).split("_1")[0] for sample in files1]
		dups_in_source = partial(list_duplicates_of, samples)
		for c in set(samples):
			if(len(dups_in_source(c)) == 2):
				samplenames.append(c)
				fastq1.append(verify_file(files[dups_in_source(c)[0]]))
				fastq2.append(verify_file(files[dups_in_source(c)[1]]))

	return(samplenames, fastq1, fastq2)

def filterExistingSamples(samples, fastq1, fastq2, output):
	'''
	Looks in the output folder to see which samples have been 
	terminally processed already and will filter them from the three vectors
	of sample names, fastqs, etc. 
	'''
	# use findIdx here to infer which samples
	# already exist
	return(samples, fastq1, fastq2)

# Process through iterator
def batch_iterator(iterator, batch_size):
	entry = True  # Make sure we loop once
	while entry:
		batch = []
		while len(batch) < batch_size:
			try:
				entry = iterator.__next__()
			except StopIteration:
				entry = None
			if entry is None:
				# End of file
				break
			batch.append(entry)
		if batch:
			yield batch


# https://stackoverflow.com/questions/1006289/how-to-find-out-the-number-of-cpus-using-python	
def available_cpu_count():
	'''
	Number of available virtual or physical CPUs on this system, i.e.
	user/real as output by time(1) when called with an optimally scaling
	userspace-only program
	'''
	# cpuset
	# cpuset may restrict the number of *available* processors
	try:
		m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$',
					  open('/proc/self/status').read())
		if m:
			res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
			if res > 0:
				return res
	except IOError:
		pass

	# Python 2.6+
	try:
		import multiprocessing
		return multiprocessing.cpu_count()
	except (ImportError, NotImplementedError):
		pass

	# http://code.google.com/p/psutil/
	try:
		import psutil
		return psutil.cpu_count()   # psutil.NUM_CPUS on old versions
	except (ImportError, AttributeError):
		pass

	# POSIX
	try:
		res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

		if res > 0:
			return res
	except (AttributeError, ValueError):
		pass

	# Windows
	try:
		res = int(os.environ['NUMBER_OF_PROCESSORS'])

		if res > 0:
			return res
	except (KeyError, ValueError):
		pass

	# jython
	try:
		from java.lang import Runtime
		runtime = Runtime.getRuntime()
		res = runtime.availableProcessors()
		if res > 0:
			return res
	except ImportError:
		pass

	# BSD
	try:
		sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
								  stdout=subprocess.PIPE)
		scStdout = sysctl.communicate()[0]
		res = int(scStdout)

		if res > 0:
			return res
	except (OSError, ValueError):
		pass

	# Linux
	try:
		res = open('/proc/cpuinfo').read().count('processor\t:')

		if res > 0:
			return res
	except IOError:
		pass

	# Solaris
	try:
		pseudoDevices = os.listdir('/devices/pseudo/')
		res = 0
		for pd in pseudoDevices:
			if re.match(r'^cpuid@[0-9]+$', pd):
				res += 1

		if res > 0:
			return res
	except OSError:
		pass

	# Other UNIXes (heuristic)
	try:
		try:
			dmesg = open('/var/run/dmesg.boot').read()
		except IOError:
			dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
			dmesg = dmesgProcess.communicate()[0]

		res = 0
		while '\ncpu' + str(res) + ':' in dmesg:
			res += 1

		if res > 0:
			return res
	except OSError:
		pass

	raise Exception('Can not determine number of CPUs on this system')
