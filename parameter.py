#####################################################################################################################
# Copyright (c) 2017, St. Jude Children's Research Hospital															#
# All rights reserved.																								#
#																													#
# Redistribution and use in source and binary forms, with or without												#
# modification, are permitted provided that the following conditions are met:										#
#																													#
# 1. Redistributions of source code must retain the above copyright notice, this 									#
#    list of conditions and the following disclaimer.																#
#																													#
# 2. Redistributions in binary form must reproduce the above copyright notice,										#
#    this list of conditions and the following disclaimer in the documentation										#
#    and/or other materials provided with the distribution.															#
#																													#
# 3. Neither the name of the copyright holder nor the names of its contributors 									#
#    may be used to endorse or promote products derived from this software 											#
#    without specific prior written permission.																		#
#																													#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 									#
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED										#
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 											#
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 										#
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 										#
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 										#
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 										#
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 									#
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 									#
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 												#
#####################################################################################################################
# Helper Class for MICA 																							#
# Author: Alireza Khatamian (akhatami@stjude.org) 																	#
# Version: 1.0.0																									#
#####################################################################################################################
import os, csv, re, math, shlex, subprocess, sys, numpy

DEFAULT_HOST_OPTION_VALUE = 'local'
DEFAULT_DELIMITER = '\t'
DEFAULT_MEM_ALLOC = 8000
DEFAULT_DMIN = 16
DEFAULT_DMAX = 19
DEFAULT_KN = list(numpy.arange(3,15)) #[3, 4, 5, 6, 7, 8, 9, 10]
DEFAULT_B = 10
DEFAULT_MODE = 'clustering' 															# [clustering, validate, b-quality, c-quality, opt-k, overlay]
DEFAULT_INTERMEDIATE_PLOTTING = True
DEFAULT_INTERMEDIATE_PLOTTING_ = False
DEFAULT_MARKER_SIZE_THRESHOLD = 1000
MINIMUM_NORMAL_MARKER_SIZE = 20 * 3
MAXIMUM_NORMAL_MARKER_SIZE = 50 * 3
MINIMUM_SPECIAL_MARKER_SIZE = 5 * 3
MAXIMUM_SPECIAL_MARKER_SIZE = 10 * 3
MINIMUM_OPT_MARKER_SIZE = 10 * 3
MAXIMUM_OPT_MARKER_SIZE = 20 * 3
DEFAULT_RUN = 50
DEFAULT_VALUE_START_COLUMN_INDEX = 2
DEFAULT_DELIMITER = '\t'
TITLE_FONTSIZE = 14 + 8
SPECIAL_TITLE_FONTSIZE = 8 + 8
AXIS_LABEL_FONTSIZE = 12 + 8
TICK_LABEL_FONTSIZE = 8 + 8
TEXT_FONTSIZE = 10 + 8
HEADER_FONTSIZE = 18 + 8
MAX_PERPLEXITY = 1500
DEFAULT_DECOMPOSITION_METHOD = 'MDS'

class Parameter(object):
	def __init__(self, args):
		self.infile = None
		self.job_name = None
		self.host = None
		self.help = False
		self.true_labels = None
		self.mica_labels = None
		self.mode = DEFAULT_MODE
		self.intermediate_plotting = DEFAULT_INTERMEDIATE_PLOTTING
		self.intermediate_plotting_ = DEFAULT_INTERMEDIATE_PLOTTING_
		self.MARKER_SIZE_THRESHOLD = DEFAULT_MARKER_SIZE_THRESHOLD
		self.run = DEFAULT_RUN
		self.expression = None
		self.biomarkers = None
		self.decomposition = DEFAULT_DECOMPOSITION_METHOD
		for i in range(0, len(args)):
			if args[i] == '-h' or args[i] == '--help':
				Parameter.print_help()
				sys.exit()
		if len(args) < 2:
			Parameter.print_help()
			sys.exit()
		for i in range(0, len(args), 2):
			if args[i] == '-m' or args[i] == '--mode':
				self.mode = args[i + 1]
			elif args[i] == '-i' or args[i] == '--infile':
				self.infile = args[i + 1]
			elif args[i] == '-J' or args[i] == '--Job-Name':
				self.job_name = args[i + 1]
			elif args[i] == '-H' or args[i] == '--Host':
				if args[i + 1].startswith('local') or args[i + 1].startswith('LSF'):
					self.host = args[i + 1]
				else:
					print '[WARN] --> [PARM] Unrecognized host option: Host will be set to \'local\'.'
					self.host = DEFAULT_HOST_OPTION_VALUE
			elif args[i] == '-T' or args[i] == '--true-labels':
				self.true_labels = args[i + 1]
			elif args[i] == '-L' or args[i] == '--MICA_labels':
				self.mica_labels = args[i + 1]
			elif args[i] == '-R' or args[i] == '--Run':
				self.run = int(args[i + 1])
			elif args[i] == '-e' or args[i] == '--expression':
				self.expression = args[i + 1]
			elif args[i] == '-bm' or args[i] == '--biomarkers':
				self.biomarkers = args[i + 1]
			elif args[i] == '-d' or args[i] == '--decomposition':
				if args[i + 1] == 'MDS' or args[i + 1] == 'LPL' or args[i + 1] == 'PCA' or args[i + 1] == 'LPCA' or args[i + 1] == 'LPCA2':
					self.decomposition = args[i + 1]
				else:
					print '[WARN] --> Unrecognized decomposition method, set to default (' + str(DEFAULT_DECOMPOSITION_METHOD) + ')'
					self.decomposition = DEFAULT_DECOMPOSITION_METHOD
			else:
				print '[WARN] --> [PARM] Unknown parameter: Ignored. ' + args[i]

	def validate(self):
		self.intermediate_plotting = True
		if self.infile == None:
			print '[EROR] --> [PARM] Input file is required.'
			return False
		else:
			tokens = self.infile.split('?')
			file_path = tokens[0]
			if not os.path.exists(file_path):
				print '[EROR] --> [PARM] Input file does not exist.'
				return False
			else:
				file_type = tokens[0].split('/')[-1].split('.')[-1]
				if file_type == 'csv' or file_type == 'tsv' or file_type == 'txt':
					if len(tokens) > 1:
						infile_options = re.split('{|}| ', tokens[1])
						delimiter = DEFAULT_DELIMITER
						dmin = DEFAULT_DMIN
						dmax = DEFAULT_DMAX
						kn = []
						bs = numpy.arange(DEFAULT_B).tolist()
						for i in range(0, len(infile_options)):
							if infile_options[i] != '':
								option_tokens = infile_options[i].split('=')
								option_name = option_tokens[0]
								option_value = option_tokens[1]
								if option_name == 'dmin':
									if option_value.isdigit():
										dmin = int(option_value)
									else:
										print '[WARN] --> [PARM] dmin option must be a number, set to default (' + str(DEFAULT_DMIN) + ')'
								elif option_name == 'dmax':
									if option_value.isdigit():
										dmax = int(option_value)
									else:
										print '[WARN] --> [PARM] dmax option must be a number, set to default (' + str(DEFAULT_DMAX) + ')'
								elif option_name == 'kn':
									value_tokens = option_value.strip('[]\"').split(',')
									for v in value_tokens:
										v_tokens = v.split(':')
										if len(v_tokens) == 1:
											if v_tokens[0].isdigit():
												if int(v_tokens[0]) not in kn:
													kn.append(int(v_tokens[0]))
											else:
												print '[EROR] --> [PARM] Invalid k, ignored.'
										elif len(v_tokens) == 2:
											if v_tokens[0].isdigit() and v_tokens[1].isdigit():
												kmin = int(v_tokens[0])
												kmax = int(v_tokens[1])
												if kmin > kmax:
													print '[EROR] --> [PARM] Invalid range for k, ignored'
												else:
													kn += [x for x in range(kmin, kmax + 1) if x not in kn]
											else:
												print '[EROR] --> [PARM] Invalid range for k, ignored'
										else:
											print '[EROR] --> [PARM] Invalid value for k, ignored'
								elif option_name == 'B':
									if option_value.isdigit():
										bs = numpy.arange(int(option_value)).tolist()
									else:
										print '[WARN] --> [PARM] Number of bootstraps must be a number, set to default (' + str(DEFAULT_B) + ')'
								elif option_name == 'delimiter':
									if file_type == 'csv':
										print '[WARN] --> [PARM] Input file extension is .csv, delimiter is set to \',\''
										delimiter = ','
									elif file_type == 'tsv':
										print '[WARN] --> [PARM] Input file extension is .tsv, delimiter is set to \'\\t\''
										delimiter = '\t'
									elif file_type == 'txt':
										if option_value == ',':
											delimiter = ','
										elif option_value == '\\t':
											delimiter = '\t'
										else:
											print '[WARN] --> [PARM] Unrecognized delimiter.' 
											return False
								else:
									print '[WARN] --> [PARM] Unrecognized input file option: Ignored.'
						if not kn and self.mode != 'validate':
							kn = None #DEFAULT_KN
							print '[WARN] --> [PARM] kn will be determined by optimal K analysis on the following range: (\"' + str(DEFAULT_KN) + '\")'
						self.infile = dict(path=tokens[0], dmin=dmin, dmax=dmax, delimiter=delimiter, kn=kn, bs=bs)
					else:
						delimiter = DEFAULT_DELIMITER
						if file_type == 'csv':
							delimiter = ','
						elif file_type == 'tsv':
							delimiter = '\t'
						elif file_type == 'txt':
							print '[EROR] --> [PARM] Input file extension is .txt, please provide delimiter option.'
							return False
						self.infile = dict(path=tokens[0], dmin=DEFAULT_DMIN, dmax=DEFAULT_DMAX, delimiter=delimiter, kn=None, bs=numpy.arange(DEFAULT_B).tolist())
				else:
					print '[EROR] --> [PARM] Unrecognized file extension.'
					return False			
		if self.job_name == None:
			file_name = '.'.join(self.infile['path'].split('/')[-1].split('.')[:-1])
			print '[INFO] --> [PARM] Job name is set to \"' + file_name + '_MICA_' + self.decomposition + '\".'
			self.job_name = file_name + '_MICA_' + self.decomposition
		if self.host == None:
			print '[INFO] --> [PARM] Host is set to default (' + str(DEFAULT_HOST_OPTION_VALUE) + ').'
			self.host = dict(host=DEFAULT_HOST_OPTION_VALUE, mem_req=None)
		else:
			tokens = self.host.split('?')
			if tokens[0] == 'local':
				if len(tokens) > 1:
					print '[WARN] --> [PARM] \"local\" option has been chosen, additional host option(s) will be ignored.'
				self.host = dict(host=DEFAULT_HOST_OPTION_VALUE, mem_req=None)
			elif tokens[0] == 'LSF':
				if len(tokens) > 1:
					host_options = re.split('{|}| ', tokens[1])
					compute_mem = -1
					queue_name = None
					for i in range(0, len(host_options)):
						if host_options[i] != '':
							option_tokens = host_options[i].split('=')
							option_name = option_tokens[0]
							option_value = option_tokens[1]
							if option_name == 'cmem':
								if option_value.isdigit():
									compute_mem = int(option_value)
								else:
									print '[EROR] --> [PARM] Allocated memory must be an integer (MB), set to default (' + str(DEFAULT_MEM_ALLOC) + ').'
									compute_mem = DEFAULT_MEM_ALLOC
							elif option_name == 'queue':
								queue_name = option_value
							else:
								print '[EROR] --> [PARM] Unknown option: Ignored.'
					if compute_mem != -1:
						self.host = dict(host='LSF', mem_req=dict(cmem=compute_mem), queue=queue_name)
					else:
						print '[EROR] --> [PARM] Missing memory option(s) required, set to default (' + str(DEFAULT_MEM_ALLOC) + ').'
						if compute_mem == -1:
							compute_mem = DEFAULT_MEM_ALLOC
						self.host = dict(host='LSF', mem_req=dict(cmem=compute_mem), queue=queue_name)
				else:
					print '[WARN] --> [PARM] \"LSF\" host has been chosen, additional options set to default.'
					split_mem = -1
					if compute_mem == -1:
						compute_mem = DEFAULT_MEM_ALLOC
					self.host = dict(host='LSF', mem_req=dict(cmem=compute_mem), queue=None)
			else:
				print '[EROR] --> [PARM] Unknown host option, set to default (' + str(DEFAULT_HOST_OPTION_VALUE) + ').'
				self.host = dict(host=DEFAULT_HOST_OPTION_VALUE, mem_req=None, queue=None)
		if self.mode == 'validate' or self.mode == 'b-quality' or self.mode == 'c-quality':
			self.intermediate_plotting = False
			if self.true_labels != None:
				if not os.path.exists(self.true_labels):
					print '[EROR] --> [PARM] True label file does not exists, ignored.'
					self.true_labels = None
					return False
			else:
				print '[EROR] --> [PARM] True label file is required.'
				return False
		if self.mode == 'validate':
			if self.mica_labels != None:
				if not os.path.exists(self.mica_labels):
					print '[EROR] --> [PARM] MICA label file does not exists, ignored.'
					self.mica_labels = None
					return False
			else:
				print '[EROR] --> [PARM] MICA label file is required.'
				return False
			self.infile['kn'] = None
			self.infile['bs'] = None
			self.infile['dmin'] = DEFAULT_DMIN
			self.infile['dmax'] = DEFAULT_DMAX
		if self.mode == 'b-quality' or self.mode == 'c-quality':
			if self.infile['kn'] == None:
				print '[WARN] --> [PARM] kn will be determined by optimal K analysis on the following range: (\"' + str(DEFAULT_KN) + '\")'
			elif len(self.infile['kn']) > 1:
				print '[EROR] --> [PARM] Multiple Ks are not supported, set to the maximum provided K: [' + str(self.infile['kn'][len(self.infile['kn']) / 2]) + '].'
				self.infile['kn'] = [self.infile['kn'][len(self.infile['kn']) / 2]]
		if self.mode == 'b-quality':
			self.infile['bs'] = self.infile['bs'] #[self.infile['bs'][-1] + 1]
		if self.mode == 'opt-k':
			self.intermediate_plotting = False
			self.infile['bs'] = None
		if self.mode == 'opt-k-2':
			self.intermediate_plotting = False
		if self.mode == 'overlay':
			self.infile['kn'] = None
			self.infile['bs'] = None
			self.infile['dmin'] = DEFAULT_DMIN
			self.infile['dmax'] = DEFAULT_DMAX
			if self.expression != None:
				tokens = self.expression.split('?')
				file_path = tokens[0]
				if not os.path.exists(file_path):
					print '[EROR] --> [PARM] Expression file does not exists, ignored.'
					self.expression = None
					return False
				else:
					file_type = tokens[0].split('/')[-1].split('.')[-1]
					if file_type == 'csv' or file_type == 'tsv' or file_type == 'txt':
						if len(tokens) > 1:
							exp_options = re.split('{|}| ', tokens[1])
							ncol = DEFAULT_VALUE_START_COLUMN_INDEX
							delimiter = DEFAULT_DELIMITER
							for i in range(0, len(exp_options)):
								if exp_options[i] != '':
									option_tokens = exp_options[i].split('=')
									option_name = option_tokens[0]
									option_value = option_tokens[1]
									if option_name == 'ncol':
										if option_value.isdigit():
											ncol = int(option_value)
										else:
											print '[WARN] --> [PARM] ncol option of expression file must be an integer greater than 0, set to default (' + str(DEFAULT_SPLIT_OPTION_VALUE) + ')'
									elif option_name == 'delimiter':
										if file_type == 'csv':
											print '[WARN] --> [PARM] Expression file extension is .csv, delimiter is set to \',\''
											delimiter = ','
										elif file_type == 'tsv':
											print '[WARN] --> [PARM] Expression file extension is .tsv, delimiter is set to \'\\t\''
											delimiter = '\t'
										elif file_type == 'txt':
											if option_value == ',':
												delimiter = ','
											elif option_value == '\\t':
												delimiter = '\t'
											else:
												print '[WARN] --> [PARM] Unrecognized delimiter.' 
												return False
									else:
										print '[WARN] --> [PARM] Unrecognized expression file option: Ignored.'
							self.expression = dict(path=tokens[0], noheader=False, ncol=ncol, delimiter=delimiter)

						else:
							delimiter = DEFAULT_DELIMITER
							if file_type == 'csv':
								delimiter = ','
							elif file_type == 'tsv':
								delimiter = '\t'
							elif file_type == 'txt':
								print '[EROR] --> [PARM] Expression file extension is .txt, please provide delimiter option.'
								return False
							self.expression = dict(path=tokens[0], noheader=False, ncol=DEFAULT_VALUE_START_COLUMN_INDEX, delimiter=delimiter)
					else:
						print '[EROR] --> [PARM] Unrecognized expression file extension.'
						return False
			else:
				print '[EROR] --> [PARM] Expression file is required.'
				return False
			if self.biomarkers == None:
				print '[EROR] --> [PARM] At least one biomarker is required.'
				return False
			else:
				self.biomarkers = self.biomarkers.strip('[]\"').split(',')
			if self.mica_labels != None:
				if not os.path.exists(self.mica_labels):
					print '[EROR] --> [PARM] MICA label file does not exists, ignored.'
					self.mica_labels = None
					return False
			else:
				print '[EROR] --> [PARM] MICA label file is required.'
				return False
		return True

	def print_parameters(self):
		print '[PRNT] --> [PARM] Running mode:\t\t' + str(self.mode)
		print '[PRNT] --> [PARM] Input file:\t\t' + str(self.infile)
		print '[PRNT] --> [PARM] Job name:\t\t' + str(self.job_name)
		#print '[PRNT] --> [PARM] Host option:\t\t' + str(self.host)
		if self.mode == 'validate' or self.mode == 'b-quality' or self.mode == 'c-quality':
			print '[PRNT] --> [PARM] True label file:\t\t' + str(self.true_labels)
		if self.mode == 'validate':
			print '[PRNT] --> [PARM] MICA label file:\t\t' + str(self.mica_labels)
		if self.mode == 'overlay':
			print '[PRNT] --> [PARM] MICA label file:\t\t' + str(self.mica_labels)
			print '[PRNT] --> [PARM] Expression file:\t\t' + str(self.expression)
			print '[PRNT] --> [PARM] Biomarkers:\t\t' + str(self.biomarkers)

	@staticmethod
	def print_help():
		usage = """-m OR --mode\t\tRunning mode of the application, default: """ + str(DEFAULT_MODE) + """.
\tPossible options: [clustering, validate, b-quality, c-quality, opt-k, overlay, opt-k-2]
-i OR --infile\t\t[Input File]\tPath to the input file, must be symmetric simmilarity matrix (required)
\tAdditional options can be added using '?' after the input file and putting each additional option into '{name=value}' format
\tAvailable additional options:
\t\tdmin\t\t[Number >= 2], minimum eigen vector dimension, default: """ + str(DEFAULT_DMIN) + """.
\t\t\t\t\tclustering, b-quality and c-quality mode: Used as the minimum number of components to do consensus clustering.
\t\t\t\t\tvalidate mode: Ignored.
\t\tdmax\t\t[Number <= input file dimension], maximum eigen vector dimension, default: """ + str(DEFAULT_DMAX) + """.
\t\t\t\t\tclustering, b-quality and c-quality mode: Used as the maximum number of components to do consensus clustering.
\t\t\t\t\tvalidate mode: Ignored.
\t\tkn\t\t[Range or list of numbers], number of clusters, default: \"""" + str(DEFAULT_KN) + """\". (ignored in validate running mode)
\t\tB\t\t[Number], number of bootstraps for k-means clustering, default: """ + str(DEFAULT_B) + """.
\t\t\t\t\tclustering and c-quality mode: Used as the number of k-means clustering.
\t\t\t\t\tb-quality mode: Used as the number of iterations.
\t\t\t\t\tvalidate and opt-k mode: Ignored.
\t\tdelimiter\t[\",\" or "\\t"], determines the delimiter in the txt input file.
\tExample:
\t\t--infile /path_to_file/infile.txt?{delimiter=","}{kn="[3,4]"}{dmin=18}{dmax=19}{B=5}
-d OR --decomposition\t[Decomposition Method]\tName of the decomposition method for dimension reduction: [PCA, MDS, LPL, LPCA, LPCA2], default: """ + str(DEFAULT_DECOMPOSITION_METHOD) + """. (optional)
-J OR --Job-Name\t[Job Name]\tName that is used for all jobs related to the process. (optional)
-T OR --true-labels\t[Path to true label file], must be in the tab-separated format. (required if mode is validate, b-quality, or c-quality)
-L OR --MICA-labels\tPath to MICA consensus clustering label file. (required if mode is validate or overlay)
-e OR --expression\tExpression matrix path. (required if mode overlay)
\tAdditional options can added using '?' after the expression file and putting each additional option into '{name=value}' format
\tAvailable additional options:
\t\tncol\t\t[Number> 0], determines the first column that includes a value, default: """ + str(DEFAULT_VALUE_START_COLUMN_INDEX) + """.
\t\tdelimiter\t[\",\" or "\\t"], determines the delimiter in the txt expression file.
\tExample:
\t\t--expression /path_to_file/infile.txt?{ncol=3}{delimiter=","}
-bm OR --biomarkers\t\tList of biomarkers which are used to overlay gene expressions on cells. (required if mode overlay)
\tExample:
\t\t--bimarkers \"[gene1, gene2, gene3]\"
"""
		print usage
