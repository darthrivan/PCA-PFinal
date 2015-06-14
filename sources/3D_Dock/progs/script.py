#!/usr/bin/env python
import sys
import re
from subprocess import Popen, CalledProcessError, PIPE
from os import getcwd

try:
	import matplotlib.pyplot as plt
	import pandas as pd
except ImportError:
	noImport = True


def run(*arg, **opt):
	proc = Popen(*arg, **opt)  # if opt is not None else Popen(c)
	retcode = proc.wait()
	if retcode != 0:
		raise CalledProcessError(retcode, proc)
	return proc


class Debugger(object):
	"""docstring for Debugger"""
	def __init__(self, arg):
		super(Debugger, self).__init__()
		self.arg = arg

	@staticmethod
	def colorize(string, color, out=sys.stdout):
		colors = {'red': 91, 'green': 92, 'yellow': 93, 'blue': 94}
		if out.isatty():
			format = '\033[%dm%%s\033[0m' % colors.get(color, 0)
		else:
			format = '%s'
		out.write(format % string)
		out.flush()

	@staticmethod
	def log(msg):
		Debugger.colorize(msg, 'blue')

	@staticmethod
	def success():
		Debugger.colorize("Done\n", 'green')

	@staticmethod
	def error(msg='', fatal=True):
		error_title = "Fatal Error" if fatal else "Error"
		if len(msg) > 0:
			Debugger.colorize(error_title + ": " + msg + "\n", 'red', sys.stderr)
		else:
			Debugger.colorize(error_title + "\n", 'red', sys.stderr)
		if fatal:
			exit(1)


class Program(object):
	"""docstring for Program"""
	def __init__(self, exec_name, path=getcwd(), args=[], input_stream=sys.stdin,
				 output_stream=sys.stdout):
		super(Program, self).__init__()
		self.exec_name = exec_name
		self.path = path
		self.args = args
		self.input_stream = input_stream
		self.output_stream = output_stream

	def get_exec_path(self):
		return self.path + '/' + self.exec_name

	def get_exec_name(self):
		return self.exec_name

	def get_arguments(self):
		return self.args

	def get_argument(self, pos):
		return self.args[pos]

	def get_input_stream(self):
		return self.input_stream

	def get_output_stream(self):
		return self.output_stream

	def set_argument(self, arg, pos=None):
		if pos is None:
			self.args.add(arg)
		else:
			self.args[pos] = arg

	def set_arguments(self, args):
		self.args = args

	def set_input_stream(self, input_stream):
		self.input_stream = input_stream

	def set_output_stream(self, output_stream):
		self.output_stream = output_stream

	def run(self):
		run([self.get_exec_path()] +
			 self.get_arguments(),
			 stdin=self.input_stream,
			 stdout=self.output_stream)


class Compiler(object):
	"""docstring for Compiler"""
	def __init__(self, code, flags=[], from_stdin=False):
		super(Compiler, self).__init__()
		self.code = code
		self.flags = flags
		self.from_stdin = from_stdin

	def compile(self, output_name):
		if self.from_stdin:
			Debugger.log("Compiling from stdin...")
			cmd = Popen(['gcc', '-o', output_name, '-x', 'c', '-'] + self.flags,
						stdin=PIPE)
			cmd.communicate(input=self.code)
			if cmd.returncode:
				Debugger.error("Compilation failed")
		else:
			Debugger.log("Compiling " + self.code + "...")
			try:
				run(['gcc', '-o', output_name, self.code] + self.flags)
			except CalledProcessError:
				Debugger.error("Compilation failed")
		Debugger.success()
		return Program(output_name)

	def set_flag(self, option, value=None):
		if value is None:
			self.flags.append(option)
		else:
			self.flags += [option, value]


class Accounter(object):
	"""docstring for Accounter"""
	def account(self):
		pass

	def get_measures(self, columns):
		pass


class InlineAccounter(Accounter):
	"""docstring for InlineAccounter"""
	def __init__(self, program, num_samples=0):
		super(InlineAccounter, self).__init__()
		self.program = program
		self.num_samples = num_samples
		self.measures = []

	def account(self):
		Debugger.log("Initiating accounting...")
		for i in range(0, self.num_samples):
			try:
				cmd = run([self.program.get_exec_path()] +
						   self.program.get_arguments(),
						   stdin=self.program.get_input_stream(),
						   stdout=self.program.get_output_stream(),
						   stderr=PIPE)
				self.__add_measure(self.__parse(cmd.stderr.read()))
			except CalledProcessError:
				Debugger.error("Accounting failed")
		Debugger.success()

	def __add_measure(self, values):
		self.measures.append(values)

	def __parse(self, str):
		def p(a):
			return (a[0], float(a[1]))
		return dict([p(line.split(':')) for line in str.splitlines() if ':' in line])

	def get_measures(self, columns):
		def f(dic):
			return dict([(k, dic[k]) for k in filter(
				lambda e: e in columns, dic.keys())])
		return map(f, self.measures)


class BintimeAccounter(Accounter):
	"""docstring for BintimeAccounter"""
	def __init__(self, program, num_samples=0):
		super(BintimeAccounter, self).__init__()
		self.program = program
		self.num_samples = num_samples
		self.measures = []

	def account(self):
		Debugger.log("Initiating accounting...")
		for i in range(0, self.num_samples):
			try:
				cmd = run(['/usr/bin/time',
						    # '-o', 'sample.txt', '-a',
						    '-f', '%U,%S,%e,%P,%I,%O,%F,%R,%W',
						    self.program.get_exec_path()
						  ] + self.program.get_arguments(),
						   stdin=self.program.get_input_stream(),
						   stdout=self.program.get_output_stream(),
						   stderr=PIPE)
				self.__add_measure(self.__parse(cmd.stderr.read()))
			except CalledProcessError:
				Debugger.error("Accounting failed")
		Debugger.success()

	def __parse(self, str):
		def p(e):
			return float(e) if e[-1] != '%' else float(e[:-1]) / 100.0
		# return [ p(elem) for elem in str[:-1].split(',') ]
		elem = str[:-1].split(',')
		return { 'User': p(elem[0]), 'System': p(elem[1]),
				 'Elapsed': p(elem[2]), 'CPU': p(elem[3]),
				 'Inputs': p(elem[4]), 'Outputs': p(elem[5]),
				 'MajFaults': p(elem[6]), 'MinFaults': p(elem[7]),
				 'Swaps': p(elem[8]) }

	def __add_measure(self, values):
		self.measures.append(values)

	def get_measures(self, columns=['User', 'System', 'Elapsed', 'CPU',
		                            'Inputs', 'Outputs', 'MajFaults',
		                            'MinFaults', 'Swaps']):
		def f(dic):
			return dict([(k, dic[k]) for k in filter(
				lambda e: e in columns, dic.keys())])
		return map(f, self.measures)


class Versioner(object):
	"""docstring for Versioner"""
	def __init__(self, file_path):
		super(Versioner, self).__init__()
		self.file_path = file_path
		self.load_versions()

	def load_versions(self):
		try:
			cmd = run(['git', 'log', '--follow', self.file_path], stdout=PIPE)
			commits = re.findall("(?<=commit )[0-9a-f]{40}", cmd.stdout.read())
		except CalledProcessError:
			Debugger.error("Could not get commits for " + self.file_path)
		self.versions = [ (com, self.read_from_commit(com)) for com in commits ]
		self.versions.reverse()

	def read_from_commit(self, commit):
		try:
			cmd = run(['git', 'show', commit + ':' + self.file_path], stdout=PIPE)
			return cmd.stdout.read()
		except CalledProcessError:
			Debugger.error("Could not get committed version")

	def get_versions(self):
		return self.versions

	def get_versions_except(self, omitted_versions):
		return filter(lambda v: v[0] not in omitted_versions, self.versions)


class Plotter(object):
	"""docstring for Plotter"""
	def __init__(self):
		super(Plotter, self).__init__()
		self.gnuplot = Popen(['gnuplot'], stdin=PIPE)

	def close(self):
		self.gnuplot.communicate('quit\n')

	def plot(self, data, save=None, xlabel=None, ylabel=None):
		self.__create_data_file(data)
		if save is not None:
			if save[-3:] == 'png':
				self.gnuplot.stdin.write('set terminal png\n')
			elif save[-3:] == 'tex':
				self.gnuplot.stdin.write('set terminal epslatex color\n')
			else:
				self.gnuplot.stdin.write('set terminal %s\n' % save[-3:])
			self.gnuplot.stdin.write('set output "%s"\n' % save)
		if xlabel is not None:
			self.gnuplot.stdin.write('set xlabel "%s"\n' % xlabel)
		if ylabel is not None:
			self.gnuplot.stdin.write('set ylabel "%s"\n' % ylabel)
		self.gnuplot.stdin.write('set style fill solid\n')
		self.gnuplot.stdin.write('plot "data.file" using 2:xtic(1) with histogram\n')
		# self.gnuplot.stdin.write('plot "data.file" using 2:xtic(1) with linespoints\n')

	def __create_data_file(self, data):
		datfile = open('data.file', 'w')
		for (commit, value) in data:
			datfile.write("%s %0.3f\n" % (commit[:5], value))
		datfile.close()


class Tester(object):
	"""docstring for Tester"""
	def __init__(self, original_program, binary=False, input_file=None):
		super(Tester, self).__init__()
		self.original_program = original_program
		self.binary = binary
		self.input_file = input_file
		self.__generate_output()

	def __generate_output(self):
		Debugger.log('Generating test output...')
		output_file = open('/tmp/output_original.out', 'wb' if self.binary else 'w')
		input_file = None
		self.original_program.set_output_stream(output_file)
		if self.input_file:
			input_file = open(self.input_file, 'rb' if self.binary else 'r')
			self.original_program.set_input_stream(input_file)
		try:
			self.original_program.run()
		except CalledProcessError:
			Debugger.error('Failed to generate output for ' +
				self.original_program.get_exec_name())
		if self.input_file:
			input_file.close()
		output_file.close()
		Debugger.success()

	def test(self, program, testFunction):
		Debugger.log('Testing %s...' % program.get_exec_name())
		output_file = open('/tmp/output_tested.out', 'wb' if self.binary else 'w')
		input_file = None
		program.set_output_stream(output_file)
		if self.input_file:
			input_file = open(self.input_file, 'rb' if self.binary else 'r')
			program.set_input_stream(input_file)
		try:
			program.run()
		except CalledProcessError:
			Debugger.error('Failed to generate output for ' +
				program.get_exec_name())
		if self.input_file:
			input_file.close()
		output_file.close()
		return testFunction('/tmp/output_original.out', '/tmp/output_tested.out')
		Debugger.success()


class EqualTester(Tester):
	"""docstring for EqualTester"""
	def test(self, program):
		if self.binary:
			def testfunction(original_path, tested_path):
				run(['cmp', '-s', original_path, tested_path])
		else:
			def testfunction(original_path, tested_path):
				run(['diff', '-q', original_path, tested_path])
		super(EqualTester, self).test(program, testfunction)


if __name__ == "__main__":
	import argparse
	argument_parser = argparse.ArgumentParser(prog=__file__,
		description="Script to load other versions of a given program and "
					"account all of them")
	argument_parser.add_argument('-f', '--file',
		help='Program to compare', required=True)
	argument_parser.add_argument('-i', '--input',
		help='Input file to program', required=False)
	argument_parser.add_argument('-u', '--uncommitted', metavar='FILE',
		help='Add uncommitted programs', nargs='+', default=[])
	argument_parser.add_argument('-t', '--testing', action='store_true',
		help='Enable output testing')
	args = vars(argument_parser.parse_args(sys.argv[1:]))

	# PARAMETERS
	OMIT_COMMITS = ['bae0dade883a99c50ba07842c51ae90c19d23dd6']
	COMPILATION_FLAGS = ['-x', 'none', '-I../../fftw-3.3.4/installation/include',
						 '-O3', '-march=native', '-pg', '-g', '-static',
						 'manipulate_structures.o', 'angles.o', 'coordinates.o',
						 'electrostatics.o', 'grid.o', 'qsort_scores.o',
						 '-L../../fftw-3.3.4/installation/lib', '-lfftw3f', '-lm']
	BINARY_OUTPUT = False
	PROGRAM_ARGUMENTS = ['-static', '../../../inputs/proteins/2pka.parsed',
						 '-mobile', '../../../inputs/proteins/5pti.parsed']
	# INPUT_FILE = 'Makefile'
	# INPUT_FILE = None
	INPUT_FILE = args.get('input')
	EXT = 'png'

	def is_output_correct(original_path, tested_path):
		def parse(line):
			values = line.split(',')
			return {'score': int(values[0]),
					'coords': (int(values[1]), int(values[2]), int(values[3])),
					'angles': (int(values[4]), int(values[5]), int(values[6]))}
		try:
			proc_original_path = open(original_path + '.prc', 'w')
			run(['sed', '-E', '-n', "/^\.?G_DATA/ s/^\.?G_DATA\s+-?[0-9]+\s+0\s+(-?[0-9]+)\s+(-?[0-9]+)\s+(-?[0-9]+)\s+(-?[0-9]+)\s+(-?[0-9]+)\s+(-?[0-9]+)\s+(-?[0-9]+)/\\1,\\2,\\3,\\4,\\5,\\6,\\7/p", original_path], stdout=proc_original_path)
			proc_original_path.close()
			proc_tested_path = open(tested_path + '.prc', 'w')
			run(['sed', '-E', '-n', "/^\.?G_DATA/ s/^\.?G_DATA\s+-?[0-9]+\s+0\s+(-?[0-9]+)\s+(-?[0-9]+)\s+(-?[0-9]+)\s+(-?[0-9]+)\s+(-?[0-9]+)\s+(-?[0-9]+)\s+(-?[0-9]+)/\\1,\\2,\\3,\\4,\\5,\\6,\\7/p", tested_path], stdout=proc_tested_path)
			proc_tested_path.close()
			with open(original_path + '.prc', 'r') as original_output:
				original = [parse(line) for line in original_output]
			with open(tested_path + '.prc', 'r') as tested_output:
				tested = [parse(line) for line in tested_output]
			# TODO compare outputs
			for i in range(0, len(original), 3):
				if original[i]['coords'] != tested[i]['coords'] or \
				   tested[i]['score'] > original[i]['score'] + 2 or \
				   tested[i]['score'] < original[i]['score'] - 2:
					return False
				if original[i + 1]['coords'] != tested[i + 1]['coords'] or \
				   tested[i + 1]['score'] > original[i + 1]['score'] + 2 or \
				   tested[i + 1]['score'] < original[i + 1]['score'] - 2:
					return False
				if original[i + 2]['coords'] != tested[i + 2]['coords'] or \
				   tested[i + 2]['score'] > original[i + 2]['score'] + 2 or \
				   tested[i + 2]['score'] < original[i + 2]['score'] - 2:
					return False
			return True
		except CalledProcessError:
			Debugger.error('Failed to compare outputs', fatal=False)

	def acc(commit, code, test=None):
		compiler  = Compiler(code, flags=COMPILATION_FLAGS, from_stdin=True)
		program   = compiler.compile('tmp_' + commit[:5])
		program.set_arguments(PROGRAM_ARGUMENTS)
		if test is not None:
			try:
				test(program, is_output_correct)
			except CalledProcessError:
				Debugger.error('Commit %s not producing same output' % commit)
		program.set_output_stream(DEVNULL)
		accounter = InlineAccounter(program, 1)
		accounter.account()
		measures  = accounter.get_measures('Elapsed')
		return sum([a['Elapsed'] for a in measures]) / len(measures)

	def acc_with_input(commit, code, test=None):
		compiler  = Compiler(code, flags=COMPILATION_FLAGS, from_stdin=True)
		program   = compiler.compile('tmp_' + commit[:5])
		program.set_arguments(PROGRAM_ARGUMENTS)
		input_file = open(INPUT_FILE, 'rb')
		program.set_input_stream(input_file)
		if test is not None:
			try:
				test(program, is_output_correct)
			except CalledProcessError:
				Debugger.error('Commit %s not producing same output' % commit)
		program.set_output_stream(DEVNULL)
		input_file.seek(0)
		# input_file.close()
		# input_file = open(INPUT_FILE, 'rb')
		program.set_input_stream(input_file)
		accounter = InlineAccounter(program, 5)
		accounter.account()
		measures  = accounter.get_measures('Elapsed')
		input_file.close()
		return sum([a['Elapsed'] for a in measures]) / len(measures)

	from os import devnull
	DEVNULL = open(devnull, 'wb')
	versioner = Versioner(args.get('file'))
	versions = versioner.get_versions_except(OMIT_COMMITS)
	for unc in args.get('uncommitted'):
		code_file = open(unc, 'r')
		versions.append(('uncommitted', code_file.read()))
		code_file.close()
	fun = None
	if args.get('testing'):
		compiler = Compiler(versions[0][1], flags=COMPILATION_FLAGS, from_stdin=True)
		program  = compiler.compile('tmp_test')
		program.set_arguments(PROGRAM_ARGUMENTS)
		tester   = Tester(program, binary=BINARY_OUTPUT, input_file=INPUT_FILE)
		fun = tester.test
	account = acc if INPUT_FILE is None else acc_with_input
	elapseds = [(commit, account(commit, code, fun)) for (commit, code) in versions]
	plotter  = Plotter()
	plotter.plot(elapseds, save='elapseds.' + EXT, xlabel='commit', ylabel='Elapsed Time')
	plotter.close()
	plotter  = Plotter()
	plotter.plot(map(lambda a: (a[0], elapseds[0][1] / a[1]), elapseds),
		save='speedups.' + EXT, xlabel='commit', ylabel='Speed Up')
	for (commit, time) in elapseds:
		print "%s %f\n" % (commit, time)
	DEVNULL.close()
