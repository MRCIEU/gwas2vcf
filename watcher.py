#!/usr/bin/env python

import os
import argparse
import subprocess

parser = argparse.ArgumentParser(description='Orchestrating pipeline.')
parser.add_argument('--directory', help='Directory in which to search for new files')
parser.add_argument('--required-files', help='List of required input files', nargs="+")
parser.add_argument('--input-flag', help='Flag to indicate whether to run. Deleted upon successful execution.', required=True)
parser.add_argument('--output-flag', help='Flag to create upon successful execution', required=False)
parser.add_argument('--command', help='Command to run')
parser.add_argument('--threads', type=int, help='FHow many threads', default=1)

args = parser.parse_args()
print(args)


# This substitutes {} for the id in a command string
# Then deletes the input flag and runs the command
# If command is successful it deletes the input flag and creats the output flag if it is specified
# If unsuccessful then recreates the input flag
def run_command(i):
	os.remove(os.path.join(args.directory, i, args.input_flag))
	cmd = args.command.replace('{}', i)
	print(cmd)
	exitcode = subprocess.call(cmd, shell=True)
	if exitcode == 0:
		if args.output_flag is not None:
			open(os.path.join(args.directory, i, args.output_flag), 'a').close()
	else:
		open(os.path.join(args.directory, i, args.input_flag), 'a').close()

# Get directories with the input flag
idlist = [ name for name in os.listdir(args.directory) if 
	os.path.isdir(os.path.join(args.directory, name)) and 
	os.path.isfile(os.path.join(args.directory, name, args.input_flag)) ]

if len(idlist) == 0:
	print("Nothing to do")
	exit()


idlistc = idlist
for req in args.required_files:
	for i in idlist:
		if not os.path.isfile(os.path.join(args.directory, i, req)):
			print("Missing: " + os.path.join(args.directory, i, req))
			idlistc.remove(i)

if idlistc is None:
	print("No IDs found")
	exit()

print("Running the following command:")
print(args.command)
print("on the following ids: " + str(idlistc))
print("using " + str(args.threads) + " threads")



# from joblib import Parallel, delayed
# Parallel(n_jobs=args.threads)(delayed(run_command)(i) for i in idlistc)

[run_command(i) for i in idlistc]
