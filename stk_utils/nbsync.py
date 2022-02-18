#! /usr/bin/python3
""" Synchronizes changes between .py UnitTest files and .ipynb Jupyter files (originally for VLASS 1.2 tests).

This script DOES NOT create new .py or new .ipynb files, and any sections added/removed should be manually managed.
"""

import os
from pathlib import Path
import re
import json

whitespace_pattern = re.compile(r"^[ \t]*[\n\r]*$")

class Section():
	""" Code section (either .ipynb or .py)
	
	Parses the raw code to be synchronized. Code/comments will be updated to match indent.

	Format:
	Starts with start_pattern and any extra comments, followed by one or more whitespace lines.
	Ends with one or more whitespace lines, followed by any number of comment lines and an end_pattern line.
	There MUST be a whitespace line after any comments at the beggining and before any comments at the end.

	Example format:
	# %% section start @
	####################

	<comments/code to be synchronized>

	##################
	# %% section end @
	"""
	start_pattern = re.compile(r"[ \t]*# %%[ \t]*(.*?)[ \t]*start.*@")
	end_pattern   = re.compile(r"[ \t]*# %%[ \t]*(.*?)[ \t]*end.*@")

	def __init__(self, raw_code, name, source_file, line_start, line_end):
		self.raw_code = raw_code
		self.name = name
		self.source_file = source_file
		self.line_start = line_start
		self.line_end = line_end

		self.sync_line_start = -1
		self.sync_line_end = -1
		self.start_code = ""
		self.sync_code = ""
		self.end_code = ""
		self.indent = ""
		self._parse_raw_code()

		self.new_sync_code = ""

		if self.sync_line_start > line_end:
			raise RuntimeError(f"Programmer error: failed to identify start of sync code in section {name} at {source_file}:{line_start}")
		if self.sync_line_end > line_end:
			raise RuntimeError(f"Programmer error: failed to identify end of sync code in section {name} at {source_file}:{line_start}")

	def _parse_raw_code(self):
		### Get source lines
		source = self.raw_code
		# merge lines into single string
		if type(source) == list: # as from a notebookd
			source = "".join(source)
		source_lines = split_lines(source)
		# don't include last empty line
		if source_lines[-1] == "":
			source_lines = source_lines[:-1]

		### Local variables
		i = -1 # current line index
		header_lines = []
		body_lines = []
		footer_lines = []

		def add_header_line(line):
			header_lines.append(line)
		def add_body_lines(lines, mark_start=False):
			body_lines.extend(lines)
			if mark_start:
				self.sync_line_start = self.line_start + i
		def add_footer_line(line, mark_end=False):
			footer_lines.append(line)
			if mark_end:
				self.sync_line_end = self.line_start + i - 1

		### Parsing state machine
		state = 'HEADER'
		is_last_whitespace = False
		for line in source_lines:
			i += 1
			line_s = line.lstrip()
			is_whitespace = whitespace_pattern.match(line_s) != None
			is_comment = (line_s != "") and (line_s[0] == "#")

			if state == 'HEADER':
				if is_whitespace: # whitespace line, next non-whitespace line after this starts the body
					add_header_line(line)
					state = 'HEADER_WHITESPACE'
				elif is_comment: # include leading comment lines in header
					add_header_line(line)
				else:
					raise RuntimeError(f"Expected a whitespace line in section at {self.source_file}:{self.line_start} before code started at line {self.line_start+i}")

			elif state == 'HEADER_WHITESPACE':
				if is_whitespace: # extra whitespace line, include in header
					add_header_line(line)
				else: # anything else starts the body
					self.indent = line.replace(line_s, "")
					add_body_lines([line], mark_start=True)
					state = 'BODY'

			elif state == 'BODY':
				if is_whitespace: # whitespace line, possibly the start of the footer
					add_footer_line(line, mark_end=True)
					state = 'FOOTER'
				else:
					add_body_lines([line])

			elif state == 'FOOTER':
				if is_whitespace: # include extra whitespace lines before the end
					if is_last_whitespace: # continuation of whitespace lines before the end
						add_footer_line(line)
					else: # last line must have been a comment line
						add_body_lines(footer_lines)
						footer_lines.clear()
						add_footer_line(line, mark_end=True)
				elif is_comment:
					add_footer_line(line)
				else:
					add_body_lines(footer_lines)
					footer_lines.clear()
					add_body_lines([line])
					state = 'BODY'

			is_last_whitespace = is_whitespace

		if state in ['HEADER','HEADER_WHITESPACE']:
			raise RuntimeError(f"Didn't find any code in section at {self.source_file}:{self.line_start}")
		if state == 'BODY':
			raise RuntimeError(f"Didn't find the end of the section at {self.source_file}:{self.line_start}. Make sure there's an empty line before the section end.")
		if len(body_lines) == 0:
			raise RuntimeError(f"Didn't find any code in the section at {self.source_file}:{self.line_start}")
		if len(footer_lines) == 0:
			raise RuntimeError(f"Expected a whitespace line in section at {self.source_file}:{self.line_start} before reaching end of section at line {self.line_start+i-1}")

		# unindent sync code lines
		body_lines = [line.replace(self.indent, "", 1) for line in body_lines]

		self.start_code = "\n".join(header_lines) + "\n"
		self.sync_code =  "\n".join(body_lines)   + "\n"
		self.end_code =   "\n".join(footer_lines) + "\n"

	def get_updated_code(self, indent_whitespace_lines=False):
		indented_sync_lines = indent_lines(split_lines(self.new_sync_code), self.indent, indent_whitespace_lines)
		indented_sync_code  = "\n".join(indented_sync_lines)
		return self.start_code + indented_sync_code + self.end_code

	def set_new_sync_code(self, new_sync_code):
		self.new_sync_code = new_sync_code

class NotebookSection(Section):
	def set_nb_vals(self, cell_idx):
		self.cell_idx = cell_idx

def split_lines(txt):
	return txt.replace("\r\n", "\n").replace("\n\r", "\n").split("\n")

def indent_lines(lines, indent, indent_whitespace_lines):
	for i in range(len(lines)):
		line = lines[i]
		if (whitespace_pattern.match(line) != None) and (not indent_whitespace_lines):
			pass # don't indent whitespace-only lines
		else:
			lines[i] = indent + line
	return lines

def find_files():
	currdir = Path(os.getcwd())
	basedir = Path(currdir).parent.absolute()

	# find the only recognized unittest .py file
	ut_name = str(Path(basedir, 'stakeholder/nb1', 'test_standard_cube_briggsbwtaper.py'))
	if not os.path.exists(ut_name):
		raise RuntimeError(f"Can't find unittest file {ut_name}")
	ut_files = [ut_name]

	# find jupyter notebook .ipynb files
	nb_names = list(currdir.glob("*.ipynb"))
	if len(nb_names) == 0:
		raise RuntimeError(f"Can't find any Jupyter notebook files at {currdir}")
	nb_files = [str(Path(currdir, nb_name)) for nb_name in nb_names]

	return ut_files, nb_files

def get_sections(file_name, file_or_lines, section_class=None):
	ret = []
	section_class = section_class if (section_class != None) else Section

	i = 0 # current line number
	in_section = False
	section_name = ""
	line_start = -1
	line_end = -1
	raw_code = ""
	for line in file_or_lines:
		i += 1
		start_match = Section.start_pattern.match(line)
		end_match = Section.end_pattern.match(line)

		if start_match != None:
			if not in_section:
				in_section = True
				section_name = start_match[1]
				line_start = i
			else:
				raise RuntimeError(f"Found unexpected section start in section \"{section_name}\" at {file_name}:{i} (previous section starts at line {line_start})")

		if in_section:
			raw_code += line

		if end_match != None:
			if not in_section:
				raise RuntimeError(f"Found unmatched end section at {file_name}:{i}")
			else: # if in_section
				end_name = end_match[1]
				if end_name != section_name:
					raise RuntimeError(f"Found unexpected end section \"{end_name}\" while in section \"{section_name}\" at {file_name}:{i} (section starts at line {line_start})")
				line_end = i			
				ret.append(section_class(raw_code, section_name, file_name, line_start, line_end))
				in_section = False
				section_name = ""
				line_start = -1
				line_end = -1
				raw_code = ""
	if in_section:
		raise RuntimeError(f"Found unmatched start section \"{section_name}\" at {file_name}:{line_start}")

	return ret

def get_ut_sections(ut_files):
	for ut_name in ut_files:
		with open(ut_name, 'r') as fin:
			sections = get_sections(ut_name, fin)
	ret = {}
	for section in sections:
		if section.name in ret:
			raise RuntimeError(f"Found duplicate section \"{section.name}\" in {section.source_file} (original in {ret[section.name].source_file})")
		ret[section.name] = section
	return ret

def get_nb_sections(nb_files):
	sections = []
	for nb_name in nb_files:
		with open(nb_name, 'r') as fin:
			parsed = json.load(fin)
		for cell_idx in range(len(parsed['cells'])):
			cell = parsed['cells'][cell_idx]
			if cell['cell_type'] != 'code':
				continue
			source = cell['source']
			if len(source) > 1:
				source[-1] = source[-1].rstrip()
			cell_sections = get_sections(nb_name, source, section_class=NotebookSection)
			for section in cell_sections:
				section.set_nb_vals(cell_idx=cell_idx)
			sections += cell_sections
	ret = {}
	for section in sections:
		if section.name in ret:
			raise RuntimeError(f"Found duplicate section \"{section.name}\" in {section.source_file} (original in {ret[section.name].source_file})")
		ret[section.name] = section
	return ret

def sync_section(ut_section, nb_section, mode='tonb'):
	if mode == 'tonb':
		nb_section.set_new_sync_code(ut_section.sync_code)
	elif mode == 'tout':
		ut_section.set_new_sync_code(nb_section.sync_code)
	return True

def update_sections_in_file(file_name, sections_to_update, mode='tonb'):
	# update sections from last to first, so that updating one section doesn't affect line indexing of later sections
	sections_to_update.sort(key=lambda s: s.line_start, reverse=True)

	def update_section(section, source_lines, section_lines):
		before_lines = [] if (section.line_start == 1)               else source_lines[:section.line_start-1]
		after_lines  = [] if (section.line_end >= len(source_lines)) else source_lines[section.line_end:]
		return before_lines + section_lines + after_lines

	if mode == 'tonb':
		with open(file_name, 'r') as fin:
			parsed = json.load(fin)
		for section in sections_to_update:
			section_lines = split_lines(section.get_updated_code(indent_whitespace_lines=True))
			section_lines = section_lines[:-1] # split_lines adds an extra line at the end
			section_lines = [line+'\n' for line in section_lines] # .ipynb json needs extra '\n' characters
			section_lines[-1] = section_lines[-1].rstrip() #        ...except for the last line
			parsed['cells'][section.cell_idx]['source'] = update_section(section, parsed['cells'][section.cell_idx]['source'], section_lines)
		outstr = json.dumps(parsed, indent=2)
		with open(file_name, 'w') as fout:
			fout.write(outstr)
	else:
		with open(file_name, 'r') as fin:
			lines = fin.readlines()
		for section in sections_to_update:
			lines = update_section(section, lines, [section.get_updated_code()])
		with open(file_name, 'w') as fout:
			fout.writelines(lines)

if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(description='Synchronizes .py UnitTest and .ipynb Jupyter files (originally for VLASS 1.2 tests)')
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument('--tonb', action='store_const', const='tonb', dest='mode', help='to NoteBook (from .py to .ipynb)')
	group.add_argument('--tout', action='store_const', const='tout', dest='mode', help='to UnitTest (from .ipynb to .py)')
	parser.add_argument('--dryrun', action='store_true', help='Dry run, don\'t modify any files')
	parser.add_argument('--file_stem', action='store', help='File stem of notebook/script to sync.')
	parser.add_argument('--verbose', '-v', action='count', default=0)

	args = parser.parse_args()
	
	# catalog all sections
	ut_files, nb_files = find_files()
	ut_sections = get_ut_sections(ut_files)
	nb_sections = get_nb_sections(nb_files)
	if (args.verbose >= 2):
		for sections in [nb_sections]:#[ut_sections, nb_sections]:
			for section in sections.values():
				print(f"Section \"{section.name}\" [{section.source_file}]:")
				print(f"Raw code range: {section.line_start}-{section.line_end}")
				print(f"Sync code range: {section.sync_line_start}-{section.sync_line_end}")
				print(f"Indent: \"{section.indent}\"")
				print(f">>>\n{section.sync_code}<<<")
				print("")
	for k in ut_sections.keys():
		if k not in nb_sections:
			print(f"Warning: UnitTest section \"{k}\" was not found in any notebooks")
	for k in nb_sections.keys():
		if k not in ut_sections:
			print(f"Warning: Notebook section \"{k}\" was not found in any unittests")

	# find the differing sections
	differing_section_names = []
	differing_sections_by_files = {}
	for k in ut_sections.keys():
		if k not in nb_sections:
			continue
		ut_section = ut_sections[k]
		nb_section = nb_sections[k]
		if ut_section.sync_code != nb_section.sync_code:
			differing_section_names.append(k)
			to_section = nb_section if (args.mode == 'tonb') else ut_section
			file_name = to_section.source_file
			if file_name not in differing_sections_by_files:
				differing_sections_by_files[file_name] = []
			differing_sections_by_files[file_name].append(to_section)

	# update the sections
	for k in differing_section_names:
		arrow = '-->' if (args.mode == 'tonb') else '<--'
		print(f"{ut_sections[k].source_file} {arrow} {nb_sections[k].source_file} [{nb_sections[k].name}]")
		sync_section(ut_sections[k], nb_sections[k], args.mode)

	# save the sections back out
	if not args.dryrun:
		for file_name in differing_sections_by_files.keys():
			update_sections_in_file(file_name, differing_sections_by_files[file_name], args.mode)