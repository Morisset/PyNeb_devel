#!/usr/bin/env python

import pyneb
import os
import subprocess

def run(command):
	try:
		proc = subprocess.Popen(command, shell=True)
		proc.communicate()    
		print('Done "{0}"'.format(command))
	except:
		pyneb.log_.warn('Failed to run "{0}"'.format(command), calling = 'setup_DOCS')
		
def run_doxygen(config_file):
	run("sed 's/PROJECT_NUMBER         = .*/PROJECT_NUMBER         = {1}/1' {0} > {0}.tmp".format(config_file, pyneb.__version__))
        run("mv -f {0}.tmp {0}".format(config_file))
        run("doxygen {0}".format(config_file))

def make_latex(dir_, manual_name):
	run("cd {0}; make".format(dir_))
        run("mv -f {0}/refman.pdf {1}".format(dir_, manual_name))

    
def create_zipfile(upload_dir, zip_name):
        run("cd {0}; zip -r --exclude=*.svn* {1}.zip .".format(upload_dir, zip_name))

def HowTo_inPDF():
	# Obsolete
	pass
	#run("/Applications/LibreOffice.app/Contents/MacOS/soffice --headless --invisible --convert-to pdf --outdir pyneb/docs pyneb/docs/PyNeb_HowTo.docx".format(pyneb.__version__))

def HandBook_inPDF():
	run("cd pyneb/doc ; ipython nbconvert --to latex --template mytemplate --post PDF  PyNeb_Handbook")
	run("cp pyneb/doc/PyNeb_Handbook.pdf pyneb/doc/html")

if __name__ == '__main__':
    run_doxygen('doxy_conf_user.txt')
    run_doxygen('doxy_conf_devel.txt')
    create_zipfile('pyneb/docs/html', '../../../dist/PyNeb_{0}_documentation'.format(pyneb.__version__))

    make_latex('pyneb/docs/user_reference_manual_latex', 'dist/PyNeb_{0}_documentation.pdf'.format(pyneb.__version__))
    make_latex('pyneb/docs/devel_reference_manual_latex', 'dist/PyNeb_{0}_documentation_devel.pdf'.format(pyneb.__version__))
    HandBook_inPDF()
    
