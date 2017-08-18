#!/usr/bin/env python
import make_documentation as md
import pyneb as pn

md.run_doxygen('doxy_conf_user.txt')
md.run_doxygen('doxy_conf_devel.txt')

#md.make_latex('pyneb/docs/user_reference_manual_latex', 'dist/PyNeb_{0}_documentation.pdf'.format(pn.__version__))
#md.make_latex('pyneb/docs/devel_reference_manual_latex', 'dist/PyNeb_{0}_documentation_devel.pdf'.format(pyneb.__version__))

#try:
#    md.HandBook_inPDF()
#except:
#    print "HandBook not transformed into PDF"

md.create_zipfile('pyneb/docs/html', '../../../dist/PyNeb_{0}_documentation'.format(pn.__version__))

#md.run('python setup.py sdist upload')

#md.run('python setup.py sdist bdist_wheel --universal')
#md.run('twine upload dist/*')

md.run('python setup.py sdist')

md.run('scp dist/PyNeb-{0}.tar.gz taranis:public_html/PyNeb/dist/PyNeb'.format(pn.__version__))
md.run('scp dist/PyNeb_{0}_documentation.zip taranis:public_html/PyNeb/'.format(pn.__version__))
md.run('scp pyneb/docs/PyNeb_Handbook.pdf taranis:public_html/PyNeb/')
md.run('ssh taranis bin/PyNeb_updates.py')

md.run('svn ci -m "This is the version {0}"'.format(pn.__version__))
