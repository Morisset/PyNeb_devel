When a new release is planned, perform the following actions:

Prepare the release
=============

* Make the list of all the changes in the CHANGE/v.xx.yy-xx-zz.txt file

Pool request to devel
==============

* Merge the current v.XX.YY branche to the devel one.
* Check that the workflows complete correctly

Pool request to master
===============

Once the devel branch is validated and passed all the tests:

* update the version removing the beta indices
* Merge it to the master branch

Make a new release tab on the github server
==============================


Publish the new version on the pypi server
=============================

* Switch to master branch
* Run the following from the root directory (where dist is) and check the tar file is created in dist:

    python setup.py sdist

* Run the following to upload the tar file to pypi server. **Update the value of the distribution in {0}:**

    twine upload dist/PyNeb-{0}.tar.gz

* Check that the new relase is the current one on pypi server: https://pypi.org/project/PyNeb/

Publish the new release on the Google group
==============================

Start a new branch with the new version
===========================



