rm -rf ../build/lib.macosx-10.5-x86_64-3.4/
cd ..
python setup.py develop -u
python setup.py develop
cd docs
rm -rf _build/html
make html
