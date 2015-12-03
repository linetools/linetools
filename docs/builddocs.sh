rm -rf ../build/lib.macosx-10.5-x86_64-3.4/
rm -rf _build/html
cd ..
python setup.py build_sphinx -w
cd docs
