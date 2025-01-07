cd python
rm -rf build dist *.egg-info
python setup.py bdist_wheel
pip uninstall glassppy -y
cd dist
ls | xargs pip install
