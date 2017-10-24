def get_package_data():
    # Installs the testing data files. Unable to get package_data
    # to deal with a directory hierarchy of files, so just explicitly list.
    return {'linetools.analysis.tests': ['files/*.fits', 'files/*.all', 'files/*.out']}
