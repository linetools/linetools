def get_package_data():
    return {
        'linetools.tests': ['files/*.fits', 'files/*.dat', 'files/*.json', 'files/*.all', 'files/*.fits.gz', 'files/*.out'],
        _ASTROPY_PACKAGE_NAME_ + '.tests': ['coveragerc']}
