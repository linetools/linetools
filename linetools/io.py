""" I/O routines useful to linetools
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb

import numpy as np
from astropy import units as u

from linetools.spectralline import EmLine
from linetools.lists.linelist import LineList

def emlines_from_alis_output(alis_file):
    """ Generate a list of emission lines from a standard ALIS output file
    Currently loads

    Parameters
    ----------
    alis_file : str

    Returns
    -------
    emlines : list
      List of EmLine objects
    """
    linelist = LineList('Galaxy')

    def get_value(istr):
        eqpos = istr.find('=')
        return istr[eqpos+1:]

    def strip_end(istr):
        for kk in range(-1, -1*len(istr),-1):
            try:
                _ = int(istr[kk])
            except ValueError:
                pass
            else:
                # Cut
                if kk == -1:
                    cutstr = istr
                else:
                    cutstr = istr[:kk+1]
                return float(cutstr)
        return

    flg_model = False
    flg_error = False
    flg_emiss = False
    flg_sigemiss = False
    line_dict = dict(continuum={}, lines=[], errors=[])
    with open(alis_file) as file_to_search:
        for line in file_to_search:
            # model
            if 'model read' in line:
                flg_model = True
            # Errors
            if 'Errors:' in line:
                flg_error = True
            if flg_model:
                # End
                if 'model end' in line:
                    flg_model = False
                    flg_emiss = False
                if 'emission' in line:
                    flg_emiss = True
            if flg_error:
                # End
                if line[0] != '#':
                    flg_error = False
                    flg_sigemiss = False
                    flg_emiss = False
                if 'emission' in line:
                    flg_sigemiss = True
                    flg_emiss = True
            # lines
            if flg_emiss:
                flg_intflux = True
                if flg_sigemiss:
                    conti_sp = line.split()[1:]
                else:
                    conti_sp = line.split()
                if len(conti_sp) == 0:
                    continue
                # Find specid
                for istr in conti_sp:
                    if 'specid' in istr:
                        specid = int(get_value(istr))
                # Add to dict
                if 'continuum=True' in line:  # Continuum
                    pass
                elif conti_sp[0] in ['gaussian']:  # Line :: Could be others
                    for istr in conti_sp:
                        if 'IntFlux=True' in istr:
                            flg_intflux = True
                        elif 'wave=' in istr:
                            wrest = float(get_value(istr))
                    # Values
                    if flg_sigemiss is False:
                        gauss_fit = {}
                        if flg_intflux:
                            gauss_fit['flux'] = strip_end(conti_sp[1])
                        else:
                            gauss_fit['amp'] = strip_end(conti_sp[1])
                        gauss_fit['z'] = strip_end(conti_sp[2])
                        gauss_fit['sigma'] = strip_end(conti_sp[3])
                        gauss_fit['type'] = 'gaussian'
                        gauss_fit['wrest'] = wrest*u.AA
                        # Init?
                        line_dict['lines'].append(gauss_fit)
                    else:
                        waves = np.array([idict['wrest'].value for idict in line_dict['lines']])
                        mt = np.where(np.abs(waves-wrest) < 1e-3)[0]
                        if len(mt) != 1:
                            raise ValueError("Bad wrest")
                        # Fill in
                        if flg_intflux:
                            line_dict['lines'][mt[0]]['sig_flux'] = strip_end(conti_sp[1])
                        else:
                            line_dict['lines'][mt[0]]['sig_amp'] = strip_end(conti_sp[1])

    # Generate EmLines
    emlines = []
    for iline in line_dict['lines']:
        obj = EmLine(iline['wrest'], z=iline['z'], linelist=linelist, closest=True)
        # Fill
        try:
            obj.attrib['flux'] = iline['flux']*u.erg/u.s
        except KeyError:
            # Amplitude
            pdb.set_trace()
        else:
            obj.attrib['sig_flux'] = iline['sig_flux']*u.erg/u.s
            obj.attrib['flag_flux'] = 1  # Measured
        # Append
        emlines.append(obj)
    # Return
    return emlines










