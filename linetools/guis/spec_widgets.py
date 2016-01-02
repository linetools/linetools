""" Module for spec widgets
"""
from __future__ import print_function, absolute_import, division, unicode_literals

from PyQt4 import QtGui
from PyQt4 import QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from linetools.guis import utils as ltgu


class ExamineSpecWidget(QtGui.QWidget):
    """ Widget to plot a spectrum and interactively
        fiddle about.  Akin to XIDL/x_specplot.pro

        12-Dec-2014 by JXP
        Parameters:
        ------------
        key_events : bool, optional
          Use key events? [True]
          Useful when coupling to other widgets
    """
    def __init__(self, ispec, parent=None, status=None, llist=None,
                 abs_sys=None, norm=True, second_file=None, zsys=None,
                 key_events=True, vlines=None, plotzero=False, exten=None):
        '''
        spec = Spectrum1D
        '''
        super(ExamineSpecWidget, self).__init__(parent)

        # Spectrum
        spec, spec_fil = ltgu.read_spec(ispec, exten=exten)
        self.orig_spec = spec  # For smoothing
        self.spec = self.orig_spec

        self.vlines = []
        if vlines is not None:
            self.vlines.extend(vlines)

        self.plotzero = plotzero

        # Other bits (modified by other widgets)
        self.continuum = None
        self.model = None
        self.bad_model = None  # Discrepant pixels in model
        self.use_event = 1

        # Abs Systems
        if abs_sys is None:
            self.abs_sys = []
        else:
            self.abs_sys = abs_sys
        self.norm = norm
        self.psdict = {}  # Dict for spectra plotting
        self.adict = {}  # Dict for analysis
        self.init_spec()
        self.xval = None  # Used with velplt

        # Status Bar?
        if not status is None:
            self.statusBar = status

        # Line List?
        if llist is None:
            self.llist = {'Plot': False, 'List': 'None', 'z': 0., 'Lists': []}
        else:
            self.llist = llist

        # zsys
        if not zsys is None:
            self.llist['z'] = zsys

        # Create the mpl Figure and FigCanvas objects.
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 150 # 150
        self.fig = Figure((8.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)

        self.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.canvas.setFocus()
        if key_events:
            self.canvas.mpl_connect('key_press_event', self.on_key)
        self.canvas.mpl_connect('button_press_event', self.on_click)

        # Make two plots
        self.ax = self.fig.add_subplot(1,1,1)
        self.fig.subplots_adjust(hspace=0.1, wspace=0.1)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.canvas)

        self.setLayout(vbox)

        # Draw on init
        self.on_draw()

    # Setup the spectrum plotting info
    def init_spec(self):
        #xy min/max
        xmin = np.min(self.spec.dispersion.value)
        xmax = np.max(self.spec.dispersion.value)
        from linetools.spectra.plotting import get_flux_plotrange
        ymin, ymax = get_flux_plotrange(self.spec.flux.value)#, mult_pos=1)
        # ymed = np.median(self.spec.flux.value)
        # ymin = 0. - 0.1*ymed
        # ymax = ymed * 1.5
        #
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        self.psdict['xmnx'] = np.array([xmin,xmax])
        self.psdict['ymnx'] = [ymin,ymax]
        self.psdict['sv_xy'] = [ [xmin,xmax], [ymin,ymax] ]
        self.psdict['tmp_xy'] = None
        self.psdict['nav'] = xxgu.navigate(0,0,init=True)
        # Analysis dict
        self.adict['flg'] = 0 # Column density flag


    # Main Driver
    def on_key(self,event):

        flg = -1

        ## NAVIGATING
        if event.key in self.psdict['nav']:
            flg = xxgu.navigate(self.psdict,event)

        ## DOUBLETS
        if event.key in ['C','M','X','4','8','B']:  # Set left
            wave = xxgu.set_doublet(self, event)
            #print('wave = {:g},{:g}'.format(wave[0], wave[1]))
            self.ax.plot( [wave[0],wave[0]], self.psdict['ymnx'], '--', color='red')
            self.ax.plot( [wave[1],wave[1]], self.psdict['ymnx'], '--', color='red')
            flg = 2 # Layer

        ## SMOOTH
        if event.key == 'S':
            self.spec = self.spec.box_smooth(2)
            flg = 1
        if event.key == 'U':
            self.spec = self.orig_spec
            flg = 1

        ## Lya Profiles
        if event.key in ['D', 'R']:
            # Set NHI
            if event.key == 'D':
                NHI = 20.3
            elif event.key == 'R':
                NHI = 19.0
            zlya = event.xdata/1215.6701 - 1.
            self.llist['z'] = zlya
            # Generate Lya profile
            lya_line = AbsLine(1215.6701*u.AA)
            lya_line.z = zlya
            lya_line.attrib['N'] = NHI
            lya_line.attrib['b'] = 30.
            self.lya_line = xspec.voigt.voigt_model(self.spec.dispersion, lya_line, Npix=3.)
            self.adict['flg'] = 4
            flg = 1

        ## ANALYSIS:  EW, AODM column density
        if event.key in ['N', 'E', '$']:
            # If column check for line list
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            if (event.key in ['N','E']) & (self.llist['List'] == 'None'):
                print('xspec: Choose a Line list first!')
                try:
                    self.statusBar().showMessage('Choose a Line list first!')
                except AttributeError:
                    pass
                self.adict['flg'] = 0
                return
            flg = 1

            if self.adict['flg'] == 0:
                self.adict['wv_1'] = event.xdata # wavelength
                self.adict['C_1'] = event.ydata # continuum
                self.adict['flg'] = 1 # Plot dot
            else:
                self.adict['wv_2'] = event.xdata # wavelength
                self.adict['C_2'] = event.ydata # continuum
                self.adict['flg'] = 2 # Ready to plot + print

                # Sort em + make arrays
                iwv = np.array(sorted([self.adict['wv_1'], self.adict['wv_2']])) * self.spec.wcs.unit
                ic = np.array(sorted([self.adict['C_1'], self.adict['C_2']]))

                # Calculate the continuum (linear fit)
                param = np.polyfit(iwv, ic, 1)
                cfunc = np.poly1d(param)
                self.spec.conti = cfunc(self.spec.dispersion)

                if event.key == '$': # Simple stats
                    pix = self.spec.pix_minmax(iwv)[0]
                    mean = np.mean(self.spec.flux[pix])
                    median = np.median(self.spec.flux[pix])
                    stdv = np.std(self.spec.flux[pix]-self.spec.conti[pix])
                    S2N = median / stdv
                    mssg = 'Mean={:g}, Median={:g}, S/N={:g}'.format(mean,median,S2N)
                else:
                    # Find the spectral line (or request it!)
                    rng_wrest = iwv / (self.llist['z']+1)
                    gdl = np.where( (self.llist[self.llist['List']].wrest-rng_wrest[0]) *
                                    (self.llist[self.llist['List']].wrest-rng_wrest[1]) < 0.)[0]
                    if len(gdl) == 1:
                        wrest = self.llist[self.llist['List']].wrest[gdl[0]]
                    else:
                        if len(gdl) == 0: # Search through them all
                            gdl = np.arange(len(self.llist[self.llist['List']]))
                        sel_widg = SelectLineWidget(self.llist[self.llist['List']]._data[gdl])
                        sel_widg.exec_()
                        line = sel_widg.line
                        #wrest = float(line.split('::')[1].lstrip())
                        quant = line.split('::')[1].lstrip()
                        spltw = quant.split(' ')
                        wrest = Quantity(float(spltw[0]), unit=spltw[1])
                    # Units
                    if not hasattr(wrest,'unit'):
                        # Assume Ang
                        wrest = wrest * u.AA

                    # Generate the Spectral Line
                    aline = AbsLine(wrest,linelist=self.llist[self.llist['List']])
                    aline.attrib['z'] = self.llist['z']
                    aline.analy['spec'] = self.spec

                    # AODM
                    if event.key == 'N':
                        # Calculate the velocity limits and load-up
                        aline.analy['vlim'] = const.c.to('km/s') * (
                            ( iwv/(1+self.llist['z']) - wrest) / wrest )

                        # AODM
                        #QtCore.pyqtRemoveInputHook()
                        #xdb.set_trace()
                        #QtCore.pyqtRestoreInputHook()
                        aline.measure_aodm()
                        mssg = 'Using '+ aline.__repr__()
                        mssg = mssg + ' ::  logN = {:g} +/- {:g}'.format(
                            aline.attrib['logN'], aline.attrib['sig_logN'])
                    elif event.key == 'E':  #EW
                        aline.analy['wvlim'] = iwv
                        aline.measure_restew()
                        mssg = 'Using '+ aline.__repr__()
                        mssg = mssg + ' ::  Rest EW = {:g} +/- {:g}'.format(
                            aline.attrib['EW'].to(mAA), aline.attrib['sig_EW'].to(mAA))
                # Display values
                try:
                    self.statusBar().showMessage(mssg)
                except AttributeError:
                    pass
                print(mssg)

                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()


        ## Velocity plot
        if event.key == 'v':
            flg = 0
            from xastropy.xguis import spec_guis as xsgui
            z=self.llist['z']
            # Check for a match in existing list and use it if so
            if len(self.abs_sys) > 0:
                zabs = np.array([abs_sys.zabs for abs_sys in self.abs_sys])
                mt = np.where( np.abs(zabs-z) < 1e-4)[0]
            else:
                mt = []
            if len(mt) == 1:
                ini_abs_sys = self.abs_sys[mt[0]]
                outfil = ini_abs_sys.absid_file
                self.vplt_flg = 0 # Old one
                print('Using existing ID file {:s}'.format(outfil))
            else:
                ini_abs_sys = None
                outfil = None
                if self.llist['List'] == 'None':
                    print('Need to set a line list first!!')
                    self.vplt_flg = -1 # Nothing to do here
                    return
                self.vplt_flg = 1 # New one

            # Outfil
            if outfil is None:
                i0 = self.spec.filename.rfind('/')
                i1 = self.spec.filename.rfind('.')
                if i0 < 0:
                    path = './ID_LINES/'
                else:
                    path = self.spec.filename[0:i0]+'/ID_LINES/'
                outfil = path + self.spec.filename[i0+1:i1]+'_z'+'{:.4f}'.format(z)+'_id.fits'
                xutils.files.ensure_dir(outfil)
                self.outfil = outfil
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()

            # Launch
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            gui = xsgui.XVelPltGui(self.spec, z=z, outfil=outfil, llist=self.llist,
                                   abs_sys=ini_abs_sys, norm=self.norm,
                                   sel_wv=self.xval*self.spec.wcs.unit)
            gui.exec_()
            if gui.flg_quit == 0: # Quit without saving (i.e. discarded)
                self.vplt_flg = 0
            else:
                # Push to Abs_Sys
                if len(mt) == 1:
                    self.abs_sys[mt[0]] = gui.abs_sys
                else:
                    self.abs_sys.append(gui.abs_sys)
                    print('Adding new abs system')
            # Redraw
            flg=1

        # Dummy keys
        if event.key in ['shift', 'control', 'shift+super', 'super+shift']:
            flg = 0

        # Draw
        if flg==1: # Default is not to redraw
            self.on_draw()
        elif flg==2: # Layer (no clear)
            self.on_draw(replot=False)
        elif flg==-1: # Layer (no clear)
            try:
                self.statusBar().showMessage('Not a valid key!  {:s}'.format(event.key))
            except AttributeError:
                pass

    # Click of main mouse button
    def on_click(self,event):
        try:
            print('button={:d}, x={:f}, y={:f}, xdata={:f}, ydata={:f}'.format(
                event.button, event.x, event.y, event.xdata, event.ydata))
        except ValueError:
            print('Out of bounds')
            return
        if event.button == 1: # Draw line
            self.xval = event.xdata
            self.ax.plot( [event.xdata,event.xdata], self.psdict['ymnx'], ':', color='green')
            self.on_draw(replot=False)

            # Print values
            try:
                self.statusBar().showMessage('x,y = {:f}, {:f}'.format(event.xdata,event.ydata))
            except AttributeError:
                return

    # ######
    def on_draw(self, replot=True, no_draw=False):
        """ Redraws the spectrum
        no_draw: bool, optional
          Draw the screen on the canvas?
        """
        #
        if replot is True:
            self.ax.clear()
            self.ax.plot(self.spec.dispersion, self.spec.flux,
                'k-',drawstyle='steps-mid')
            try:
                self.ax.plot(self.spec.dispersion, self.spec.sig, 'r:')
            except ValueError:
                pass
            self.ax.set_xlabel('Wavelength')
            self.ax.set_ylabel('Flux')

            # Continuum?
            if self.continuum is not None:
                self.ax.plot(self.continuum.dispersion, self.continuum.flux,
                    color='purple')

            # Model?
            if self.model is not None:
                self.ax.plot(self.model.dispersion, self.model.flux,
                    color='cyan')
                if self.bad_model is not None:
                    self.ax.scatter(self.model.dispersion[self.bad_model],
                        self.model.flux[self.bad_model],  marker='o',
                        color='red', s=3.)


            # Spectral lines?
            if self.llist['Plot'] is True:
                ylbl = self.psdict['ymnx'][1]-0.2*(self.psdict['ymnx'][1]-self.psdict['ymnx'][0])
                z = self.llist['z']
                wvobs = np.array((1+z) * self.llist[self.llist['List']].wrest)
                gdwv = np.where( (wvobs > self.psdict['xmnx'][0]) &
                                 (wvobs < self.psdict['xmnx'][1]))[0]
                for kk in range(len(gdwv)):
                    jj = gdwv[kk]
                    wrest = self.llist[self.llist['List']].wrest[jj].value
                    lbl = self.llist[self.llist['List']].name[jj]
                    # Plot
                    self.ax.plot(wrest*np.array([z+1,z+1]), self.psdict['ymnx'], 'b--')
                    # Label
                    self.ax.text(wrest*(z+1), ylbl, lbl, color='blue', rotation=90., size='small')

            # Abs Sys?
            if not self.abs_sys is None:
                ylbl = self.psdict['ymnx'][0]+0.2*(self.psdict['ymnx'][1]-self.psdict['ymnx'][0])
                clrs = ['red', 'green', 'cyan', 'orange', 'gray', 'purple']*10
                ii=-1
                for abs_sys in self.abs_sys:
                    ii+=1
                    lines = abs_sys.list_of_abslines()
                    #QtCore.pyqtRemoveInputHook()
                    #xdb.set_trace()
                    #QtCore.pyqtRestoreInputHook()
                    wrest = Quantity([line.wrest for line in lines])
                    wvobs = wrest * (abs_sys.zabs+1)
                    gdwv = np.where( ((wvobs.value+5) > self.psdict['xmnx'][0]) &  # Buffer for region
                                    ((wvobs.value-5) < self.psdict['xmnx'][1]))[0]
                    #for kk in range(len(gdwv)):
                    for jj in gdwv:
                        if lines[jj].analy['do_analysis'] == 0:
                            continue
                        # Paint spectrum red
                        wvlim = wvobs[jj]*(1 + lines[jj].analy['vlim']/const.c.to('km/s'))
                        pix = np.where( (self.spec.dispersion > wvlim[0]) & (self.spec.dispersion < wvlim[1]))[0]
                        self.ax.plot(self.spec.dispersion[pix], self.spec.flux[pix], '-',drawstyle='steps-mid',
                                     color=clrs[ii])
                        # Label
                        lbl = lines[jj].analy['name']+' z={:g}'.format(abs_sys.zabs)
                        self.ax.text(wvobs[jj].value, ylbl, lbl, color=clrs[ii], rotation=90., size='x-small')
            # Analysis? EW, Column
            if self.adict['flg'] == 1:
                self.ax.plot(self.adict['wv_1'], self.adict['C_1'], 'go')
            elif self.adict['flg'] == 2:
                self.ax.plot([self.adict['wv_1'], self.adict['wv_2']],
                             [self.adict['C_1'], self.adict['C_2']], 'g--', marker='o')
                self.adict['flg'] = 0
            # Lya line?
            if self.adict['flg'] == 4:
                self.ax.plot(self.spec.dispersion,
                    self.lya_line.flux, color='green')

        # Reset window limits
        self.ax.set_xlim(self.psdict['xmnx'])
        self.ax.set_ylim(self.psdict['ymnx'])

        if self.plotzero:
            self.ax.axhline(0, lw=0.3, color='k')

        for line in self.vlines:
            self.ax.axvline(line, color='k', ls=':')

        # Draw
        if not no_draw:
            self.canvas.draw()

    # Notes on usage
    def help_notes():
        doublets = [ 'Doublets --------',
                     'C: CIV',
                     'M: MgII',
                     'O: OVI',
                     '8: NeVIII',
                     'B: Lyb/Lya'
                     ]
        analysis = [ 'Analysis --------',
                     'N/N: Column density (AODM)',
                     'E/E: EW (boxcar)',
                     '$/$: stats on spectrum'
                     ]

