# Added by Sara Saray and Shailesh Appukuttan to run validation tests and register the results in the HBP Validation Framework

import os
from hippounit.utils import ModelLoader
from neuron import h
import numpy

class Poirazi_2003_CA1(ModelLoader):
    # model UUID on HBP Model Catalog
    model_uuid = "d49c98fd-745f-4478-b3c6-4173c28b2387"
    model_version = "1.0"

    def __init__(self):
        # path to mod files
        model_path = os.path.dirname(os.path.abspath(__file__))
        print "model_path = ", model_path
        mod_files_path = model_path + "/mechanism/"

        #Load cell model
        ModelLoader.__init__(self, mod_files_path = mod_files_path)

	super(Poirazi_2003_CA1, self).__init__(mod_files_path=mod_files_path)

        # outputs will be saved in subfolders named like this:
        self.name="Poirazi_et_al_2003"

        # path to hoc file
        # the model must not display any GUI!!
        self.hocpath = model_path + "/experiment/cluster-dispersion/main_model.hoc"

        # If the hoc file doesn't contain a template, this must be None (the default value is None)
        self.template_name = None

        # model.SomaSecList_name should be None, if there is no Section List in the model for the soma, or if the name of the soma section is given by setting model.soma (the default value is None)
        self.SomaSecList_name = None
        # if the soma is not in a section list or to use a specific somatic section, add its name here:
        self.soma = 'soma[0]'

        # For the PSP Attenuation Test, and Back-propagating AP Test a section list containing the trunk sections is needed
        self.TrunkSecList_name = 'apical_trunk_list'
        # For the Oblique Integration Test a section list containing the oblique dendritic sections is needed
        self.ObliqueSecList_name = 'oblique_dendrites'

        self.v_init = -70
        self.celsius = 34

        #setting synapse parameters
        self.AMPA_name = 'GLU'
        self.NMDA_name = 'NMDA' 
        self.AMPA_NMDA_ratio = 1/2.5

    def set_multiple_ampa_nmda(self, dend_loc, number, AMPA_weight):
        """Used in ObliqueIntegrationTest"""

        ndend, xloc, loc_type = dend_loc

        h("Deadtime_GLU = 0.025")
        h("Deadtime_NMDA = 0.025")

        exec("self.dendrite=h." + ndend)

        for i in range(number):

            if self.AMPA_name: # if this is given, the AMPA model defined in a mod file is used, else the built in Exp2Syn
                exec("self.ampa_list[i] = h."+self.AMPA_name+"(xloc, sec=self.dendrite)")

                self.ampa_list[i].gmax = AMPA_weight
                # self.ampa_list[i].Deadtime = 0.025
            else: 
                self.ampa_list[i] = h.Exp2Syn(xloc, sec=self.dendrite)
                self.ampa_list[i].tau1 = self.AMPA_tau1
                self.ampa_list[i].tau2 = self.AMPA_tau2
                print 'The built in Exp2Syn is used as the AMPA component. Tau1 = ', self.AMPA_tau1, ', Tau2 = self.AMPA_tau2.'

            exec("self.nmda_list[i] = h."+self.NMDA_name+"(xloc, sec=self.dendrite)")
            self.nmda_list[i].gmax = AMPA_weight/self.AMPA_NMDA_ratio
            # self.nmda_list[i].Deadtime = 0.025

        self.ndend = ndend
        self.xloc = xloc

    def set_pointers(self, interval, number, AMPA_weight):
        """Used in ObliqueIntegrationTest"""

        self.fakecell = h.Section(name='fakecell')

        for i in range(number):
            self.epsp_ic[i] = h.IClamp(self.fakecell(0.5))
            self.epsp_ic[i].amp = 1
            self.epsp_ic[i].dur = 0.05  # let it be shorter than the interval bw synapse activation 
            self.epsp_ic[i].delay = self.start + (i*interval)
            # print self.ampa_list[i].gmax

            # h.setpointer(self.ampa_list[i].pre, self.epsp_ic[i].i)
            # h.setpointer(self.nmda_list[i].pre, self.epsp_ic[i].i)

            h.setpointer(self.epsp_ic[i]._ref_i, 'pre', self.ampa_list[i])
            h.setpointer(self.epsp_ic[i]._ref_i, 'pre', self.nmda_list[i])
        # print "pointer is set"



    def run_multiple_syn(self, dend_loc, interval, number, weight):
        """Used in ObliqueIntegrationTest"""

        self.ampa_list = [None] * number 
        self.nmda_list = [None] * number 
        self.ns_list = [None] * number

        self.epsp_ic = [None] * number


        self.initialise()

        if self.cvode_active:
            h.cvode_active(1)
        else:
            h.cvode_active(0)

        self.set_multiple_ampa_nmda(dend_loc, number, weight)

        self.set_pointers(interval, number, weight)


        exec("self.sect_loc=h." + str(self.soma)+"("+str(0.5)+")")

        # initiate recording
        rec_t = h.Vector()
        rec_t.record(h._ref_t)

        rec_v = h.Vector()
        rec_v.record(self.sect_loc._ref_v)

        rec_v_dend = h.Vector()
        rec_v_dend.record(self.dendrite(self.xloc)._ref_v)

        # rec_i = h.Vector()
        # rec_i.record(self.epsp_ic[0]._ref_i)

        h.stdinit()

        dt = 0.025
        h.dt = dt
        h.steps_per_ms = 1 / dt
        h.v_init = self.v_init #-80

        h.celsius = self.celsius
        h.init()
        h.tstop =500
        h.run()

        # get recordings
        t = numpy.array(rec_t)
        v = numpy.array(rec_v)
        v_dend = numpy.array(rec_v_dend)
        # i=numpy.array(rec_i)

        # print 'vdend max', numpy.max(v_dend)
        # print 'i max', numpy.max(i)
        # print 'i min', numpy.min(i)
        return t, v, v_dend
