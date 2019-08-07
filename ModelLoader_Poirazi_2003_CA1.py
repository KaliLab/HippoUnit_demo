# 

from hippounit.utils import ModelLoader
from neuron import h
import numpy

class ModelLoader_Poirazi_2003_CA1(ModelLoader):

    def __init__(self, name="model", mod_files_path = None):

        #Load cell model
        ModelLoader.__init__(self, name="model", mod_files_path=mod_files_path)

	super(ModelLoader_Poirazi_2003_CA1, self).__init__(name=name, mod_files_path=mod_files_path)

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
