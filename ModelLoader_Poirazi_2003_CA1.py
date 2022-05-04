# Use this class (instead of the ModelLoader class of hippounit.utils) to test the Poirazi et al. 2003 model using its own receptor models in the Oblique Integration Test. This new version of the synapse functions of the ModelLoader class of HippoUnit can deal with the different (a bit outdated way using pointers) activation of the receptor models (point processes). This child class inherits from the ModelLoader class.

from __future__ import division

from builtins import range
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
                #print 'The built in Exp2Syn is used as the AMPA component. Tau1 = ', self.AMPA_tau1, ', Tau2 = self.AMPA_tau2.'

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
        h.steps_per_ms = 1/dt
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
        
        
    def set_ampa_nmda_pathway(self, dend_loc, pathway, AMPA_weight):
        """Used in PathwayInteractionTest"""

        ndend, xloc = dend_loc
        
        h("Deadtime_GLU = 0.025")
        h("Deadtime_NMDA = 0.025")

        exec("self.dendrite=h." + ndend)


        if self.AMPA_name: # if this is given, the AMPA model defined in a mod file is used, else the built in Exp2Syn
            exec("self.ampa = h."+self.AMPA_name+"(xloc, sec=self.dendrite)")
            self.ampa.gmax = AMPA_weight
        else:
            self.ampa = h.Exp2Syn(xloc, sec=self.dendrite)
            self.ampa.tau1 = self.AMPA_tau1
            self.ampa.tau2 = self.AMPA_tau2
            #print 'The built in Exp2Syn is used as the AMPA component. Tau1 = ', self.AMPA_tau1, ', Tau2 = ', self.AMPA_tau2 , '.'

        if self.NMDA_name: # if this is given, the NMDA model defined in a mod file is used, else the default NMDA model of HippoUnit
            exec("self.nmda= h."+self.NMDA_name+"(xloc, sec=self.dendrite)")
            #self.nmda.gmax = AMPA_weight/self.AMPA_NMDA_ratio
            if pathway == 'PP':
                self.nmda.gmax = AMPA_weight/self.PP_AMPA_NMDA_ratio
            elif pathway == 'SC':
                self.nmda.gmax = AMPA_weight/self.SC_AMPA_NMDA_ratio
            else:
                self.nmda.gmax = AMPA_weight/self.AMPA_NMDA_ratio
        else:
            try:
                exec("self.nmda = h."+self.default_NMDA_name+"(xloc, sec=self.dendrite)")
            except:
                h.nrn_load_dll(self.default_NMDA_path + self.libpath)
                exec("self.nmda = h."+self.default_NMDA_name+"(xloc, sec=self.dendrite)")

        self.ndend = ndend
        self.xloc = xloc
        
    def set_single_pointer(self):
        """Used in ObliqueIntegrationTest"""

        self.fakecell = h.Section(name='fakecell')


        self.epsp_ic = h.IClamp(self.fakecell(0.5))
        self.epsp_ic.amp = 1
        self.epsp_ic.dur = 0.05  # let it be shorter than the interval bw synapse activation
        self.epsp_ic.delay = self.start


        h.setpointer(self.epsp_ic._ref_i, 'pre', self.ampa)
        h.setpointer(self.epsp_ic._ref_i, 'pre', self.nmda)
        # print "pointer is set"
        
    def run_syn_pathway(self, dend_loc, weight, pathway):
        """Used in PathwayInteractionTest"""

        # self.ampa_list = [None] * number
        # self.nmda_list = [None] * number
        # self.ns_list = [None] * number
        # self.ampa_nc_list = [None] * number
        # self.nmda_nc_list = [None] * number

        ndend, xloc = dend_loc

        dend_num = ndend.split('[')[1]  # to get the number of the dendrite (eg. 80 from dendrite[80])
        dend_num = int(dend_num[:-1])
        # print dend_num

        self.initialise()

        exec("self.dendrite=h." + ndend)

        if self.cvode_active:
            h.cvode_active(1)
        else:
            h.cvode_active(0)


        self.set_ampa_nmda_pathway(dend_loc, pathway, weight)

        self.set_single_pointer()



        exec("self.sect_loc=h." + str(self.soma)+"("+str(0.5)+")")

        # initiate recording
        rec_t = h.Vector()
        rec_t.record(h._ref_t)

        rec_v = h.Vector()
        rec_v.record(self.sect_loc._ref_v)

        rec_v_dend = h.Vector()
        # rec_v_dend.record(self.shead[0](0.5)._ref_v)
        rec_v_dend.record(self.dendrite(self.xloc)._ref_v)


        h.stdinit()

        dt = 0.025
        h.dt = dt
        h.steps_per_ms = 1 / dt
        h.v_init = self.v_init #-80

        h.celsius = self.celsius
        h.init()
        h.tstop =650
        h.run()

        # get recordings
        t = numpy.array(rec_t)
        v = numpy.array(rec_v)
        v_dend = numpy.array(rec_v_dend)

        return t, v, v_dend


    def set_ampa_nmda_multiple_loc_theta(self, dend_loc, pathway, AMPA_weight, num_trains, num_stimuli_in_train):
        """Used in PathwayInteractionTest"""

        # ndend, xloc, loc_type = dend_loc

        # exec("self.dendrite=h." + ndend)
        
        h("Deadtime_GLU = 0.025")
        h("Deadtime_NMDA = 0.025")

        for j in range(num_trains):
            for k in range(num_stimuli_in_train):
                for i in range(len(dend_loc)):

                    ndend, xloc = dend_loc[i]
                    exec("self.dend=h." + ndend)


                    if self.AMPA_name: # if this is given, the AMPA model defined in a mod file is used, else the built in Exp2Syn
                        exec("self.synapse_lists[\'ampa_list_"+ pathway + "\'][j][k][i] = h."+self.AMPA_name+"("+str(xloc)+", sec=self.dend)")
                        self.synapse_lists['ampa_list_'+pathway][j][k][i].gmax=AMPA_weight
                    else: 
                        print("Give the name of the NMDA receptor of the Poirazi model or use the general ModelLoader class of HippoUnit to use the default AMPA receptor model")

                    if self.NMDA_name: # if this is given, the NMDA model defined in a mod file is used, else the default NMDA model of HippoUnit
                        # exec("self.nmda_list[i] = h."+self.NMDA_name+"(0.5, sec=self.shead[i])")
                        exec("self.synapse_lists[\'nmda_list_"+ pathway + "\'][j][k][i] = h."+self.NMDA_name+"("+str(xloc)+", sec=self.dend)")
                        #self.synapse_lists['nmda_list_'+pathway][i].gmax=AMPA_weight/self.AMPA_NMDA_ratio
                        
                        if pathway == 'PP':
                            self.synapse_lists['nmda_list_'+pathway][j][k][i].gmax=AMPA_weight/self.PP_AMPA_NMDA_ratio
                        elif pathway == 'SC':
                            self.synapse_lists['nmda_list_'+pathway][j][k][i].gmax=AMPA_weight/self.SC_AMPA_NMDA_ratio
                        else:
                            self.synapse_lists['nmda_list_'+pathway][j][k][i].gmax=AMPA_weight/self.AMPA_NMDA_ratio
                    else:
                        print("Give the name of the NMDA receptor of the Poirazi model or use the general ModelLoader class of HippoUnit to use the default NMDA receptor model")

        # self.ndend = ndend
        # self.xloc = xloc
        
    def set_pointers_multiple_loc_theta(self, dend_loc, AMPA_weight, pathway, interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train):
        """Used in PathwayInteractionTest"""
        
        self.fakecell = h.Section(name='fakecell')


        for j in range(num_trains):
            for k in range(num_stimuli_in_train):
                for i in range(len(dend_loc)):

                    self.synapse_lists['epsp_ic_list_'+pathway][j][k][i] = h.IClamp(self.fakecell(0.5))
                    self.synapse_lists['epsp_ic_list_'+pathway][j][k][i].amp = 1
                    self.synapse_lists['epsp_ic_list_'+pathway][j][k][i].dur = 0.05  # let it be shorter than the interval bw synapse activation
                    self.synapse_lists['epsp_ic_list_'+pathway][j][k][i].delay = self.start + (j * interval_bw_trains) + (k * interval_bw_stimuli_in_train)
                    
                    h.setpointer(self.synapse_lists['epsp_ic_list_'+pathway][j][k][i]._ref_i, 'pre', self.synapse_lists['ampa_list_'+pathway][j][k][i])
                    h.setpointer(self.synapse_lists['epsp_ic_list_'+pathway][j][k][i]._ref_i, 'pre', self.synapse_lists['nmda_list_'+pathway][j][k][i])
                #print('delay: ', self.synapse_lists['epsp_ic_list_'+pathway][j][k][i].delay)



    def activate_theta_stimuli(self, dend_loc, AMPA_weight, pathway, interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train):


        # self.ampa_list = [None] * len(dend_loc)
        # self.nmda_list = [None] * len(dend_loc)
        # self.ns_list = [None] * len(dend_loc)
        # self.ampa_nc_list = [None] * len(dend_loc)
        # self.nmda_nc_list = [None] * len(dend_loc)
        # self.ampa_nc_list = [[None]*len(dend_loc) for i in range(num_of_trains)]
        # self.nmda_nc_list = [[None]*len(dend_loc) for i in range(num_of_trains)]
        # self.ns_list = [[None]*len(dend_loc) for i in range(num_of_trains)]
        self.synapse_lists.update({'ampa_list_' + pathway : [[[None]*len(dend_loc) for i in range(num_stimuli_in_train)] for j in range(num_trains)],
                            'nmda_list_' + pathway : [[[None]*len(dend_loc) for i in range(num_stimuli_in_train)] for j in range(num_trains)],
                            'epsp_ic_list_' + pathway : [[[None]*len(dend_loc) for i in range(num_stimuli_in_train)] for j in range(num_trains)]
                            })  # if synapses of one of the pathways exist already, the dictionary shouldn't be overwritten, but new items are added, therefore 'update' is used. 

        # self.block_Na()

        self.set_ampa_nmda_multiple_loc_theta(dend_loc, pathway, AMPA_weight, num_trains, num_stimuli_in_train)
        self.set_pointers_multiple_loc_theta(dend_loc, AMPA_weight, pathway, interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train)
