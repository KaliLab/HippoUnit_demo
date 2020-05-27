# Added by Sara Saray and Shailesh Appukuttan to run validation tests and register the results in the HBP Validation Framework

from __future__ import print_function
import os
from hippounit.utils import ModelLoader

class Migliore_et_al_2011(ModelLoader):
    # model UUID on HBP Model Catalog
    model_uuid = "263d3d47-9241-46e6-a1fd-4ada01757f49"
    model_version = "1.0"

    def __init__(self):
        # path to mod files
        model_path = os.path.dirname(os.path.abspath(__file__))
        print("model_path = ", model_path)
        mod_files_path = model_path + "/"

        #Load cell model
        ModelLoader.__init__(self, mod_files_path = mod_files_path)

        # outputs will be saved in subfolders named like this:
        self.name="Migliore_et_al_2011"

        # path to hoc file
        # the model must not display any GUI!!
        self.hocpath = model_path + "/main_model.hoc"

        # If the hoc file doesn't contain a template, this must be None (the default value is None)
        self.template_name = None

        # model.SomaSecList_name should be None, if there is no Section List in the model for the soma, or if the name of the soma section is given by setting model.soma (the default value is None)
        self.SomaSecList_name = None
        # if the soma is not in a section list or to use a specific somatic section, add its name here:
        self.soma = 'soma[0]'

        # For the PSP Attenuation Test, and Back-propagating AP Test a section list containing the trunk sections is needed
        self.TrunkSecList_name = 'trunk_sec_list'
        # For the Oblique Integration Test a section list containing the oblique dendritic sections is needed
        self.ObliqueSecList_name = 'oblique_dendrites'

        #This will be argument to those tests, where dendritic locatins are selected according to distances. If not set, the end of the above given soma section will be used as reference point for distance determination
        self.trunk_origin = ['soma[1]', 1] 

        self.v_init = -65
        self.celsius = 35

        #setting synapse parameters
        self.AMPA_tau1 = 0.4
        self.AMPA_tau2 = 1
