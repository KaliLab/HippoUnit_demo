# Added by Sara Saray and Shailesh Appukuttan to run validation tests and register the results in the HBP Validation Framework

from __future__ import print_function
import os
from hippounit.utils import ModelLoader

class Golding_et_al_2001_weak_fig8A(ModelLoader):
    # model UUID on HBP Model Catalog
    model_uuid = "ca3008af-d251-4691-a95a-ac071ca10b43"
    model_version = "fig8A"

    def __init__(self):
        # path to mod files
        model_path = os.path.dirname(os.path.abspath(__file__))
        print("model_path = ", model_path)
        mod_files_path = model_path + "/fig08/"

        #Load cell model
        ModelLoader.__init__(self, mod_files_path = mod_files_path)

        # outputs will be saved in subfolders named like this:
        self.name="Golding_2001_weak_prop_fig8A"

        # path to hoc file
        # the model must not display any GUI!!
        self.hocpath = model_path + "/fig08/main_model_weak_prop_fig8A.hoc"

        # If the hoc file doesn't contain a template, this must be None (the default value is None)
        self.template_name = None

        # model.SomaSecList_name should be None, if there is no Section List in the model for the soma, or if the name of the soma section is given by setting model.soma (the default value is None)
        self.SomaSecList_name = None
        # if the soma is not in a section list or to use a specific somatic section, add its name here:
        self.soma = 'somaA'

        # For the PSP Attenuation Test, and Back-propagating AP Test a section list containing the trunk sections is needed
        self.TrunkSecList_name = 'primary_apical_list'
        # For the Oblique Integration Test a section list containing the oblique dendritic sections is needed
        self.ObliqueSecList_name = 'oblique_dendrites'

        self.v_init = -70
        self.celsius = 35

class Golding_et_al_2001_strong_fig8B(ModelLoader):
    # model UUID on HBP Model Catalog
    model_uuid = "ca3008af-d251-4691-a95a-ac071ca10b43"
    model_version = "fig8B"

    def __init__(self):
        # path to mod files
        model_path = os.path.dirname(os.path.abspath(__file__))
        print("model_path = ", model_path)
        mod_files_path = model_path + "/fig08/"

        #Load cell model
        ModelLoader.__init__(self, mod_files_path = mod_files_path)

        # outputs will be saved in subfolders named like this:
        self.name="Golding_2001_strong_prop_fig8B"

        # path to hoc file
        # the model must not display any GUI!!
        self.hocpath = model_path + "/fig08/main_model_strong_prop_fig8B.hoc"

        # If the hoc file doesn't contain a template, this must be None (the default value is None)
        self.template_name = None

        # model.SomaSecList_name should be None, if there is no Section List in the model for the soma, or if the name of the soma section is given by setting model.soma (the default value is None)
        self.SomaSecList_name = None
        # if the soma is not in a section list or to use a specific somatic section, add its name here:
        self.soma = 'somaA'

        # For the PSP Attenuation Test, and Back-propagating AP Test a section list containing the trunk sections is needed
        self.TrunkSecList_name = 'primary_apical_list'
        # For the Oblique Integration Test a section list containing the oblique dendritic sections is needed
        self.ObliqueSecList_name = 'oblique_dendrites'

        self.v_init = -70
        self.celsius = 35

class Golding_et_al_2001_weak_fig9B(ModelLoader):
    # model UUID on HBP Model Catalog
    model_uuid = "ca3008af-d251-4691-a95a-ac071ca10b43"
    model_version = "fig9B"

    def __init__(self):
        # path to mod files
        model_path = os.path.dirname(os.path.abspath(__file__))
        print("model_path = ", model_path)
        mod_files_path = model_path + "/fig09/"

        #Load cell model
        ModelLoader.__init__(self, mod_files_path = mod_files_path)

        # outputs will be saved in subfolders named like this:
        self.name="Golding_2001_weak_prop_fig9B"

        # path to hoc file
        # the model must not display any GUI!!
        self.hocpath = model_path + "/fig09/main_model_weak_prop_fig9B.hoc"

        # If the hoc file doesn't contain a template, this must be None (the default value is None)
        self.template_name = None

        # model.SomaSecList_name should be None, if there is no Section List in the model for the soma, or if the name of the soma section is given by setting model.soma (the default value is None)
        self.SomaSecList_name = None
        # if the soma is not in a section list or to use a specific somatic section, add its name here:
        self.soma = 'somaA'

        # For the PSP Attenuation Test, and Back-propagating AP Test a section list containing the trunk sections is needed
        self.TrunkSecList_name = 'primary_apical_list'
        # For the Oblique Integration Test a section list containing the oblique dendritic sections is needed
        self.ObliqueSecList_name = 'oblique_dendrites'

        self.v_init = -70
        self.celsius = 35

