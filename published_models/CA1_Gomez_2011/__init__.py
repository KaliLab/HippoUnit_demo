# Added by Sara Saray and Shailesh Appukuttan to run validation tests and register the results in the HBP Validation Framework

from __future__ import print_function
import os
from hippounit.utils import ModelLoader

class CA1_Gomez_2011_n123(ModelLoader):
    # model UUID on HBP Model Catalog
    model_uuid = "42cb7f65-22e6-44b9-8a00-d70b7aa13b75"
    model_version = "n123_morph"

    def __init__(self):
        # path to mod files
        model_path = os.path.dirname(os.path.abspath(__file__))
        print("model_path = ", model_path)
        mod_files_path = model_path + "/mechanism/"

        #Load cell model
        ModelLoader.__init__(self, mod_files_path = mod_files_path)

        # outputs will be saved in subfolders named like this:
        self.name="Gomez_Gonzalez_2011_n123_morph"

        # path to hoc file
        # the model must not display any GUI!!
        self.hocpath = model_path + "/experiment/junio/main_model_n123_morph.hoc"

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

        self.v_init = -66
        self.celsius = 34


        #setting synapse parameters
        self.NMDA_name = 'NMDA'
        self.AMPA_name = 'GLU'
        self.AMPA_NMDA_ratio = 1/0.396


class CA1_Gomez_2011_n125(ModelLoader):
    # model UUID on HBP Model Catalog
    model_uuid = "42cb7f65-22e6-44b9-8a00-d70b7aa13b75"
    model_version = "n125_morph"

    def __init__(self):
        # path to mod files
        model_path = os.path.dirname(os.path.abspath(__file__))
        print("model_path = ", model_path)
        mod_files_path = model_path + "/mechanism/"

        #Load cell model
        ModelLoader.__init__(self, mod_files_path = mod_files_path)

        # outputs will be saved in subfolders named like this:
        self.name="Gomez_Gonzalez_2011_n125_morph"

        # path to hoc file
        # the model must not display any GUI!!
        self.hocpath = model_path + "/experiment/junio/main_model_n125_morph.hoc"

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

        self.v_init = -66
        self.celsius = 34


        #setting synapse parameters
        self.NMDA_name = 'NMDA'
        self.AMPA_name = 'GLU'
        self.AMPA_NMDA_ratio = 1/0.396


class CA1_Gomez_2011_n128(ModelLoader):
    # model UUID on HBP Model Catalog
    model_uuid = "42cb7f65-22e6-44b9-8a00-d70b7aa13b75"
    model_version = "n128_morph"

    def __init__(self):
        # path to mod files
        model_path = os.path.dirname(os.path.abspath(__file__))
        print("model_path = ", model_path)
        mod_files_path = model_path + "/mechanism/"

        #Load cell model
        ModelLoader.__init__(self, mod_files_path = mod_files_path)

        # outputs will be saved in subfolders named like this:
        self.name="Gomez_Gonzalez_2011_n128_morph"

        # path to hoc file
        # the model must not display any GUI!!
        self.hocpath = model_path + "/experiment/junio/main_model_n128_morph.hoc"

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

        self.v_init = -66
        self.celsius = 34


        #setting synapse parameters
        self.NMDA_name = 'NMDA'
        self.AMPA_name = 'GLU'
        self.AMPA_NMDA_ratio = 1/0.396


class CA1_Gomez_2011_n129(ModelLoader):
    # model UUID on HBP Model Catalog
    model_uuid = "42cb7f65-22e6-44b9-8a00-d70b7aa13b75"
    model_version = "n129_morph"

    def __init__(self):
        # path to mod files
        model_path = os.path.dirname(os.path.abspath(__file__))
        print("model_path = ", model_path)
        mod_files_path = model_path + "/mechanism/"

        #Load cell model
        ModelLoader.__init__(self, mod_files_path = mod_files_path)

        # outputs will be saved in subfolders named like this:
        self.name="Gomez_Gonzalez_2011_n129_morph"

        # path to hoc file
        # the model must not display any GUI!!
        self.hocpath = model_path + "/experiment/junio/main_model_n129_morph.hoc"

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

        self.v_init = -66
        self.celsius = 34


        #setting synapse parameters
        self.NMDA_name = 'NMDA'
        self.AMPA_name = 'GLU'
        self.AMPA_NMDA_ratio = 1/0.396


class CA1_Gomez_2011_n130(ModelLoader):
    # model UUID on HBP Model Catalog
    model_uuid = "42cb7f65-22e6-44b9-8a00-d70b7aa13b75"
    model_version = "n130_morph"

    def __init__(self):
        # path to mod files
        model_path = os.path.dirname(os.path.abspath(__file__))
        print("model_path = ", model_path)
        mod_files_path = model_path + "/mechanism/"

        #Load cell model
        ModelLoader.__init__(self, mod_files_path = mod_files_path)

        # outputs will be saved in subfolders named like this:
        self.name="Gomez_Gonzalez_2011_n130_morph"

        # path to hoc file
        # the model must not display any GUI!!
        self.hocpath = model_path + "/experiment/junio/main_model_n130_morph.hoc"

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

        self.v_init = -66
        self.celsius = 34


        #setting synapse parameters
        self.NMDA_name = 'NMDA'
        self.AMPA_name = 'GLU'
        self.AMPA_NMDA_ratio = 1/0.396
