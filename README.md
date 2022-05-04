# HippoUnit demo

The purpose of this repository is to provide examples on how to run the validation tests of [HippoUnit](https://github.com/KaliLab/hippounit) on models of the hippocampal CA1 pyramidal cell available in the literature and on [ModelDB](https://senselab.med.yale.edu/modeldb/)

## Journal Paper

Sára Sáray, Christian A. Rössert, Shailesh Appukuttan, Rosanna Migliore, Paola Vitale, Carmen A. Lupascu, Luca L. Bologna, Werner Van Geit, Armando Romani, Andrew P. Davison, Eilif Muller, Tamás F. Freund, Szabolcs Káli: Systematic comparison and automated validation of detailed models of hippocampal neurons, PLOS Computational Biology; 
https://doi.org/10.1371/journal.pcbi.1008114  

## Models tested

Models from literature published on ModelDB typically implement their own simulations and plots to make it possible for users and readers to reproduce and visualize the results that have been shown in the corresponding paper. Therefore to be able to test the models described below using HippoUnit, we needed to create standalone versions of them. These standalone versions should not display any GUI, or contain any built in simulations and run-time modifications, but their behavior should be identical to the published version of the models. We also added section lists of the radial oblique and the trunk dendritic sections to those models where this was not done yet, as some of the tests require the implementation of these. 

Tested models (containing the modified hoc files to be tested) are available in the *published_models* folder of this repository (https://github.com/KaliLab/HippoUnit_demo/tree/master/published_models)


The **Golding et al. (2001) model**  was developed to show the dichotomy of the back-propagation efficacy and the amplitudes of the back-propagating action potentials at distal trunk regions in CA1 pyramidal cells and to make predictions on the possible causes of this behaviour. It contains only the most important ion channels (Na, KDR,KA) needed to reproduce the generation and propagation of action potentials. Here we tested three different versions of the model: the ones belonging to the Figure 8 A, Figure 8 B and Figure 9 B of the corresponding paper.

(Golding NL, Kath WL, Spruston N (2001) Dichotomy of action-potential backpropagation in CA1 pyramidal neuron dendrites. J Neurophysiol 86:2998-3010; ModelDB accession number: 64167)

Modified hoc files used for running the validation tests of HippoUnit on the model: 
* HippoUnit_demo/published_models/Golding_et_al_2001_dichotomy/fig08/main_model_strong_prop_fig8B.hoc
* HippoUnit_demo/published_models/Golding_et_al_2001_dichotomy/fig08/main_model_weak_prop_fig8A.hoc
* HippoUnit_demo/published_models/Golding_et_al_2001_dichotomy/fig09/main_model_weak_prop_fig9B.hoc


The **Katz et al. (2009) model** is based on the Golding et al. (2001) model and was built to investigate the functional consequences of the distribution of strength and density of synapses on the apical dendrites, that they observed experimentally, on the mode of dendritic integration.

(Katz Y, Menon V, Nicholson DA, Geinisman Y, Kath WL, Spruston N (2009) Synapse distribution suggests a two-stage model of dendritic integration in CA1 pyramidal neurons. Neuron 63:171-7; ModelDB accession number: 127351)

Modified hoc file used for running the validation tests of HippoUnit on the model: 
* HippoUnit_demo/published_models/Katz_et_al_2009_2stageintegration_code/main_model.hoc


The **Migliore et al. (2011) model** was used to study schizophrenic behaviour, and is based on models of the same modeling group, which were used to investigate the initiation and propagation of action potentials in oblique dendrites, and have been validated against a number of different electrophysiological data.

(Migliore M, De Blasi I, Tegolo D, Migliore R (2011) A modeling study suggesting how a reduction in the context-dependent input on CA1 pyramidal neurons could generate schizophrenic behavior. Neural Netw 24:552-9; ModelDB accession number: 138205)

Modified hoc file used for running the validation tests of HippoUnit on the model: 
* HippoUnit_demo/published_models/Migliore_et_al_2011_Schizophr/main_model.hoc


The **Bianchi et al (2012) model** was designed to show the mechanisms behind depolarization block observed experimentally in the somatic spiking behavior of CA1 pyramidal cells. It has been developed by combining and modifying the Shah et al. (2008) and the Poirazi et al. (2003) models. The former was developed to show the significance of axonal M-type potassium channels.

(Bianchi D, Marasco A, Limongiello A, Marchetti C, Marie H, Tirozzi B, Migliore M (2012) On the mechanisms underlying the depolarization block in the spiking dynamics of CA1 pyramidal neurons. J Comput Neurosci 33:207-25; ModelDB accession number: 143719)

Modified hoc file used for running the validation tests of HippoUnit on the model: 
* HippoUnit_demo/published_models/Ca1_Bianchi_2012/experiment/main_model.hoc


The **Poirazi et al. (2003) model** was designed to clarify the issues about the integrative properties of thin apical dendrites, that may arise from the different and sometimes conflicting interpretation of available experimental data. This is a quite complex model in the sense that it contains a larger number of different types of ion channels, whose properties were adjusted to fit in vitro experimental data, and it also contains four types of synaptic receptors.

(Poirazi P, Brannon T, Mel BW (2003) Arithmetic of subthreshold synaptic summation in a model CA1 pyramidal cell. Neuron 37:977-87; Poirazi P, Brannon T, Mel BW (2003) Pyramidal neuron as two-layer neural network. Neuron 37:989-99; ModelDB accession number: 20212)

Modified hoc file used for running the validation tests of HippoUnit on the model:
* HippoUnit_demo/published_models/Poirazi_et_al_2003_CA1_multi/experiment/cluster-dispersion/main_model.hoc


The **Gómez González et al. (2011) model** is based on the Poirazi et al. (2003) model and it was modified to replicate the experimental data of Losonczy and Magee (2006) (A. Losonczy and J. C. Magee (2006) Neuron, 50: 291-307) on the nonlinear signal integration of radial oblique dendrites when the inputs arrive in a short time window. The model was adjusted to five different detailed morphologies.

(Gómez González JF, Mel BW, Poirazi P (2011) Distinguishing Linear vs. Non-Linear Integration in CA1 Radial Oblique Dendrites: It's about Time. Front Comput Neurosci 5:44; ModelDB accession number: 144450)

Modified hoc files used for running the validation tests of HippoUnit on the model:
* HippoUnit_demo/published_models/CA1_Gomez_2011/experiment/junio/main_model_n123_morph.hoc
* HippoUnit_demo/published_models/CA1_Gomez_2011/experiment/junio/main_model_n125_morph.hoc
* HippoUnit_demo/published_models/CA1_Gomez_2011/experiment/junio/main_model_n128_morph.hoc
* HippoUnit_demo/published_models/CA1_Gomez_2011/experiment/junio/main_model_n129_morph.hoc
* HippoUnit_demo/published_models/CA1_Gomez_2011/experiment/junio/main_model_n130_morph.hoc

## Target features

The *target_features* folder of this repository contains the target experimental data to each tests of HippoUnit in JSON files to which the models' behaviour are compared quantitatively. The target data always consists of the mean and the standard deviation values of the features of interest, that are derived from the statistical analysis of the results of experimental measurements on several different cells.

* **feat_CA1_pyr_cACpyr_more_features.json**: contains target feature values for the **Somatic Features Test** of HippoUnit obtained from sharp electrode recordings in adult rat CA1 pyramidal cells made at UCL. (Migliore et al. (2018) PLOS Computational Biology, 14:9; https://doi.org/10.1371/journal.pcbi.1006423)

* *feat_CA1_int_bAC_more_features.json, feat_CA1_int_cAC_more_features.json, feat_CA1_int_cNAC_more_features.json*: contain target feature values for the **Somatic Features Test** of HippoUnit obtained from sharp electrode recordings in adult rat CA1 interneurons made at UCL. (Migliore et al., PLOS Computational Biology, vol. 14, no. 9, 2018; https://doi.org/10.1371/journal.pcbi.1006423) These data are not used for the models tested here.

* **feat_rat_CA1_JMakara_more_features.json**: contains target feature values for the **Somatic Features Test** of HippoUnit obtained from patch clamp recordings in rat CA1 pyramidal cells made by Judit Makara.

* **feat_PSP_attenuation_target_data.json**: contains target feature values for the **PSP Attenuation Test** of HippoUnit obtained from Figure 1 E and Figure 2 B of the paper Magee & Cook 2000 (Nat Neurosci 3:895-903; DOI: 10.1038/78800) using the [DigitizeIt](https://www.digitizeit.de/) software.

* **feat_backpropagating_AP_target_data.json**: contains target feature values for the **Backpropagating AP Test** of HippoUnit obtained from Figure 1 B of the paper Golding et al. 2001 (J Neurophysiol 86:2998-3010; DOI: 10.1152/jn.2001.86.6.2998) using the [DigitizeIt](https://www.digitizeit.de/) software.

* **depol_block_target_data.json**: contains target feature values for the **Depolarization Block Test** of HippoUnit obtained from the paper Bianchi et al. 2012 (J Comput Neurosci, 33: 207-225; doi: 10.1007/s10827-012-0383-y)

* **oblique_target_data.json**: contains target feature values for the **Oblique Integration Test** of HippoUnit obtained from the paper A. Losonczy and J. C. Magee 2006 (Neuron, 50: 291-307; DOI: 10.1016/j.neuron.2006.03.016)

* **pathway_interaction_target_data.json**: contains target feature values for the **Pathway Interaction Test** of HippoUnit obtained from the paper Takahashi and Magee 2009 (Neuron, 62: 102–111; DOI: 10.1016/j.neuron.2009.03.007)

The **Examples_on_creating_JSON_files** subfolder contains some examples on how the JSON files were created.

*  The **Somatic_Features** subfolder contains an example script on how to convert the output (features.json) of [BluePyEfe](https://github.com/BlueBrain/BluePyEfe) to a stimulus and a feature JSON file that are compatible with the Somatic Features Test of HippoUnit.

*  The **Magee2000-PSP_att** subfolder contains data digitized from the figures of the paper and a script (extract.py) that converts these data for the feature and stimulus JSON files of the PSP Attenuation Test of HippoUnit.


## Jupyter notebooks

The *jupyter_notebooks* folder of this repository contains Jupyter notebooks to each of the tested models to show how to run the validation tests of HippoUnit on these models.


## Validation results

The *published_models_validation_results* folder of this repository contains all the outputs of the tests of HippoUnit that were obtained by running the Jupyter notebooks available in the *jupyter_notebooks* folder to test the models in the *published_models* folder against the target data available in the *target_features* folder of this repository.



## Acknowledgments

This research has been partially supported by the Human Brain Project, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1 and SGA2).

SUPPORTED BY THE ÚNKP-19-3-III NEW NATIONAL EXCELLENCE PROGRAM OF THE MINISTRY FOR INNOVATION AND TECHNOLOGY.  

This research has been partially supported by the European Union, co-financed by the European Social Fund (EFOP-3.6.3-VEKOP- 16-2017-00002 ).

