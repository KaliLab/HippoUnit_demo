import json
import collections

""" !!! EDIT THIS PART !!! """
# Load data

with open('/mnt/extra31/Modellezo_csapat/Sara/features_for_validation/ALL_features_newBluePyEfe/FINALLY_USED/cNAC-noljp_ALL_features_new_output_nanmean_cell_False/features.json') as f:
    data = json.load(f, object_pairs_hook=collections.OrderedDict)

#num_cells = 31
#num_cells = 4
#num_cells = 3
num_cells = 5   # only those feature will be taken into account, that could and have been evaluated for all the cells in the given data set

with open('/mnt/extra31/Modellezo_csapat/Sara/features_for_validation/feature_selection.json') as f:
    feature_selection = json.load(f, object_pairs_hook=collections.OrderedDict)


valid_data = collections.OrderedDict()
#valid_data["features"] = collections.OrderedDict()

#main_key = data.keys()[0]
stim_list = []

for key, stim_dict in data.iteritems():
	split_key = key.split('_')
	#print split_key
	stim = split_key[0].capitalize() + split_key[1]
	#print stim
	stim_list.append(float(split_key[1]))
	for i in range(len(data[key]["soma"])):
		
		feat = data[key]["soma"][i]['feature']
		value = data[key]["soma"][i]['val']
		
		if data[key]["soma"][i]['n'] == num_cells: 

			if float(split_key[1]) > 0 and (feat in feature_selection["all_features"]["first_problematic"] or feat in feature_selection["all_features"]["only_positive"] or feat in feature_selection["all_features"]["negative_and_positive"]):

				valid_data[feat + "." + stim]= collections.OrderedDict({"Std" : str(value[1]), "Mean": str(value[0]), "Weight": "1", "Stimulus": stim, "Type": feat})

			elif float(split_key[1]) < 0 and (feat in feature_selection["all_features"]["only_negative"] or feat in feature_selection["all_features"]["negative_and_positive"]):
				valid_data[feat + "." + stim]= collections.OrderedDict({"Std" : str(value[1]), "Mean": str(value[0]), "Weight": "1", "Stimulus": stim, "Type": feat})

'''
sorted_stim_list = sorted(stim_list)
stim_data = {"stimuli" : collections.OrderedDict()}

for i in range(len(sorted_stim_list)):
	stim_data["stimuli"].update({"Step"+str(sorted_stim_list[i]) : 
	    {"Delay": "500",
            "RecLocationX": "0.5",
            "TTX": "0.0",
            "StimLocationX": "0.5",
            "Threshold": "-20.0",
            "Amplitude": str(sorted_stim_list[i]),
            "Duration": "800.0",
            "HypAmp": "0.0",
            "Type": "SquarePulse",
            "StimSectionName": "soma[0]",
            "RecSectionName": "soma[0]"}})
'''
	

ordered_valid_data = {}
ordered_valid_data = collections.OrderedDict(sorted(valid_data.items()))
#print valid_data
print len(ordered_valid_data.keys())

"""Check if there are differences"""
'''
for key in sample_data.keys():
	if key not in ordered_valid_data["features"].keys():
		print key
for key in ordered_valid_data["features"].keys():
	if key not in sample_data.keys():
		print key
'''
		
file_name_valid_data= '/mnt/extra31/Modellezo_csapat/Sara/features_for_validation/ALL_features_newBluePyEfe/FINALLY_USED/validation_json_files/feat_CA1_int_cNAC_more_features.json'
json.dump(ordered_valid_data, open(file_name_valid_data, "wb"), indent=4)

