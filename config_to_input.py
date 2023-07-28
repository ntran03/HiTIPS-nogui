import pandas as pd
import os
import sys


class Converter():
    def __init__(self, input_path, experiment_name):
        #convert to df and then to a dict
        self.experiment_name = experiment_name
        self.input_path = input_path
        self.config_file = pd.read_csv(os.path.join(input_path, "analysis_configuration.csv"))
        config_names = self.config_file[self.config_file.columns[0]].tolist()
        config_values = self.config_file[self.config_file.columns[1]].tolist()
        self.config_dict = {config_names[i]: config_values[i] for i in range(len(config_names))}
    def convert_to_input(self):
        #search folder for input_template
        input_template = pd.read_csv('./input_template.csv')
        #print(self.config_dict)
        for i in range(len(input_template['Name'].tolist())):
            name = input_template['Name'].iloc[i]
            #print(name)
            if name in self.config_dict:
                #print(self.config_dict[name])
                input_template.at[i,'Value'] = self.config_dict[name]
        
        #write to csv file
        print(os.path.join(self.input_path, self.experiment_name))
        input_template.to_csv(os.path.join(self.input_path, self.experiment_name))
        return input_template
                

converter = Converter(sys.argv[1],sys.argv[2])
converter.convert_to_input()
print("done")