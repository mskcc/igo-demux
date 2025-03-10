# function to parse bionano json file and create a file/push to lims

import json
import sys

class bionano_report:
    def __init__(self, filepath):
        json_data = self.read_json(filepath)
        self.job_id = int(json_data["job"]["value"]["jobpk"]["value"].replace(",", ""))
        self.project_id = json_data["job"]["value"]["name"]["value"]
        self.igo_id = json_data["job"]["value"]["samplename"]["value"]
        self.instrument = json_data["job"]["value"]["instrument"]["value"]
        self.chip_serial_number = json_data["job"]["value"]["serialnumber"]["value"]
        self.reference = json_data["job"]["value"]["reference"]["value"]
        self.max_mol_len = json_data["mqr"]["value"]["max_mollen"]["value"].split(" ")[0].replace(",", "")
        self.total_DNA= json_data["mqr"]["value"]["quantity"]["value"].split(" ")[0].replace(",", "")
        self.mol_n_50= json_data["mqr"]["value"]["mol_n_50"]["value"].split(" ")[0].replace(",", "")
        self.total_DNA_150kb_site9= json_data["mqr"]["value"]["quantity_ge_150kb_site9"]["value"].split(" ")[0].replace(",", "")
        self.mol_n_50_ge_150kb_site9= json_data["mqr"]["value"]["mol_n_50_ge_150kb_site9"]["value"].split(" ")[0].replace(",", "")
        self.avg_label_density= json_data["mqr"]["value"]["avg_label_density"]["value"].split(" ")[0].replace(",", "")
        self.effective_coverage= json_data["mqr"]["value"]["coverage"]["value"]
        self.map_rate= json_data["mqr"]["value"]["map_rate"]["value"].replace("%", "")
        self.base_pairs_per_pixel= json_data["mqr"]["value"]["bpp"]["value"]
        self.plv= json_data["mqr"]["value"]["fp_rate"]["value"]
        self.nlv= json_data["mqr"]["value"]["fn_rate"]["value"]
        self.mol_integrity_num= json_data["mqr"]["value"]["integrity_num"]["value"]

    def read_json(self, filename):
        # Open the JSON file
        with open(filename, 'r') as f:
            # Load the JSON data into a Python dictionary
            data = json.load(f)[0]
        return data

if __name__ == '__main__':
    filename = sys.argv[1]
    bionano = bionano_report(filename)
    for property_name, value in vars(bionano).items():
        print(f"{property_name}: {value}")
        