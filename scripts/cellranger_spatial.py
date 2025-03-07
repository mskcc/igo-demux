import os
import json
import os.path
import requests
import shutil
import scripts.cellranger_config as CONFIG

# sample_id can be get from sample sheet, will be the part in front of _IGO_
class Spatial_sample:
    def __init__(self, sample, project_id):
        self.sample_name = sample.split("_IGO_")[0]
        self.IGO_ID = sample.split("_IGO_")[1]
        self.chip_position = "EMPTY"
        self.chip_id = "EMPTY"
        self.preservation = "EMPTY"
        self.tiff_image = "EMPTY"
        self.json = "EMPTY"
        self.HE_tiff_image = "EMPTY"
        self.get_info_from_LIMS()
        self.copy_tiff(project_id)
        self.copy_json(project_id)

    def get_info_from_LIMS(self):
        response = requests.get(CONFIG.VISIUM_ENDPOINT + self.IGO_ID , auth = ("pms", "tiagostarbuckslightbike"), verify = False)
        response_data = json.loads(response.text.encode("utf8"))
        self.chip_position = response_data["chipPosition"]
        self.chip_id = response_data["chipID"]
        self.preservation = response_data["preservation"]
        self.cytAssist = response_data["cytAssist"]
    
    def copy_tiff(self, project_id):
        # project_id format as Project_12345
        source_loc_dir = CONFIG.original_tiff_images_directory + project_id
        destination_loc = CONFIG.tiff_images_directory + project_id + "/CytAssist"
        destination_file = destination_loc + "/" + self.sample_name + ".tif"
        destination_HE_loc = CONFIG.tiff_images_directory + project_id + "/Microscope"
        destination_HE_file = destination_HE_loc + "/HE_" + self.sample_name + ".tif"       
        # create TIFF_images director if not exists
        if not os.path.exists(destination_loc):
            os.makedirs(destination_loc)
        # create microscope image director if not exists
        if not os.path.exists(destination_HE_loc):
            os.makedirs(destination_HE_loc)

        # copy image file per sample
        original_tiff_image = source_loc_dir + "/" + self.sample_name + ".tif"
        if os.path.isfile(original_tiff_image):
            shutil.copy(original_tiff_image, destination_file)
            self.tiff_image = destination_file
            print("copy {} to {}".format(original_tiff_image, destination_file))
        else:
            print("tif file is not in proper format for sample {}, please check".format(self.IGO_ID))

        # copy HE file per sample if exists
        original_HE_tiff_image = source_loc_dir + "/Microscope/HE_" + self.sample_name + ".tif"
        if os.path.isfile(original_HE_tiff_image):
            shutil.copy(original_HE_tiff_image, destination_HE_file)
            self.HE_tiff_image = destination_HE_file
            print("copy {} to {}".format(original_HE_tiff_image, destination_HE_file))
        else:
            print("HE tif file does not exist for sample {}, please check".format(self.IGO_ID))
    
    # copy json file if exists
    def copy_json(self, project_id):
        # project_id format as Project_12345
        source_loc = CONFIG.original_tiff_images_directory + project_id + "/json/" + self.sample_name + ".json"
        destination_loc = CONFIG.tiff_images_directory + project_id + "/json"
        destination_file = destination_loc + "/" + self.sample_name + ".json"

        # create director if not exists
        if not os.path.exists(destination_loc):
            os.makedirs(destination_loc)
        
        if os.path.isfile(source_loc):
            shutil.copy(source_loc, destination_file)
            self.json = destination_file
            print("copy {} to {}".format(source_loc, destination_file))
        else:
            print("json file does not exist for {}".format(self.sample_name))            
