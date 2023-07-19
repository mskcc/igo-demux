import pandas as pd
import sys
import os
import json
import os.path
import requests
import shutil
import glob


ENDPOINT = "https://igolims.mskcc.org:8443/LimsRest/getConfig?igoId="
original_tiff_images_directory = "/skimcs/mohibullahlab/IGO_Pipeline_Results/Single_Cell/10X_Genomics/TIFF_Images/"
tiff_images_directory = "/igo/work/igo/TIFF_Images/"


# sample_id can be get from sample sheet, will be the part in front of _IGO_
class Spatial_sample:
    def __init__(self, sample, project_id):
        self.sample_name = sample.split("_IGO_")[0]
        self.IGO_ID = sample.split("_IGO_")[1]
        self.chip_position = "EMPTY"
        self.chip_id = "EMPTY"
        self.preservation = "EMPTY"
        self.tiff_image = "EMPTY"
        self.get_info_from_LIMS()
        self.copy_tiff(project_id)

    def get_info_from_LIMS(self):
        response = requests.get(ENDPOINT + self.IGO_ID , auth = ("pms", "tiagostarbuckslightbike"), verify = False)
        response_data = json.loads(response.text.encode("utf8"))
        self.chip_position = response_data["chipPosition"]
        self.chip_id = response_data["chipID"]
        self.preservation = response_data["preservation"]
    
    def copy_tiff(self, project_id):
        # project_id format as Project_12345
        source_loc_dir = original_tiff_images_directory + project_id
        destination_loc = tiff_images_directory + project_id + "/" + self.sample_name + ".tif"
        # create TIFF_images director if not exists
        if not os.path.exists(source_loc_dir):
            os.makedirs(source_loc_dir)

        # copy all the image files using rsync?
        original_tiff_image = glob.glob(source_loc_dir + "/" + self.sample_name + "*")
        if len(original_tiff_image) != 0 or "tif" not in original_tiff_image[0]:
            print("tif file is not in proper format for sample {}, please check".format(self.IGO_ID))
        else:
            shutil.copy(original_tiff_image[0], destination_loc)
            self.tiff_image = destination_loc

