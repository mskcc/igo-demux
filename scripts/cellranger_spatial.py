import pandas as pd
import sys
import os
import json
import os.path
import requests
import shutil


ENDPOINT = "https://igolims.mskcc.org:8443/LimsRest/getConfig?igoId="
original_tiff_images_directory = "/skimcs/mohibullahlab/IGO_Pipeline_Results/Single_Cell/10X_Genomics/TIFF_Images/"
tiff_images_directory = "/igo/work/igo/TIFF_Images/"


# sample_id can be get from sample sheet, will be the part in front of _IGO_
class Spatial_sample:
    def __init__(self, sample):
        self.sample_name = sample.split("_IGO_")[0]
        self.IGO_ID = sample.split("_IGO_")[1]
        self.chip_position = "EMPTY"
        self.chip_id = "EMPTY"
        self.preservation = "EMPTY"
        # TODO need to be confirmed
        self.tiff_image = tiff_images_directory
        self.get_info_from_LIMS()

    def get_info_from_LIMS(self):
        response = requests.get(ENDPOINT + self.IGO_ID , auth = ("pms", "tiagostarbuckslightbike"), verify = False)
        response_data = json.loads(response.text.encode("utf8"))
        self.chip_position = response_data["chipPosition"]
        self.chip_id = response_data["chipID"]
        self.preservation = response_data["preservation"]
    
    # rule for tiff file name?
    def get_tiff_location(self):
        sample_tiff_image = [tiff_image for tiff_image in tiff_images if sample_tiff_id in tiff_image][0]
        sample_tiff_directory_and_image = "{}/{}".format(tiff_images_directory, sample_tiff_image)

def copy_tiff(project_ID):
    source_loc = original_tiff_images_directory + project_ID
    destination_loc = tiff_images_directory + project_ID
	# create TIFF_images director if not exists
    if not os.path.exists(destination_loc):
        os.makedirs(destination_loc)

    # copy all the image files using rsync?
	tiff_images = os.listdir(source_loc)
	for tiff_image in tiff_images:
			original_tiff_image = "{}/{}".format(source_loc, tiff_image)
			print(original_tiff_image)
			shutil.copy(original_tiff_image, destination_loc)

