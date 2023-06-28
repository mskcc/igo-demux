import tempfile
import os
import organise_fastq_split_by_lane

def test_create_fastq_folders_normal():
    temp_dir = tempfile.TemporaryDirectory()

    # Create the folder inside the temporary directory
    folder_name = "Project_07323_R"
    project_folder_path = os.path.join(temp_dir.name, folder_name)
    os.mkdir(project_folder_path)

    # Create a fastq inside the folder
    file_name1 = "IR116-RNA5_IGO_07323_R_5_S5_R1_001.fastq.gz"
    file_path1 = os.path.join(project_folder_path, file_name1)
    with open(file_path1, 'w') as file1:
        file1.write("Example fastq")

       # Create a fastq inside the folder
    file_name2 = "IR122-RNA6_IGO_07323_R_6_S6_R1_001.fastq.gz"
    file_path2 = os.path.join(project_folder_path, file_name2)
    with open(file_path2, 'w') as file2:
        file2.write("Example fastq")

    # Print the paths for verification
    print("Temporary Directory:", temp_dir.name)
    print("Folder Path:", project_folder_path)
    print("File Path:", file_path1)

    organise_fastq_split_by_lane.create_fastq_folders(temp_dir.name)
    final_path = os.path.join(project_folder_path, "Sample_IR116-RNA5_IGO_07323_R_5", file_name1)
    print("FINAL PATH:", final_path)
    assert(os.path.exists(final_path))

# TODO add test case when sample name is the same but IGO ID is different


    