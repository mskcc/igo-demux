import scripts.deliver_pipeline as deliver_pipeline
import os
import tempfile


def test_find_bams():
    tmpdir = tempfile.mkdtemp()
    run_path = os.path.join(tmpdir, "DIANA_0479_BHM2NVDSX3")
    os.mkdir(run_path)
    bam_file = os.path.join(run_path, "DIANA_0479_BHM2NVDSX3___P12785_H___GA28_ot_IGO_12785_H_1.bam")
    bam = open(bam_file, "w")
    bam.write("SOME BAM INFORMATION FOR YOU")
    bam.close()

    print("Created {}".format(bam_file))
    
    bamdict = deliver_pipeline.find_bams("12785_H", tmpdir)
    print(bamdict)
    assert(len(bamdict["12785_H_1"]) == 1)

    deliver_pipeline.write_bams_to_share(bamdict, tmpdir +"/abcdefg")

    assert(os.path.isfile(tmpdir + "/abcdefg/GA28_ot_IGO_12785_H_1.bam"))


def test_merge_bams_command():
    tmpdir = tempfile.mkdtemp()
    bamdict = {"12785_H_1": ["DIANA_0479_BHM2NVDSX3___P12785_H___GA28_ot_IGO_12785_H_1.bam", "DIANA_0480_BHM2NVDSX3___P12785_H___GA28_ot_IGO_12785_H_1.bam"]}
    bsub_commands = deliver_pipeline.write_bams_to_share(bamdict,tmpdir)
    assert(bsub_commands[0].count("--INPUT") == 2)
