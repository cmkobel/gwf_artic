from gwf import *



gwf = Workflow(defaults={
    "mail_user": "kobel@pm.me",
    "mail_type": "FAIL",
    "account": "clinicalmicrobio",
})

environment_base = """source /home/cmkobel/miniconda3/etc/profile.d/conda.sh"""


titles = {'0': "20200722",
          '1': "20200825",
          '2': "20200904",
          '3': "20200914", #3
          '4': "20201013", #4
          '4b': "20201013b"} # premature test for #4

# Select title
title = titles['4b']

path = ""
samples = {}




def sanify(*input):
    """ Makes sure that the name of the gwf target is not illegal. """
    input = ''.join(input)
    output = []

    
    for i in str(input):
        
        ascii = ord(i)
        if (ascii >= 48 and ascii <= 57) or (ascii >= 65 and ascii <= 90) or (ascii >= 97 and ascii <= 122) or ascii == 95:
            output.append(i)
        else:
            output.append('_')

    return ''.join(output)



with open(f"input/{title}.tab") as input_file:
    for i_, line in enumerate(input_file):
        if line[0] == "#" or line == "/n":
            continue
        
        line = line.strip("\n")
        if line.startswith("path\t"):
            path = line.split("\t")[1]

        elif line.startswith("sample\t"):
            line_splitted = line.split("\t")
            if not line_splitted[1].startswith("NB"):
                raise Exception(f"Barcodes ought to start with 'NB' for sample ({line_splitted[2]}) on line {i_+1} in the input file. ")
            samples[line_splitted[1]] = line_splitted[2]


print("title:   ", title)
print(" path:    ", path)
#print(" samples: ", samples)




gwf.target(sanify('ar1_bsecal_', title),
    inputs = [path],
    outputs = [f"output/{title}/guppy_basecaller/pass/barcode{i[2:4]}" for i in samples.keys()],
    cores = 4,
    memory = '32gb',
    walltime = '24:00:00',
    queue = 'gpu', 
    gres = 'gpu:1') << f"""

        #software/ont-guppy/bin/guppy_basecaller --help > help.txt

        mkdir -p output/{title}/guppy_basecaller/

        software/ont-guppy/bin/guppy_basecaller \
        --input_path {path} \
        --recursive \
        --save_path output/{title}/guppy_basecaller \
        -q 0 \
        --gpu_runners_per_device 32 \
        --cpu_threads_per_caller 4 \
        --flowcell FLO-MIN106 --kit SQK-LSK109 \
        --barcode_kits EXP-NBD104 EXP-NBD114 \
        --require_barcodes_both_ends \
        --min_score 60.000000 \
        --qscore_filtering --min_qscore=7 \
        --device "auto" \
         1> output/{title}/guppy_basecaller/pipeline.log 2> output/{title}/guppy_basecaller/pipeline.err 

        touch output/{title}/guppy/done.flag


         """

# Demultiplexing is performed in the previous step (ar1)
#gwf.target(sanify('ar2_demux_', title),
#    inputs = [f"output/{title}/guppy_basecaller/"],
#    outputs = [f"output/{title}/guppy_barcoder/"], 
#    cores = 2, 
#    memory = '8g',
#    walltime = '01:00:00',
#    queue = 'gpu',
#    gres = 'gpu:1') << f"""
#    
#        #software/ont-guppy/bin/guppy_barcoder --help >> guppy_barcoder.txt
#        #exit 0
#
#        software/ont-guppy/bin/guppy_barcoder --require_barcodes_both_ends -i output/{title}/guppy_basecaller -s output/{title}/guppy_barcoder/ --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
#
#
#        """


for barcode, sample_name in samples.items():
    barcode_short = barcode[2:4]

    print("  " + barcode, barcode_short, sample_name, sep = "\t")


    # We use a length filter here of between 400 and 700 to remove obviously chimeric reads.
    gwf.target(sanify('ar2_gpplex_', barcode, '_', title),
        inputs = [f"output/{title}/guppy_basecaller/pass/barcode{barcode_short}"], 
        outputs = [f"output/{title}/artic_guppyplex/{sample_name}.fastq"], 
        cores = 2, 
        memory = '8g',
        walltime = '01:00:00',
        queue = 'gpu', 
        gres = 'gpu:1') << f"""
        
            {environment_base}
            conda activate artic-ncov2019


            artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory output/{title}/guppy_basecaller/pass/barcode{barcode_short} --prefix prefix_{sample_name} --output output/{title}/artic_guppyplex/{sample_name}.fastq


            """

            
    # Run the MinION pipeline
    gwf.target(sanify('ar3_articm_', barcode, '_', title),
            inputs = [f"output/{title}/artic_guppyplex/{sample_name}.fastq"],
            outputs = [f"output/{title}/artic/{sample_name}"], 
            cores = 4, 
            memory = '8g',
            walltime = '01:00:00') << f"""
            
                {environment_base}
                conda activate artic-ncov2019


                artic minion --normalise 200 --threads 4 --scheme-directory artic-ncov2019/primer_schemes --read-file output/{title}/artic_guppyplex/{sample_name}.fastq" --fast5-directory output/{title}/guppy_basecaller/pass --sequencing-summary output/{title}/sequencing_summary.txt nCoV-2019/V3 {sample_name}

                touch

                """



# # Collect consensus-files
# gwf.target(sanify('a4_collect_', title),
#         inputs = [], #[f"output/{title}/guppy_basecaller/pass/barcode{barcode_short}/"],
#         outputs = [],
#         cores = 2, 
#         memory = '8g',
#         walltime = '01:00:00',
#         queue = 'gpu', 
#         gres = 'gpu:1') << f"""
#         
#             {environment_base}
#             conda activate artic-ncov2019
# 
# 
#             artic minion --normalise 200 --threads 4 --scheme-directory ~/artic-ncov2019/primer_schemes --read-file run_name_barcode03.fastq --fast5-directory path_to_fast5 --sequencing-summary path_to_sequencing_summary.txt nCoV-2019/V3 samplename
# 
#             """

print()











