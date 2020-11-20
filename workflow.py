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
          '3': "20200914",
          '4': "20201013",
          '5': "20201030",
          '6': "20201104",
          '7': "20201104_alt"} # I copied over the basecalled data from the basement, and then the #6 job started....


# Select title
title = titles['6']



# Initialize empty variables
path = ""
samples = {}

input_file = f"input/{title}.tab"


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

check_barcode_set = set()
check_sample_name_set = set()
check_n_samples = 0

with open(input_file) as input_file_path:
    for i_, line in enumerate(input_file_path):
        if line[0] == "#" or line == "/n":
            continue
        
        line = line.strip("\n")
        if line.startswith("path\t"):
            path = line.split("\t")[1]

        elif line.startswith("sample\t"):
            line_splitted = line.split("\t")
            
            barcode = line_splitted[1]
            sample_name = line_splitted[2]
            
            check_n_samples += 1
            check_barcode_set.add(barcode)
            check_sample_name_set.add(sample_name)

            if not barcode.startswith("NB"):
                raise Exception(f"Barcodes ought to start with 'NB' for sample ({line_splitted[2]}) on line {i_+1} in the input file {input_file}. ")
            samples[line_splitted[1]] = line_splitted[2]

# Check that the labels are unique:
if len(check_barcode_set) != len(check_barcode_set) or len(check_sample_name_set) != check_n_samples:
    raise Exception(f"The number of unique barcodes, and the number of unique sample_names do not match in the input file {input_file}.")

# Debug info
print("title:   ", title)
print(" path:    ", path)




gwf.target(sanify('ar1_bsecal_', title),
    inputs = [path],
    outputs = [f"output/{title}/guppy_basecaller/sequencing_summary.txt"],
    cores = 8,
    memory = '16gb',
    walltime = '24:00:00',
    queue = 'gpu', 
    gres = 'gpu:1') << f"""

        #software/ont-guppy/bin/guppy_basecaller --help > help.txt


        mkdir -p output/{title}/guppy_basecaller/

        cp {input_file} output/{title}/{title}.samples.txt



        #nvidia-smi > smi.out 2> smi.err


        # todo: run the quality-control here.

        # 3rd idea, split into steps, better for debugging.
        software/ont-guppy/bin/guppy_basecaller \
            --num_callers 8 \
            --gpu_runners_per_device 32 \
            --chunks_per_runner 1024 \
            --chunk_size 2000 \
            --cpu_threads_per_caller 1 \
            -c dna_r9.4.1_450bps_hac.cfg \
            -i {path} \
            --save_path output/{title}/guppy_basecaller \
            -x auto \
            -r \
            --qscore_filtering --min_qscore=7 \
            2> output/{title}/guppy_basecaller/std.err 1> output/{title}/guppy_basecaller/std.out



        #touch output/{title}/guppy_basecaller/done.flag

        # Kan filtreringen være problemet?
        # --qscore_filtering --min_qscore=7 \



         """

# Demultiplexing is not necessarily performed in the previous step (ar1)
gwf.target(sanify('ar2_barcod_', title),
    inputs = [f"output/{title}/guppy_basecaller/sequencing_summary.txt"],
    outputs = [f"output/{title}/guppy_barcoder/barcode{barcode[2:4]}" for barcode in samples.keys()], 
    cores = 2, 
    memory = '8g',
    walltime = '01:00:00',
    queue = 'gpu') << f"""
    
        # software/ont-guppy/bin/guppy_barcoder --help > guppy_barcoder_help.txt
        

        mkdir -p output/{title}/guppy_barcoder

        software/ont-guppy/bin/guppy_barcoder \
            --require_barcodes_both_ends \
            -i output/{title}/guppy_basecaller/pass \
            -s output/{title}/guppy_barcoder \
            --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg" \
            --min_score 60.000000 \
            --compress_fastq \
            1> output/{title}/guppy_barcoder/std.out 2> output/{title}/guppy_barcoder/std.err

        # This is an intermediate file, which could be deleted.. 

        """


for barcode, sample_name in samples.items():
    barcode_short = barcode[2:4]

    print("  " + barcode, barcode_short, sample_name, sep = "\t")


    # We use a length filter here of between 400 and 700 to remove obviously chimeric reads.
    gwf.target(sanify('ar3_filter_', barcode, '_', title),
        inputs = [f"output/{title}/guppy_barcoder/barcode{barcode_short}"], 
        outputs = [f"output/samples/{sample_name}/filtered/{title}_{barcode}_{sample_name}.fastq"], 
        cores = 2, 
        memory = '8g',
        walltime = '01:00:00') << f"""
        
            {environment_base}
            conda activate artic-ncov2019


            mkdir -p output/samples/{sample_name}/filtered

            artic guppyplex \
                --skip-quality-check \
                --min-length 400 \
                --max-length 700 \
                --directory output/{title}/guppy_barcoder/barcode{barcode_short} \
                --prefix prefix_{sample_name} \
                --output output/samples/{sample_name}/filtered/{title}_{barcode}_{sample_name}.fastq


            """

            
    # Run the MinION pipeline
    gwf.target(sanify('ar4_articm_', barcode, '_', title),
            inputs = [path,
                      f"output/samples/{sample_name}/filtered/{title}_{barcode}_{sample_name}.fastq"],
            outputs = [f"output/samples/{sample_name}/genome/{sample_name}.consensus.fasta"], 
            cores = 8, 
            memory = '8g',
            walltime = '08:00:00') << f"""
            
                {environment_base}
                conda activate artic-ncov2019

                mkdir -p output/samples/{sample_name}/genome
                cd output/samples/{sample_name}/genome

                artic minion \
                    --normalise 200 \
                    --threads 8 \
                    --scheme-directory ../../../../artic-ncov2019/primer_schemes \
                    --read-file ../filtered/{title}_{barcode}_{sample_name}.fastq \
                    --fast5-directory ../../../../{path} \
                    --sequencing-summary ../../../{title}/guppy_basecaller/sequencing_summary.txt \
                    nCoV-2019/V3 \
                    {sample_name}
                          

                """



# Collect consensus-files
gwf.target(sanify('a4_cllect_', title),
        inputs = [f"output/samples/{sample_name}/genome/{sample_name}.consensus.fasta" for sample_name in samples.values()], #[f"output/{title}/guppy_basecaller/pass/barcode{barcode_short}/"],
        outputs = [f"output/{title}/{title}.all.fasta"],
        cores = 1, 
        memory = '1g',
        walltime = '00:10:00') << f"""
        
            # {environment_base}
            # conda activate artic-ncov2019

            cat {' '.join([f"output/samples/{sample_name}/genome/{sample_name}.consensus.fasta" for sample_name in samples.values()])} > output/{title}/{title}.all.fasta
            
            """





# Apply pangolin 
gwf.target(sanify('a5.1_pgln_', title),
    inputs = [f"output/{title}/{title}.all.fasta"],
    outputs = [f"output/{title}/{title}.pangolin.csv"],
    cores = 8,
    memory = '8g',
    walltime = '02:00:00') << f"""
        
        

        singularity run --bind output/{title}:/seq \
            docker://staphb/pangolin:2.0.4 pangolin /seq/{title}.all.fasta \
                --threads 8 \
                --outdir /seq

        mv output/{title}/lineage_report.csv output/{title}/{title}.pangolin.csv


        """


# Apply nextclade
gwf.target(sanify('a5.2_nxcl_,', title),
    inputs = [f"output/{title}/{title}.all.fasta"],
    outputs = [f"output/{title}/{title}.nextclade.csv"],
    cores = 4,
    memory = '8g',
    walltime = '02:00:00') << f"""

        # Courtesy of Marc

        singularity run --bind output/{title}:/seq \
            docker://neherlab/nextclade:0.7.5 nextclade.js \
                --input-fasta /seq/{title}.all.fasta \
                --output-csv /seq/{title}.nextclade.csv


        """


print()











