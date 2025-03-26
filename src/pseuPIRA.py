#!/usr/bin/env python

"""
pseuPIRA.py by Rohan Maddamsetti.

themisto, and minimap2 must be in the $PATH.

"""

import subprocess
import os
import gzip
from Bio import SeqIO
from os.path import basename, exists
import pprint
import polars as pl
import HTSeq ## for filtering fastq multireads.
import numpy as np ## for matrix multiplications for running PIRA.
import argparse ## for command-line interface.

"""
TODO list:

1) fix my polars dataframe code style to use parentheses, and start lines with ".join" and so forth.
Example:

filtered_df = (
    df.groupby("AnnotationAccession")
      .agg([
          (pl.col("SufficientReadCount").all()).alias("AllTrue")
      ])
      .filter(pl.col("AllTrue"))
      .join(df, on="AnnotationAccession", how="inner")
      .drop("AllTrue")
)

3) clean up code to be consistent throughout in the use of
AnnotationAccessions, RefSeq_IDs, and/or 'AnnotationAccession_genomic' as directory names.


"""

################################################################################
## Functions.

def make_replicon_fasta_references(gbk_gz_path, fasta_outdir):
    print("reading in as input:", gbk_gz_path)
    ## open the input reference genome file.
    with gzip.open(gbk_gz_path, 'rt') as gbk_gz_fh:
        SeqID = None
        SeqType = None
        for i, record in enumerate(SeqIO.parse(gbk_gz_fh, "genbank")):
            SeqID = record.id
            if "chromosome" in record.description or i == 0:
                ## IMPORTANT: we assume here that the first record is a chromosome.
                SeqType = "chromosome"
            elif "plasmid" in record.description:
                SeqType = "plasmid"
            else:
                continue
            ## replace spaces with underscores in the replicon annotation field.
            replicon_description = record.description.replace(" ","_")
            header = ">" + "|".join(["SeqID="+SeqID,"SeqType="+SeqType,"replicon="+replicon_description])
            my_replicon_fastafile = SeqID + ".fna"
            my_replicon_outfilepath = os.path.join(fasta_outdir, my_replicon_fastafile)
            with open(my_replicon_outfilepath, "w") as outfh:
                outfh.write(header + "\n")
                outfh.write(str(record.seq) + "\n")
    return


def make_fasta_replicon_list_file(genome_id, fasta_outdir):
    """
    IMPORTANT: the replicons in the fasta reference list file need to be sorted by length.
    Directly measure the length of each of the replicons here, and then use this to sort
    the list of replicon FASTA paths, before writing to file.
    """

    ## IMPORTANT: exclude any FASTA file that has the genome_id as a name.
    ## (this is a file containing all the replicons, created and used in the PIRA stages downstream).
    replicon_fasta_filelist = [x for x in os.listdir(fasta_outdir) if x.endswith(".fna") and not x.startswith(genome_id)]
    replicon_fasta_pathlist = [os.path.join(fasta_outdir, x) for x in replicon_fasta_filelist]

    ## Decorate-sort-undecorate by fasta sequence length.
    decorated_replicon_fasta_pathlist = list()
    for fastapath in replicon_fasta_pathlist:
        my_replicon = SeqIO.read(fastapath, "fasta")
        ## Get the length of the replicon
        replicon_length = len(my_replicon.seq)
        ## append a tuple of (fastapath, replicon_length)
        my_tuple = (fastapath, replicon_length)
        decorated_replicon_fasta_pathlist.append(my_tuple)
    ## sort the decorated list in place by replicon_length, in descending order from
    ## from largest to smallest replicon.
    decorated_replicon_fasta_pathlist.sort(key=lambda x: x[1], reverse=True)
    ## undecorate the path list, which is now in sorted order.
    sorted_replicon_fasta_pathlist = [x[0] for x in decorated_replicon_fasta_pathlist]

    ## write the path list to file for themisto.
    replicon_listfile = os.path.join(fasta_outdir, genome_id + ".txt")
    with open(replicon_listfile, "w") as fastatxtfile_fh:
        for fasta_path in sorted_replicon_fasta_pathlist:
            fastatxtfile_fh.write(fasta_path + "\n")
    return


def make_replicon_fasta_refs(genome_id, refgenome_gzpath, fasta_outdir):
    ## this function makes a genome directory for each genome.
    ## each directory contains separate fasta files for each replicon.

    assert refgenome_gzpath.endswith("gbff.gz")
    
    ## make the fasta output directory if it does not exist.
    if not exists(fasta_outdir):
        os.mkdir(fasta_outdir)
    
    make_replicon_fasta_references(refgenome_gzpath, fasta_outdir)
    make_fasta_replicon_list_file(genome_id, fasta_outdir)
    return


def run_command_with_retry(command_string, tempdir=None, max_retries=3, timeout=20):
    ## This code handles a bug in themisto build-- sometimes randomly hangs, have to delete temp files
    ## and restart and then it usually works.
    retries = 0
    while retries < max_retries:
        process = subprocess.Popen(command_string, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        timer = threading.Timer(timeout, process.kill)  ## Kill process if it exceeds timeout

        try:
            timer.start()
            stdout, stderr = process.communicate()
        finally:
            timer.cancel()

        if process.returncode == 0:
            print("Command succeeded:", stdout.decode())
            return stdout.decode()
        else:
            print(f"*********COMMAND FAILED (attempt {retries + 1}):", stderr.decode())
            if tempdir is not None: ## remove temporary files from the failed run.
                print(f"removing {tempdir}")
                subprocess.run(f"rm -rf {tempdir}", shell=True)
                print(f"remaking {tempdir} before restarting")
                os.mkdir(tempdir)
            retries += 1
            time.sleep(0.1)  ## Small delay before retrying
    
    print("Command failed after maximum retries.")
    return


def build_genome_index(genome_id, ref_fasta_dir, genome_index_dir, nthreads="4"):
    
    ## make sure that this path is real.
    assert os.path.isdir(ref_fasta_dir)

    ## make the genome index directory if it does not exist.
    if not exists(genome_index_dir):
        os.mkdir(genome_index_dir)
    
    index_input_filelist = os.path.join(ref_fasta_dir, genome_id + ".txt")

    ## set the index_prefix to write index files into the themisto_index_dir.
    index_prefix = os.path.join(genome_index_dir, genome_id)

    ## make the temp directory if it doesn't exist.
    tempdir = os.path.join(genome_index_dir, "temp")
    if not exists(tempdir):
        os.mkdir(tempdir)

    themisto_build_args = ["themisto", "build", "-k","31", "-i", index_input_filelist, "--index-prefix", index_prefix, "--temp-dir", tempdir, "--mem-gigas", "4", "--n-threads", nthreads, "--file-colors"]
    themisto_build_string = " ".join(themisto_build_args)

    print(themisto_build_string)
    ## if themisto build hangs for a long time, then kill and restart.
    run_command_with_retry(themisto_build_string, tempdir)

    return


def pseudoalign_reads(genome_id, index_dir, readpath_list, pseudoalignment_dir, nthreads="4"):

    assert os.path.isdir(index_dir) ## make sure the index directory exists.

    ## make the output directory if it does not exist.
    if not exists(pseudoalignment_dir):
        os.mkdir(pseudoalignment_dir)

    ## make a file containing the paths to the sequencing read data
    ## for themisto pseudoalign.
    seqdata_listfile = os.path.join(pseudoalignment_dir, genome_id + "_reads.txt")
    with open(seqdata_listfile, "w") as seqdata_txtfile_fh:
        for readpath in readpath_list:
            seqdata_txtfile_fh.write(readpath + "\n")

    ## make the temp directory if it doesn't exist.
    tempdir = os.path.join(pseudoalignment_dir, "temp")
    if not exists(tempdir):
        os.mkdir(tempdir)

    ## we make corresponding pseudoalignment output files for each sequencing dataset.
    ## This list goes into the output listfile.
    output_listfile = os.path.join(pseudoalignment_dir, genome_id + "_pseudoalignments.txt")
    with open(output_listfile, "w") as output_listfile_fh:
        for readpath in readpath_list:
            read_filename = os.path.basename(readpath).split(".fastq")[0]
            output_filename = os.path.join(pseudoalignment_dir, read_filename + "_pseudoalignment.txt")
            output_listfile_fh.write(output_filename + "\n")

    ## make the arguments relating to the paths to the index files for this genome.
    my_index_prefix = os.path.join(index_dir, genome_id)
            
    ## now run themisto pseudoalign.
    themisto_pseudoalign_args = ["themisto", "pseudoalign", "--query-file-list", seqdata_listfile, "--index-prefix", my_index_prefix, "--temp-dir", tempdir, "--out-file-list", output_listfile, "--n-threads", nthreads, "--threshold", "0.7"]
    themisto_pseudoalign_string = " ".join(themisto_pseudoalign_args)
    subprocess.run(themisto_pseudoalign_string, shell=True)
    return


def map_themisto_IDs_to_replicon_metadata(genome_id, replicon_ref_dir):
    ## now, let's map the themisto replicon ID numbers to a (SeqID, SeqType, SeqLength) tuple.
    themisto_ID_to_seq_metadata_dict = dict()
    my_themisto_replicon_ID_mapping_file = os.path.join(replicon_ref_dir, genome_id + ".txt")
    with open(my_themisto_replicon_ID_mapping_file, "r") as replicon_ID_mapping_fh:
        for i, fasta_file_path in enumerate(replicon_ID_mapping_fh):
            fasta_file_path = fasta_file_path.strip()
            ## get the header from this fasta file.
            with open(fasta_file_path, "r") as fasta_fh:
                my_header = fasta_fh.readline().strip()
                my_seq = fasta_fh.readline().strip()
            fields = my_header.strip(">").split("|")
            SeqID = fields[0].split("=")[-1]
            SeqType = fields[1].split("=")[-1]
            SeqLength = len(my_seq)
            themisto_ID_to_seq_metadata_dict[i] = (SeqID, SeqType, SeqLength)
    return themisto_ID_to_seq_metadata_dict


def summarize_pseudoalignment_results(genome_id, themisto_ID_to_seq_metadata_dict, replicon_ref_dir, pseudoalignment_dir, themisto_results_csvfile_path):
    
    with open(themisto_results_csvfile_path, "w") as output_csv_fh:
        ## first, write the output header.
        output_header = "AnnotationAccession,SeqID,SeqType,ReadCount"
        output_csv_fh.write(output_header + "\n")

        ## initialize a dictionary to store the pseudoalignment counts.
        pseudoalignment_read_count_dict = dict()
        pseudoalignment_filepaths = [os.path.join(pseudoalignment_dir, x) for x in os.listdir(pseudoalignment_dir) if x.endswith("_pseudoalignment.txt")]
        ## now, summarize the read counts for this genome.
        for my_pseudoalignment_filepath in pseudoalignment_filepaths:
            with open(my_pseudoalignment_filepath) as my_pseudoalignment_fh:
                for line in my_pseudoalignment_fh:
                    ## handle odd behavior in themisto: we need to sort the themisto replicon ID numbers ourselves.
                    ## IMPORTANT: sort the themisto ids (strings) by their numerical value.
                    replicon_set_string = " ".join(sorted(line.strip().split()[1:], key=lambda x: int(x)))
                    if replicon_set_string in pseudoalignment_read_count_dict:
                        pseudoalignment_read_count_dict[replicon_set_string] += 1
                    else:
                        pseudoalignment_read_count_dict[replicon_set_string] = 1

        ## now write the pseudoalignment counts to file.
        for replicon_set_string in sorted(pseudoalignment_read_count_dict.keys()):
            read_count = pseudoalignment_read_count_dict[replicon_set_string]
            replicon_ID_list = replicon_set_string.split()
            if len(replicon_ID_list) == 0:
                SeqID = "NA"
                SeqType = "NA"
            elif len(replicon_ID_list) == 1:
                my_replicon_ID = replicon_ID_list[0]
                ## the SeqLength parameter in themisto_ID_to_seq_metadata_dict is not used here.
                SeqID, SeqType, _ = themisto_ID_to_seq_metadata_dict[int(my_replicon_ID)]
            else:
                SeqID = "&".join([themisto_ID_to_seq_metadata_dict[int(replicon_ID)][0] for replicon_ID in replicon_ID_list])
                SeqType = "multireplicon_sequence"
            ## now write to file.
            rowdata = ",".join([genome_id, SeqID, SeqType, str(read_count)])
            output_csv_fh.write(rowdata + "\n")
    return


def make_replicon_metadata_dataframe(genome_id, themisto_ID_to_seq_metadata_dict):
    ## reformat metadata into a polars DataFrame.
    annotation_accession_col  = list()
    themisto_id_col = list()
    SeqID_col = list()
    SeqType_col = list()
    replicon_length_col = list() ## using this name instead of "SeqLength" to match current PIRA code.
    
    for replicon_ID in sorted(themisto_ID_to_seq_metadata_dict.keys()):
        SeqID, SeqType, replicon_length = themisto_ID_to_seq_metadata_dict[replicon_ID]
        annotation_accession_col.append(genome_id)
        themisto_id_col.append(replicon_ID)
        SeqID_col.append(SeqID)
        SeqType_col.append(SeqType)
        replicon_length_col.append(replicon_length)
        
    replicon_metadata_dict = {
        'AnnotationAccession' : annotation_accession_col,
        'ThemistoID' : themisto_id_col,
        'SeqID' : SeqID_col,
        'SeqType' : SeqType_col,
        'replicon_length' : replicon_length_col}
    
    replicon_metadata_df = pl.DataFrame(replicon_metadata_dict)
    return replicon_metadata_df


def naive_PCN_estimation(replicon_metadata_df, themisto_results_csv_file, naive_PCN_csv_file):
    ## This function simply ignores multireplicon reads when estimating PCN.
    ##Also note that this function omits replicons with zero mapped reads.
    print("running naive themisto PCN estimation (ignoring multireplicon reads)")
    ## import the data as polars dataframes.
    
    naive_themisto_read_count_df = pl.read_csv(themisto_results_csv_file).filter(
        (pl.col("SeqType") == "chromosome") | (pl.col("SeqType") == "plasmid")).join(
            replicon_metadata_df, on = "SeqID").with_columns(
                (pl.col("ReadCount") / pl.col("replicon_length")).alias("SequencingCoverage"))
    
    ## make a second dataframe containing just the sequencing coverage for the longest replicon for each genome.
    ## to do so, first group by AnnotationAccession and compute maximum replicon_length within each group.
    longest_replicon_df = naive_themisto_read_count_df.group_by(
        "AnnotationAccession").agg(pl.col("replicon_length").max()).join(
            ## now join with the original DataFrame to filter for rows with the maximum replicon_length
            naive_themisto_read_count_df, on=["AnnotationAccession", "replicon_length"], how="inner").select(
                pl.col("AnnotationAccession", "SequencingCoverage")).with_columns(
                    pl.col("SequencingCoverage").alias("LongestRepliconCoverage")).select(
                        pl.col("AnnotationAccession", "LongestRepliconCoverage"))

    ## now normalize SequencingCoverage by LongestRepliconCoverage for each genome to calculate PCN.
    naive_themisto_PCN_df = naive_themisto_read_count_df.join(
        longest_replicon_df, on = "AnnotationAccession", coalesce=True).with_columns(
            (pl.col("SequencingCoverage") / pl.col("LongestRepliconCoverage")).alias("CopyNumber"))

    ## now write the naive PCN estimates to file.
    print(naive_themisto_PCN_df)
    naive_themisto_PCN_df.write_csv(naive_PCN_csv_file)    
    return


def assign_multireplicon_reads(genome_df):
    ## create a new data frame with just chromosome and plasmid data.
    updated_df = genome_df.filter(
        (pl.col("SeqType") == "chromosome") | (pl.col("SeqType") == "plasmid")).with_columns(
            ## cast the ReadCount to floats, and
            ## replace null values in the ReadCount columns with floating-point zeros.
            pl.col("ReadCount").cast(pl.Float64).fill_null(pl.lit(0.0)))

    ## get the multireplicon data.
    multireplicon_df = genome_df.filter(pl.col("SeqType") == "multireplicon_sequence")
    ## iterate over each multireplicon read set.
    for row_dict in multireplicon_df.iter_rows(named=True):
        seq_id_list = row_dict["SeqID"].split("&")
        read_count = row_dict["ReadCount"]
        num_reads_to_assign_to_each_replicon = float(read_count) / float(len(seq_id_list))
        ## iterate over each replicon in this multireplicon read set.
        for seq_id in seq_id_list:
            ## update the relevant values in the dataframe
            old_seqID_readcount = updated_df.filter(updated_df["SeqID"] == seq_id)["ReadCount"][0]
            new_seqID_readcount = old_seqID_readcount + num_reads_to_assign_to_each_replicon
            ## a kludgy hacky solution, but works:
            ## create a new column called temp, that has the new_seqID_readcount for the given SeqID,
            ## and zeros elsewhere in the column.
            ## then make a new column called updatedReadCount that takes the maximum of ReadCount and temp.
            ## then drop ReadCount and temp and rename updatedReadCount as ReadCount.
            updated_df = updated_df.with_columns(
                "ReadCount",
                pl.when(pl.col("SeqID") == seq_id).then(new_seqID_readcount).otherwise(pl.lit(0.0)).alias("temp")
            ).with_columns(pl.max_horizontal("ReadCount", "temp").alias("updatedReadCount")).select(
                pl.col("*").exclude(["ReadCount", "temp"])).rename({"updatedReadCount":"ReadCount"})
    return updated_df


def filter_fastq_files_for_multireads(multiread_data_dir, pseudoalignment_dir, readpath_list):
    ## writes filtered multireads in vanilla .fastq format.
    
    ## make the directory for filtered multireads  if it does not exist.
    if not exists(multiread_data_dir):
        os.mkdir(multiread_data_dir)

    for my_readpath in readpath_list:
        ## before sorting the pseudoalignment results, check whether the filtered multireads
        ## for this pseudoalignment already exist on disk. If so, then skip--
        ## don't repeat work that has already been done.

        ## construct the path to the filtered fastq file.
        ## note: multiread_genome_dir is only created if filtered multireads exist.
        ## this means that this path only exists on disk if filtered multireads for this genome
        ## were already written to disk in a previous call of this function.

        ## NOTE: input may be either *.fastq or gzipped *.fastq.gz data,
        ## but filtered reads are always just *.fastq format.
        my_readfile = basename(my_readpath)
        my_filtered_fastq_file = "multireads_" + my_readfile.replace(".gz", "")
        my_filtered_fastq_path = os.path.join(multiread_data_dir, my_filtered_fastq_file)
        ## now check to see if the filtered fastq file already exists.
        if exists(my_filtered_fastq_path): continue ## if so, then don't repeat the work.

        my_pseudoalignment_file = basename(my_readpath).split(".fastq")[0] + "_pseudoalignment.txt"
        my_pseudoalignment_path = os.path.join(pseudoalignment_dir, my_pseudoalignment_file)
        
        list_of_multiread_tuples = list()
        with open(my_pseudoalignment_path, "r") as pseudoalignment_fh:
            for line in pseudoalignment_fh:
                line = line.strip() ## remove whitespace
                fields = line.split()
                if len(fields) <= 2: continue ## skip reads that map to zero or one reference sequences.
                read_index = int(fields[0])
                matched_reference_tuple = tuple(fields[1:])
                multiread_tuple = (read_index, matched_reference_tuple)
                list_of_multiread_tuples.append(multiread_tuple)
        ## sort the list of multiread tuples in-place by the read_index.
        list_of_multiread_tuples.sort(key=lambda x: x[0])
        ## if there are multireads, then filter the original fastq file for multireads.
        if len(list_of_multiread_tuples):
            ## if there are multireads, then make multiread_data_dir for this genome.
            if not exists(multiread_data_dir):
                os.mkdir(multiread_data_dir)
            ## IMPORTANT: read indices in the themisto pseudoalignments are zero-based (first index is 0).
            multiread_indices = {x[0] for x in list_of_multiread_tuples} ## this is a set
            print(f"filtering {my_readfile} for multireads. multireads written into {my_filtered_fastq_path}")
            ## write out the filtered reads
            with open(my_filtered_fastq_path, "w") as filtered_fastq_fh:
                my_fastq_reader = HTSeq.FastqReader(my_readpath)
                for i, read in enumerate(my_fastq_reader):
                    ## skip reads that aren't in the set of multireads
                    if i not in multiread_indices: continue
                    read.write_to_fastq_file(filtered_fastq_fh)
    return


def make_fasta_reference_genomes_for_minimap2(genome_id, replicon_ref_dir):
    
    ## get the paths to each replicon FASTA file.
    ## IMPORTANT: the order of replicons in this file MATTERS.
    ## This corresponds to the zero-indexed replicon assignment done by themisto pseudoalignment.
    replicon_fasta_pathlist = []
    my_replicon_fasta_listfilepath = os.path.join(replicon_ref_dir, genome_id + ".txt")
    with open(my_replicon_fasta_listfilepath, "r") as replicon_listfile_fh:
        for line in replicon_listfile_fh:
            line = line.strip()
            replicon_fasta_pathlist.append(line)

    my_fasta_reference_genome_outfile = genome_id + ".fna"
    my_fasta_reference_genome_outpath = os.path.join(replicon_ref_dir, my_fasta_reference_genome_outfile)

    with open(my_fasta_reference_genome_outpath, "w") as fasta_outfh:
        for themisto_replicon_num, replicon_fasta_path in enumerate(replicon_fasta_pathlist):
            with open(replicon_fasta_path, "r") as my_fasta_infh:
                for i, line in enumerate(my_fasta_infh):
                    ## add the replicon ID number assigned by themisto to the FASTA header
                    if i == 0:
                        header_string = line.lstrip(">")
                        replicon_id_string = "ThemistoRepliconID=" + str(themisto_replicon_num)
                        updated_header = ">" +  replicon_id_string + "|" + header_string
                        fasta_outfh.write(updated_header)
                    else:
                        fasta_outfh.write(line)
    return


def align_multireads_with_minimap2(genome_id, replicon_ref_dir, multiread_data_dir, multiread_alignment_dir):
    ## Example: minimap2 -x sr ref.fa read1.fq read2.fq > aln.paf
    ## Details of PAF format are here: https://github.com/lh3/miniasm/blob/master/PAF.md

    ## make the directory for multiread alignments  if it does not exist.
    if not exists(multiread_alignment_dir):
        os.mkdir(multiread_alignment_dir)

    ## make sure there are no weird files like .DS_Store in this list!
    multiread_data_pathlist = [os.path.join(multiread_data_dir, x) for x in os.listdir(multiread_data_dir) if x.endswith(".fastq")]
    multiread_data_pathlist.sort() ## sort the filtered data (say if there are paired read files, etc.)

    if len(multiread_data_pathlist) == 0: ## no multireads to align.
        print("no multireads for PIRA. Naive estimates == PIRA estimates.")
        print("minimap2 alignment step skipped.")
        return

    ref_genome_fasta_file = genome_id + ".fna"
    reference_genome_path = os.path.join(replicon_ref_dir, ref_genome_fasta_file)

    ## run single-end read alignment on each fastq file separately.
    for my_index, cur_multiread_data_path in enumerate(multiread_data_pathlist):
        my_alignment_file = basename(cur_multiread_data_path).split(".fastq")[0] + ".paf"
        my_multiread_alignment_outpath = os.path.join(multiread_alignment_dir, my_alignment_file)

        minimap2_cmd_string = " ".join(["minimap2 -x sr", reference_genome_path, cur_multiread_data_path, ">", my_multiread_alignment_outpath])
        print(minimap2_cmd_string)
        subprocess.run(minimap2_cmd_string, shell=True)
    return


def parse_read_alignments(genome_dir):
    ## make a dictionary from reads to the multiset of replicons that the read maps to.
    paf_alignment_files = [x for x in os.listdir(genome_dir) if x.endswith(".paf")]
    paf_alignment_paths = [os.path.join(genome_dir, x) for x in paf_alignment_files]
    
    read_mapping_dict = dict()
    for paf_path in paf_alignment_paths:
        with open(paf_path, "r") as paf_fh:
            for line in paf_fh:
                line = line.strip()
                fields = line.split("\t")
                ## see PAF format documentation here:
                ## https://github.com/lh3/miniasm/blob/master/PAF.md
                read_name = fields[0]
                themisto_replicon_ID = int(fields[5].split("|")[0].split("=")[-1]) ## IMPORTANT: turn into an integer
                
                if read_name in read_mapping_dict:
                    read_mapping_dict[read_name].append(themisto_replicon_ID)
                else:
                    read_mapping_dict[read_name] = [themisto_replicon_ID]
    return read_mapping_dict


def initialize_GenomeDataFrame(themisto_ID_to_seq_metadata_dict):
    """ I use themisto_ID_to_seq_metadata_dict to initialize the dataframe for a genome.
     themisto_ID_to_seq_metadata_dict contains metadata for all replicons in the genome,
     regardless of whether any reads mapped to that replicon.
     To handle the cases that:
     1) a replicon only contains multireads, therefore it's naive PCN == 0, or
     2) a replicon does not have any multireads, therefore, it is not present in additional_replicon_reads_dict, or
     3) no replicons have any unireads mapping to them, so my_naive_themisto_PCN_df is empty.
     we use themisto_ID_to_seq_metadata_dict to initialize a basic polars dataframe containing
     all replicons in our genome.
    """
    themisto_id_col = sorted(themisto_ID_to_seq_metadata_dict.keys())
    replicon_seq_id_col = list()
    replicon_seq_type_col = list()
    replicon_length_col = list()

    for themisto_id in themisto_id_col:
        ## get the SeqID and SeqType given the Themisto ID.
        seq_id, seq_type, seq_length = themisto_ID_to_seq_metadata_dict[themisto_id]
        replicon_seq_id_col.append(seq_id)
        replicon_seq_type_col.append(seq_type)
        replicon_length_col.append(seq_length)

    genome_df = pl.DataFrame(
        {"ThemistoID" : themisto_id_col,
         "SeqID" : replicon_seq_id_col,
         "SeqType" : replicon_seq_type_col,
         "replicon_length" : replicon_length_col})

    return genome_df


def make_PIRAGenomeDataFrame(
        additional_replicon_reads_dict,
        themisto_ID_to_seq_metadata_dict,
        my_naive_themisto_PCN_df):
    """ Make the DataFrame with the data needed for PIRA on a given genome.
        We have to update the results of the Naive PCN estimates from Themisto
        (results of stage 16) by adding the additional replicon reads found by re-aligning multireads
        with minimap2.
    """

    ## Turn the additional_replicon_reads_dict into a polars dataframe.
    ## First initialize the dataframe to contain all replicons in the genome.
    additional_replicon_reads_df = initialize_GenomeDataFrame(themisto_ID_to_seq_metadata_dict)
    ## initialize the AdditionalReadCount Column with zeros for each row in the initialized DataFrame
    additional_readcount_col = [0 for i in range(additional_replicon_reads_df.shape[0])]        
    ## now update the values in additional_readcount_col using additional_replicon_reads_dict.
    for themisto_id, readcount in additional_replicon_reads_dict.items():
        additional_readcount_col[themisto_id] = readcount
    
    ## now add the AdditionalReadCount Column to the DataFrame.
    additional_readcount_series = pl.Series("AdditionalReadCount", additional_readcount_col)
    additional_replicon_reads_df = additional_replicon_reads_df.with_columns([additional_readcount_series])

    ## IMPORTANT: the previous code ensures that additional_replicon_reads_df
    ## contains ALL replicons in the genome, even if no reads pseudoaligned or aligned to that replicon.
    ## this may happen for short contigs that are annotated as plasmids or plasmid fragments in the
    ## reference genome (example: genome GCF_017654585.1_ASM1765458v1).

    ## To fill in missing values in the AnnotationAccession column after the merge with additional_replicon_reads_df,
    ## get the unique value in the AnnotationAccession column in my_naive_themisto_PCN_df, and use this
    ## to fill in the missing values.
    my_AnnotationAccession = my_naive_themisto_PCN_df.select(pl.col('AnnotationAccession').unique())
    assert len(my_AnnotationAccession) == 1 ## make sure this value is unique!
    
    ## merge the DataFrames containing the ReadCounts.
    merged_readcount_df = additional_replicon_reads_df.join(
        ## IMPORTANT: my_naive_themisto_PCN_df may not contain rows for replicons that didn't have
        ## any reads pseudoalign to it. therefore, we need to left_join my_native_themisto_PCN_df to
        ## additional_replicon_reads_df, which contains the data for ALL replicons in the genome,
        ## even if the Count is zero.
        my_naive_themisto_PCN_df, on="SeqID", how="left", coalesce=True).with_columns(
            ##fill in missing values in the AnnotationAccession column after the merge.
            pl.col("AnnotationAccession").fill_null(my_AnnotationAccession)).with_columns(
                    ## set missing values in the InitialReadCount column to 0.
                    pl.col("InitialReadCount").fill_null(strategy="zero")).with_columns(
                        ## sum those ReadCounts,
                        (pl.col("InitialReadCount") + pl.col("AdditionalReadCount")).alias("ReadCount")).with_columns(
                            ## and re-calculate SequencingCoverage,
                            (pl.col("ReadCount") / pl.col("replicon_length")).alias("SequencingCoverage"))
    
    ## The following is a hack, following code in the function native_themisto_PCN_estimation(),
    ## to recalculate LongestRepliconCoverage and CopyNumber with the additional reads.

    ## find the length of the longest replicon
    max_replicon_length = merged_readcount_df.select(pl.col("replicon_length").max()).item()

    ## filter the DataFrame to get the row for the longest replicon
    longest_replicon_row_df = merged_readcount_df.filter(
        pl.col("replicon_length") == max_replicon_length).with_columns(
            ## and define the LongestRepliconCoverage column for merging back.
            pl.col("SequencingCoverage").alias("LongestRepliconCoverage")).select(
            pl.col("AnnotationAccession", "LongestRepliconCoverage"))

    ## now normalize SequencingCoverage by LongestRepliconCoverage for each genome to calculate PCN.
    PIRAGenomeDataFrame = merged_readcount_df.join(
        longest_replicon_row_df, on = "AnnotationAccession", coalesce=True).with_columns(
        (pl.col("SequencingCoverage") / pl.col("LongestRepliconCoverage")).alias("InitialCopyNumberEstimate")).sort(
            ## and sort by the ThemistoID column,
            "ThemistoID").select(
                ## super annoying, not sure why the AnnotationAccession_right columns
                ## are kept. This select command removes these redundancies.
                ['ThemistoID', 'SeqID', 'SeqType', 'replicon_length',
                 'AdditionalReadCount', 'AnnotationAccession', 'InitialReadCount',
                 'ReadCount', 'SequencingCoverage', 'LongestRepliconCoverage', 'InitialCopyNumberEstimate']
            )
    
    """
    The use of the function make_fasta_replicon_list_file(genome_id, fasta_outdir)
    ensures that Themisto Replicon IDs are sorted by replicon length in descending order
    (i.e., Replicon ID == 0 corresponds to the longest chromosome).
    """
    return PIRAGenomeDataFrame


def initializePIRA(multiread_mapping_dict, themisto_ID_to_seq_metadata_dict, my_naive_themisto_PCN_df):
    ## Iterate over the multiread_mapping_dict, and split into two data structures:
    ## 1) reads that map to a single replicon are counted up in a dictionary.
    additional_replicon_reads_dict = dict()
    ## 2) reads that map to multiple replicons are stored in a list of lists,
    ## that will then be turned into the numpy array to store the match matrix M.
    match_matrix_list_of_rows = list()

    for read, replicon_list in multiread_mapping_dict.items():
        replicon_set = set(replicon_list)
        if len(replicon_set) == 0: ## the read does not map to any replicons -- should not occur.
            raise AssertionError(f"READ {read} DID NOT ALIGN TO ANY REPLICON")
        elif len(replicon_set) == 1: ## the read maps to a single replicon.
            my_replicon = replicon_set.pop()
            if my_replicon in additional_replicon_reads_dict:
                additional_replicon_reads_dict[my_replicon] += 1
            else:
                additional_replicon_reads_dict[my_replicon] = 1
        else: ## the read maps to multiple replicons
            ## initialize a row of zeros, based on the number of replicons in this genome.
            match_matrix_rowlist = [0 for k in themisto_ID_to_seq_metadata_dict.keys()]
            for replicon_index in replicon_list:
                match_matrix_rowlist[replicon_index] += 1
            match_matrix_list_of_rows.append(match_matrix_rowlist)

    ## now set up the data structures for PIRA on this genome.
    MatchMatrix = np.array(match_matrix_list_of_rows)
    
    """ Generate the DataFrame containing the ReadCounts,  replicon lengths, and initial PCN estimates.
    We update the results of the Naive PCN estimates from Themisto (results of stage 16)
    by adding the additional replicon reads found by re-aligning multireads
    with minimap2.
    """
    PIRAGenomeDataFrame = make_PIRAGenomeDataFrame(
        additional_replicon_reads_dict, themisto_ID_to_seq_metadata_dict, my_naive_themisto_PCN_df)
    ## return inputs requires for PIRA for this genome.
    return (MatchMatrix, PIRAGenomeDataFrame)


def run_PIRA(M, PIRAGenomeDataFrame, epsilon = 0.00001):
    ## Run PIRA for a genome, assuming that the zero-th index of the match matrix M is the chromosome for normalization.
    print("RUNNING PIRA.")
    print(M)
    print(M.shape)
    print()

    print(PIRAGenomeDataFrame)
    print()
    print(PIRAGenomeDataFrame.glimpse(max_items_per_column=100))
    print()

    """
    Stage 12 calls generate_replicon_fasta_reference_list_file_for_themisto(fasta_outdir), which
    ensures that Themisto Replicon IDs are sorted by replicon length in descending order
    (i.e., Replicon ID == 0 corresponds to the longest chromosome).

    Therefore, we can  extract the InitialCopyNumber Column as a numpy vector as the initial PCN estimate guess for PIRA.
    """
    
    v = PIRAGenomeDataFrame["InitialCopyNumberEstimate"].to_numpy()
    readcount_vector = PIRAGenomeDataFrame["ReadCount"].to_numpy()
    replicon_length_vector = PIRAGenomeDataFrame["replicon_length"].to_numpy()
    
    """
    M may be empty, in the case that all multireads uniquely mapped to chromosome or plasmid.
    In this case, the PIRAGenomeDataFrame has already incorporated all available information,
    so we should just return v, the vector containing the InitialCopyNumberEstimate.

    I handle this logic by only running the PIRA loop if M is non-empty.
    """ 
    
    if M.shape[0] > 0: ## only run PIRA if the Match Matrix M has rows.
        convergence_error = 10000.0
        ## Iterate PIRA until error converges to zero.
        while convergence_error > epsilon:
            print(f"current convergence error: {convergence_error}")
            print(f"current PCN estimate vector: {v}")
            ## Weight M by PCN guess v-- need to turn v into a diagonal matrix first.
            diagonal_v = np.diag(v)
            weighted_M = np.matmul(M, diagonal_v)
            print(weighted_M)

            ## Normalize rows of weighted_M to sum to 1: this the probabilistic read assignment.
            ## Compute the sum of each row
            weighted_M_row_sums = weighted_M.sum(axis=1)
            ## Normalize each row by its sum
            normalized_weighted_M = weighted_M / weighted_M_row_sums[:, np.newaxis]

            print(normalized_weighted_M)
            ## sum over the rows of the normalized and weighted M matrix to generate the
            ## multiread vector.
            multiread_vector = normalized_weighted_M.sum(axis=0)
            print(multiread_vector)

            ## update the PCN estimate vector v using the multiread vector
            unnormalized_v = (multiread_vector + readcount_vector) / replicon_length_vector
            normalization_factor = unnormalized_v[0] ## coverage for the longest chromosome.
            updated_v = unnormalized_v / normalization_factor

            ## update the error estimate
            convergence_error = abs(sum(updated_v - v))
            ## and update the PCN estimate vector v
            v = updated_v
        print(f"final convergence error: {convergence_error}")

    print(f"final PCN estimate vector: {v}")
    return v


def run_PIRA_test_suite():
    ## Test 1: a simple test to check for convergence from very bad initial PCN estimates.
    Test1_DataFrame = pl.DataFrame(
        {
            "replicon_length" : [1000000, 100000, 1000], ## log 10: 6, 5, 3
            "ReadCount" : [1000000, 1000000, 1000000], ## log 10: 6, 6, 6
            "InitialCopyNumberEstimate" : [1, 1, 1] ## should converge to: 1, 10, 1000
        }
    )

    ## Test 1: the Match Matrix is negligible.
    multi_read_row1 = [1, 1, 0]
    test1_match_matrix_list_of_rows = [multi_read_row1]
    Test1_M = np.array(test1_match_matrix_list_of_rows)
    print("*"*80)
    print("PIRA TEST 1: check convergence from bad initial PCN estimates given negligible match matrix")
    run_PIRA(Test1_M, Test1_DataFrame)
    print()


    ## Test 2-- check what happens if the initial PCN vector contains zeros-- can estimates be updated
    ## properly?
    Test2_DataFrame = pl.DataFrame(
        {
            "replicon_length" : [1000000, 100000, 1000], ## log 10: 6, 5, 3
            "ReadCount" : [1000000, 1000000, 1000000], ## log 10: 6, 6, 6
            "InitialCopyNumberEstimate" : [1, 0, 0] ## should converge to: 1, 10, 1000
        }
    )

    ## Test 2: the Match Matrix is negligible, and the initial PCN estimate contains zeros
    print("*"*80)
    print("PIRA TEST 2: check convergence from bad initial PCN estimates with zeros, given negligible match matrix")
    run_PIRA(Test1_M, Test2_DataFrame)
    print()
    
    ## Test 3: Case that Match Matrix is super important for correct PCN estimation.
    Test3_DataFrame = pl.DataFrame(
        {
            "replicon_length" : [1000000, 100000, 1000], ## log 10: 6, 5, 3
            "ReadCount" : [1000000, 0, 0], ## log 10: 6, 6, 6
            "InitialCopyNumberEstimate" : [1, 1, 1] ## should converge to: 1, 10, 1000
        }
    )

    ## Test 3: the Match Matrix is super important for accurate PCN estimation.
    multi_read_row3 = [0, 1, 1]
    test3_match_matrix_list_of_rows = [multi_read_row3] * 1000000
    Test3_M = np.array(test3_match_matrix_list_of_rows)
    print("*"*80)
    print("PIRA TEST 3: check convergence when match matrix needed for accurate PCN estimation")
    run_PIRA(Test3_M, Test3_DataFrame)
    print()
    
    return


def run_PIRA_on_the_genome(genome_ID, multiread_alignment_dir, replicon_ref_dir, naive_themisto_PCN_csv_file, PIRA_PCN_csv_file):
 
    ## only run PIRA if the genome has multireads.
    if len(os.listdir(multiread_alignment_dir)) == 0:
        print("no multiread alignments for PIRA. Naive estimates == PIRA estimates.")
        print("PIRA step skipped. check the naive PCN estimates file for the final PCN estimates.")
        return

    ## map the themisto replicon ID numbers to a (SeqID, SeqType, SeqLength) tuple.
    themisto_ID_to_seq_metadata_dict = map_themisto_IDs_to_replicon_metadata(genome_ID, replicon_ref_dir)
    ## import the results of the Naive PCN estimates.

    ## IMPORTANT: This DataFrame will NOT contain rows for replicons that didn't have
    ## any reads pseudoalign to it, and in some cases, may be completely empty.
    ##therefore, we need to join it to an initialized DataFrame
    ## which contains the data for ALL replicons in the genome, even if the Count is zero.

    naive_themisto_PCN_estimates_df = pl.read_csv(naive_themisto_PCN_csv_file).select(
        ## select only the columns we need.
        ['AnnotationAccession', 'SeqID', 'SeqType', 'ReadCount', 'replicon_length']
    ).rename({"ReadCount": "InitialReadCount"}) ## rename ReadCount to InitialReadCount

    ## Initialize a dataframe containing SeqID, SeqType, replicon_length metadata for all
    ## replicons in this genome, 
    my_naive_PCN_df = initialize_GenomeDataFrame(themisto_ID_to_seq_metadata_dict)

    ## Add a column for AnnotationAccession.
    annotation_accession_col = [genome_ID for i in range(my_naive_PCN_df.shape[0])]
    AnnotationAccession_series = pl.Series("AnnotationAccession", annotation_accession_col)
    my_naive_PCN_df = my_naive_PCN_df.with_columns([AnnotationAccession_series])
    
    ## and now merge with the naive themisto PCN estimates.
    my_naive_PCN_df = my_naive_PCN_df.join(
        naive_themisto_PCN_estimates_df, on="SeqID", how="left", coalesce=True).with_columns(
            ## set missing values in the InitialReadCount column to 0.
            pl.col("InitialReadCount").fill_null(strategy="zero")).select(
                ## super annoying, not sure why the AnnotationAccession_right columns
                ## are kept. This select command removes these redundancies.
                ["ThemistoID", "AnnotationAccession", "SeqID", "SeqType",
                 "replicon_length", "InitialReadCount"])
    
    ## make a dictionary mapping reads to Themisto replicon IDs.
    multiread_mapping_dict = parse_read_alignments(multiread_alignment_dir)

    ## initialize the data structures for PIRA.
    MatchMatrix, PIRAGenomeDataFrame = initializePIRA(
        multiread_mapping_dict, themisto_ID_to_seq_metadata_dict, my_naive_PCN_df)
    
    ## now run PIRA for this genome.
    PIRA_PCN_estimate_vector = run_PIRA(MatchMatrix, PIRAGenomeDataFrame)
    print(f"PIRA PCN estimate vector for genome {genome_ID} is: {PIRA_PCN_estimate_vector}")
    print("*****************************************************************************************")

    ## now add the PIRA estimates as a column to the PIRAGenomeDataFrame.
    ## First convert the NumPy array to a Polars Series
    PIRA_PCN_estimate_series = pl.Series("PIRA_CopyNumberEstimate", PIRA_PCN_estimate_vector)
    ## Then add the Polars Series of PIRA estimates  to the DataFrame with initial data
    PIRA_PCN_estimates_DataFrame = PIRAGenomeDataFrame.with_columns([PIRA_PCN_estimate_series])

    ## arrange the columns of all_PIRA_estimates_DataFrame in a nice fashion.
    PIRA_PCN_estimates_DataFrame = PIRA_PCN_estimates_DataFrame.select(
        pl.col("AnnotationAccession", "SeqID", "SeqType",
               "ThemistoID", "replicon_length", "InitialReadCount", "AdditionalReadCount", "ReadCount",
               "SequencingCoverage", "LongestRepliconCoverage", "InitialCopyNumberEstimate", "PIRA_CopyNumberEstimate"))
    ## now save PIRA_PCN_estimates_DataFrame to disk.
    PIRA_PCN_estimates_DataFrame.write_csv(PIRA_PCN_csv_file)
    print()
    print(PIRA_PCN_estimates_DataFrame)
    return


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run pseuPIRA for estimating plasmid copy number from microbial genome sequencing data.",
        usage="python pseuPIRA.py -r reference.gbff.gz reads1.fastq [reads2.fastq...]"
    )

    ## Required option for reference file
    parser.add_argument(
        '-r', '--reference',
        required=True,
        help="File containing reference sequence in gzipped GenBank format. (REQUIRED)"
    )

    ## Optional arguments
    parser.add_argument(
        '-j', '--num-processors',
        type=int,
        default=1,
        help="Number of processors to use in multithreaded steps (DEFAULT=1)"
    )
    parser.add_argument(
        '-o', '--output',
        default='../results/test-run',
        help="Path to pseuPIRA output (DEFAULT=../results/test-run)"
    )
    parser.add_argument('-q',
                        '--quick',
                        action='store_true',
                        help="quick mode runs pseudoalignment and skips PIRA"
    )

    ## Positional arguments (input FASTQ files)
    parser.add_argument(
        'reads',
        nargs='+',
        help="FASTQ read files"
    )

    return parser.parse_args()
################################################################################

def main():
    args = parse_args()

    ## Make the output directory if it doesn't exist.
    if not exists(args.output):
        os.mkdir(args.output)

    ## this is used as a standard prefix and ID for various files.
    genome_id = basename(args.reference).split(".gbff.gz")[0].strip("_genomic")
        
    replicon_ref_dir = os.path.join(args.output, "themisto_replicon_references/")
    ## Stage 1: Make FASTA input files for Themisto.
    ## Write out separate fasta files for each replicon in each genome, in a directory for each genome.
    ## Then, write out a text file that contains the paths to the FASTA files of the genomes, one file per line.
    ## See documentation here: https://github.com/algbio/themisto.
    print('making fasta input files...')
    make_replicon_fasta_refs(genome_id, args.reference, replicon_ref_dir)

    ## Stage 2: Build the themisto index.
    index_dir = os.path.join(args.output, "themisto_index/")
    print('building themisto index...')
    build_genome_index(genome_id, replicon_ref_dir, index_dir)

    ## Stage 3: Pseudoalign reads against the index.
    print('pseudoaligning reads...')
    pseudoalignment_dir = os.path.join(args.output, "themisto_pseudoalignments/")
    pseudoalign_reads(genome_id, index_dir, args.reads, pseudoalignment_dir)

    ## Stage 4: generate a CSV file summarizing the themisto pseudoalignment read counts
    print('summarizing pseudoalignment read counts...')
    ## but first map the themisto replicon ID numbers to a (SeqID, SeqType, SeqLength) tuple.
    ## this metadata is needed for calculating PCN estimates.
    themisto_ID_to_seq_metadata_dict = map_themisto_IDs_to_replicon_metadata(genome_id, replicon_ref_dir)
    themisto_results_csvfile_path = os.path.join(args.output, "themisto-replicon-read-counts.csv")
    summarize_pseudoalignment_results(genome_id, themisto_ID_to_seq_metadata_dict, replicon_ref_dir, pseudoalignment_dir, themisto_results_csvfile_path)

    ## Stage 5: estimate plasmid copy numbers using the themisto read counts.
    ## Naive PCN calculation, ignoring multireplicon reads.
    naive_PCN_csv_file = os.path.join(args.output, "naive-PCN-estimates.csv")
    ## make a polars dataframe of the replicon metadata.
    replicon_metadata_df = make_replicon_metadata_dataframe(genome_id, themisto_ID_to_seq_metadata_dict)
    naive_PCN_estimation(replicon_metadata_df, themisto_results_csvfile_path, naive_PCN_csv_file)
    if args.quick: ## then we don't run PIRA in this case.
        quit() 
    
    #############################################################################
    ## Use Probabilistic Iterative Read Assignment (PIRA) to improve PCN estimates.
    #############################################################################

    ## Stage 6: filter fastq reads for multireads.
    multiread_data_dir = os.path.join(args.output, "filtered_multireads") ## directory for filtered multireads.
    filter_fastq_files_for_multireads(multiread_data_dir, pseudoalignment_dir, args.reads)            
    
    ## Stage 7: make FASTA reference genomes with Themisto Replicon IDs
    ## for multiread mapping with minimap2.
    make_fasta_reference_genomes_for_minimap2(genome_id, replicon_ref_dir)
    
    ## Stage 8: for each genome, align multireads to the replicons with minimap2.
    ## directory for multiread alignments constructed with minimap2.
    multiread_alignment_dir = os.path.join(args.output, "multiread_alignments/")
    align_multireads_with_minimap2(genome_id, replicon_ref_dir, multiread_data_dir, multiread_alignment_dir)    
    ## Stage 9: Run PIRA.
    ## parse minimap2 results to form the match matrix and refine the initial PCN guesses.
    ## This file contains PIRA estimates for the genomes that have multireads called by themisto.
    PIRA_PCN_csv_file = os.path.join(args.output, "PIRA-PCN-estimates.csv")
    run_PIRA_on_the_genome(genome_id, multiread_alignment_dir, replicon_ref_dir, naive_PCN_csv_file, PIRA_PCN_csv_file)
        
    return


if __name__ == "__main__":
    main()

