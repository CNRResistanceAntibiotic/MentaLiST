# -*- coding: utf-8 -*-
import logging


import argparse
import collections
import os
from Bio import SeqIO


def main(param):
    logger = logging.getLogger()

    logging.basicConfig(level=param.loglevel, format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s:%(message)s', datefmt='%I:%M:%S %p')

    # Opening novel alleles:
    logger.info("Opening the novel alleles file ...")
    novel = collections.defaultdict(list)
    for seq_record in SeqIO.parse(param.novel, "fasta"):
        novel[seq_record.id].append(seq_record)

    # create folder if it does not exist:
    if not os.path.isdir(param.output):
        os.makedirs(param.output)
    # Open mlst
    mlst = {}
    logger.info("Opening the MLST schema and adding novel alleles ...")
    for f in os.listdir(param.pathDB):
        f_path = os.path.join(param.pathDB, f)
        logger.debug("Opening file %s ..." % f)
        file_no_ext, ext = os.path.splitext(f_path)
        locus = os.path.basename(f_path).split(".")[0]
        record_list = [seq_record for seq_record in SeqIO.parse(f_path, "fasta")]
        # if there are novel alleles for this locus, add:

        if len(novel[locus]) > 0:
            # find maximum id present, novel alleles gets next;
            id_list = [int(s.id.split("_")[-1]) for s in record_list]
            next_id = max(id_list) + 1
            # append novels:
            for record in novel[locus]:
                record.id += "_%d" % next_id
                record.name = record.description = ""
                next_id += 1
                record_list.append(record)
        # save:
        SeqIO.write(record_list, os.path.join(param.output, os.path.basename(f_path)), "fasta")
    logger.info("Done.")


def run():
    parser = argparse.ArgumentParser(description="Adds novel alleles to an existing MLST scheme.")
    parser.add_argument("-n", "--novel", type=str, help="FASTA with novel alleles.")
    parser.add_argument("-o", "--output", type=str, help="Output folder for new scheme.")
    parser.add_argument("-i", "--id", type=int, default=1000, help="Start numbering new alleles on this value, later will implement from last allele id +1.")
    # parser.add_argument("-t", "--threads", type=int, default=4, help="number of threads")
    parser.add_argument("-db", "--pathDB", type=str, help="MLST Fasta Database Directory")
    parser.add_argument('-ll', '--loglevel', type=str, default="INFO", choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], help='Set the logging level')
    param = parser.parse_args()

    main(param)


if __name__ == '__main__':
    run()