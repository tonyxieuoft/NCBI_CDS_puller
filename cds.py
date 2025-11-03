import os
import time
from os import path

from Bio import Entrez


def string_empty(str):
    return str != ""

def extract_cds_from_accession(file, out_dirpath):

    out_path_name = path.splitext(path.basename(file))[0]

    file_num = 1
    out_path = path.join(out_dirpath, out_path_name + ".fas")
    while path.isfile(out_path):

        out_path = path.join(out_dirpath, out_path_name +
                             " (" + str(file_num) + ").fas")
        file_num += 1

    output_f = open(out_path, "w")

    accessions = filter(string_empty, open(file).read().split())

    for acc in accessions:

        # get taxid
        summary_handle = Entrez.esummary(db='nuccore', id=acc)
        summaries = Entrez.parse(summary_handle)
        taxa = "0"
        for summary in summaries:
            taxa = int(summary["TaxId"])

        # get name from taxid
        tax_fetch = Entrez.efetch(db="taxonomy", id=taxa)
        tax = Entrez.parse(tax_fetch)
        first, last = None, None
        for t in tax:
            first, last = t["ScientificName"].split()[:2]

        # get CDS range
        success, attempts = False, 0
        while not success and attempts < 5:
            try:
                # gets the fasta file in .txt format
                sequence_handle = Entrez.efetch(db='nuccore', id=acc,
                                                rettype='ft',
                                                retmode='text')
                success = True
            except:
                time.sleep(1)
                attempts += 1
                print("attempting to fetch sequence ")

        if attempts == 5:
            raise Exception("Failure to get feature table for accession: " + acc)

        output = str(sequence_handle.read()).split()
        start = None
        end = None
        ge = False
        for i in range(len(output)):
            if (i - 2 >= 0) and output[i] == "CDS" and \
                    output[i-1].isnumeric() and \
                    output[i-2].isnumeric() and \
                    int(output[i-2]) < int(output[i-1]):
                start = output[i-2]
                end = output[i-1]
            if (i - 2 >= 0) and output[i] == "CDS" and \
                    len(output[i-1]) > 0 and \
                    output[i-1][1:].isnumeric() and \
                    output[i-2].isnumeric() and \
                    int(output[i-2]) < int(output[i-1][1:]) and \
                    output[i-1][0] == ">":
                start = output[i-2]
                end = output[i-1][1:]
            if (i -2 >= 0) and output[i] == "CDS" and \
                    len(output[i-2]) > 0 and \
                    output[i-2][1:].isnumeric() and \
                    output[i-1].isnumeric() and \
                    int(output[i-2][1:]) < int(output[i-1]) and \
                    output[i-2][0] == "<":
                start = output[i-2][1:]
                end = output[i-1]

        # get sequence
        if start is None:
            print("Issue: no CDS for accession " + acc)
        else:
            success, attempts = False, 0
            while not success and attempts < 5:
                try:
                    # gets the fasta file in .txt format
                    sequence_handle = Entrez.efetch(db='nuccore', id=acc,
                                                    rettype='fasta',
                                                    retmode='text',
                                                    seq_start = start,
                                                    seq_stop = end)
                    success = True
                except:
                    time.sleep(1)
                    attempts += 1
                    print("attempting to fetch sequence " + acc)

            if attempts == 5:
                raise Exception("Failure to retrieve sequence with accession: " + acc)

            output_f.write(">" + first[0] + last[:3] + "\n" + "\n".join(sequence_handle.read().split("\n")[1:]))
            #output_f.write(sequence_handle.read())
            print("success for accession: " + acc)


if __name__ == "__main__":

    Entrez.email = input("Please enter your email. This is required to access "
                         "NCBI Genbank via the API.")

    in_dir = input("Please enter the complete path to a directory of .txt files "
                               "containing the accessions. The files should only have one accession per "
                               "line.")

    out_dir = input("Please enter the complete path to a "
                    "directory you wish the files to be outputted "
                    "to. It is recommended that you create a new directory "
                    "first, then specify its path here. ")

    while not path.isdir(out_dir):
        accession_dirpath = input("Invalid directory path, please enter again.")

    for file_name in os.listdir(in_dir):
        file_path = os.path.join(in_dir, file_name)
        extract_cds_from_accession(file_path, out_dir)

    print("done")
