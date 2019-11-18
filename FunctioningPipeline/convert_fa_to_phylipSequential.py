# convert .fa to .phylip for partition finder 2

import glob, os
from Bio import AlignIO

FOLDER = 'sampledMSAs'

for f in glob.glob(os.path.join(FOLDER, '*.fa')):
    input_handle = open(f, "rU")
    output_handle = open(f.replace('fa','phy'), "w")

    align = AlignIO.read(input_handle, "fasta")
    AlignIO.write(align, output_handle, "phylip-sequential")

    output_handle.close()
    input_handle.close()