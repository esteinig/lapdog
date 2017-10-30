#!/usr/bin/env python

from lapdog.lapdog import Lapdog


def main():

    ld = Lapdog(outdir="woof", force=True)

    ld.run(pe1="DAR_4145_1.fastq.gz", pe2="DAR_4145_2.fastq.gz", fasta="DAR4145.mecA.fasta")


main()



