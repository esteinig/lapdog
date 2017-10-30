#!/usr/bin/env python

import os
import shutil
import logging
import subprocess

from Bio import SeqIO


class Lapdog:

    """

    Based on process outlined in IMAGE algorithm for gap closing, but hijacked for extension of
    reads in iterative mapping against target genes. Scouty the Chihuahua is responsible for the
    working name of this code.

    https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-4-r41

    """

    def __init__(self, outdir="lapdog", database="mecA.fasta", aligner="bwa_mem",
                 assembler="spades", force=True):

        # Need to convert print statements to logs. Track processing better.
        self.logger = logging.Logger("lapdog.log", level="DEBUG")

        self.outdir = outdir
        self.database = database

        self.aligner = aligner
        self.assembler = assembler

        self.force = force

        if self.force:
            shutil.rmtree(outdir, ignore_errors=True)
        else:
            if os.path.exists(outdir):
                print("Output path exists, use the --force to overwrite.")
            else:
                os.makedirs(outdir)

    def run(self, pe1, pe2, fasta="mecA.fasta", max_rounds="auto", max_length=10000, flash=False):

        """ Run the pipeline on paired-end Illumina reads mapping against a single sequence (.fasta) """

        running = True

        iteration = 1
        previous_length = 0

        # Check for last iteration if not force to continue:

        # Start iterations of target fasta assembly with BWA and SPAdes
        # I think this needs to be encoded in Snakemake for parallel assembly iterations
        # when targeting an entire database of target genes. Currently each iteration takes a few minutes
        # so perhaps we could look into replacing the mapper with something like miniasm or miniasm2.
        while running:

            # Make iteration directories
            alignment_path, assembly_path, iteration_path = self._make_iteration_dir(iteration)

            # Initiate mapper and run BWA
            mapper = Mapper(outdir=alignment_path)
            mapper.run_bwa(pe1=pe1, pe2=pe2, fasta=fasta, bam=True, sort=True)
            # Extract flanking reads unmapped, mapped, both
            mapper.extract_flanking()
            # Convert to FASTQ for assembly
            fastq = mapper.to_fastq()

            # Initiate assembler and run SPAdes (essentially default settings):

            # Error in this process: after a few rounds in the mecA test data
            # there is an odd number of reads in the input FASTQ. There might
            # be a mistake in how I extract the reads above, it should be two reads
            # per ID.
            assembler = Assembler(outdir=assembly_path)
            scaffold = assembler.run_spades(fastq=fastq)

            # Check scaffold length, select largest contig as next iteration's target for mapping:
            scaffold_length, scaffold_file = self._check_scaffold(scaffold, iteration, iteration_path)

            # Assign it to be target fasta explicitly:
            fasta = scaffold_file

            # Terminate while loop conditions:
            if type(max_rounds) == int and iteration == max_rounds:
                print("Reached maximum rounds for mapping and assembly.")
                running = False
            else:
                # Auto mode:
                if scaffold_length >= max_length:
                    print("Assembly is longer than maximum allowed length.")
                    running = False
                if previous_length < scaffold_length < max_length:
                    print("Scaffold is smaller than maximum allowed length but larger than previous iteration.")

            previous_length = scaffold_length
            iteration += 1

        # Need check identity of locally assembled regions  against reference genome, see how well it performed
        # Speed is an issue and definitely need to distribute in Snakemake for multi-fasta input.

        # Compare on how to predict context in complete assemblies (speed, completeness, prediction of controls)
        # In the ST772-MRSA reference genome for example, the contig breaks seem to be associated with repeated
        # insertion sequences (IS). If we do not have to resolve an immediate conflict of IS seqeunces in our
        # local assemblies, perhaps we can extend our assembly region more easily without breaking from repeated
        # IS along the whole genome.

    def _make_iteration_dir(self, iteration):

        iteration_dir = os.path.join(self.outdir, "round_" + str(iteration))

        if not os.path.exists(iteration_dir):
            os.makedirs(iteration_dir)

        return os.path.join(iteration_dir, "alignment"), os.path.join(iteration_dir, "assembly"), iteration_dir

    @staticmethod
    def _check_scaffold(scaffold, iteration, iteration_path):

        records = list(SeqIO.parse(scaffold, format="fasta"))

        count = len(records)

        if count == 1:
            print("There is only one scaffold, proceeding with file:", scaffold)
        elif count <= 1:
            print("There are no scaffolds in", scaffold)
            exit(1)
        else:
            print("There are", count, "scaffolds in", scaffold)
            print("Selecting longest scaffold as reference for next iteration...")

        lengths = [len(record.seq) for record in records]

        # Get longest scaffold, needs fix in case of ties:
        max_length = max(lengths)
        max_index = lengths.index(max_length)

        longest_scaffold = [record for i, record in enumerate(records) if i == max_index]

        scaffold_file = os.path.join(iteration_path, "scaffold_" + str(iteration) + ".fasta")

        SeqIO.write(longest_scaffold, scaffold_file, format="fasta")

        print("Longest scaffold is", max_length, "base pairs.")
        print("New reference file for next iteration is", scaffold_file)

        return max_length, scaffold_file


class Program:

    def __init__(self, outdir, name, version, log, verbose, force):

        self.log = log
        self.name = name
        self.force = force
        self.version = version
        self.verbose = verbose

        self.command = None

        self.outdir = os.path.abspath(outdir)

        if os.path.exists(self.outdir) and not force:
            print("Could not generate output directory, use the --force.")
        else:
            os.makedirs(self.outdir)

    def run(self, command=None):

        """ Run executable of classes that inherit from Program """

        if command is None:
            command = self.command

        # Because shell = True transform to string:
        command = " ".join(command)

        print(command)

        try:
            with open(self.log, "a") as log_file:
                subprocess.call(command, stdout=log_file, stderr=log_file, shell=True)
        except subprocess.CalledProcessError:
            print("Could not run", self.name, self.version, "please see", self.log)


# Minimal wrappers around programs used in the pipeline:
class Assembler(Program):

    def __init__(self, outdir, name="Assembler", version="0.1", spades="spades.py", log="assembler.log",
                 verbose=True, force=False):

        Program.__init__(self, outdir, name, version, log, verbose, force)

        self.spades = spades

        self.scaffold = None

    def run_spades(self, fastq):

        """ Run SPADES on final extracted mapped reads (either mates) - single FASTQ """

        # Made input single-end library even though paired-end due to triplicate reads (need to check);
        # documentation also recommends when using FLASh
        spades_command = [self.spades, "-s", fastq, "-o", self.outdir]

        self.run(spades_command)

        self.scaffold = os.path.join(self.outdir, "contigs.fasta")

        return self.scaffold


class Mapper(Program):

    """ Merges BWA and samtools into mapping and extraction program. """

    def __init__(self, outdir, name="Mapper", version="0.1", bwa="bwa", samtools="samtools", bedtools="bedtools",
                 log="mapper.log", verbose=True, force=False):

        Program.__init__(self, outdir, name, version, log, verbose, force)

        self.samtools = samtools
        self.bedtools = bedtools
        self.bwa = bwa

        self.mapped = None
        self.extracted = None
        self.fastq = None

    def run_bwa(self, pe1, pe2, fasta, output="bwa", bam=True, sort=True):

        """
        Indexes and maps paired-end Illumina reads to sequence (.fasta).

        Optionally exports BAM or sorted BAM (default ON).

        """

        # Output path for SAM / BAM
        output = os.path.join(self.outdir, output)
        # Copy fasta into output directory, then index:
        out_fasta = os.path.join(self.outdir, os.path.basename(fasta))
        shutil.copy(fasta, out_fasta)
        # Index fasta in output directory:
        index_command = [self.bwa, "index", out_fasta]
        # Run indexing:
        self.run(index_command)

        # Make command for mapping:
        map_command = [self.bwa, "mem", out_fasta, pe1, pe2]
        if not bam:
            self.mapped = output + ".sam"
            map_command += [">", self.mapped]
        else:
            map_command += ["|", self.samtools, "view", "-b", "-"]

            if not sort:
                self.mapped = output + ".bam"
                map_command += [">", self.mapped]
            else:
                self.mapped = output + ".sorted.bam"
                map_command += ["|", self.samtools, "sort", "-", ">", self.mapped]
        # Run mapping:
        self.run(map_command)

        return self.mapped

    def extract_flanking(self, bam_file=None):

        """ Extract and merge flanking reads (mates mapped / unmapped) and mapped reads (both mapped) from BAM """

        if bam_file is None:
            bam_file = self.mapped

        file_path, extension = os.path.splitext(bam_file)

        um_name = file_path + "_um" + ".bam"
        mu_name = file_path + "_mu" + ".bam"
        mm_name = file_path + "_mm" + ".bam"

        # Do not output alignments where read is unmapped, do output alignments where mate is unmapped
        # i.e. output all alignments where read is mapped, and mate is unmapped:
        mu = [self.samtools, "view", bam_file, "-b", "-F4", "-f 8", ">", mu_name]
        # Do not output alignments where mate is unmapped, do output alignments where read is unmapped
        # i.e. output all alignments where read is unmapped, and mate is mapped:
        um = [self.samtools, "view", bam_file, "-b", "-F8", "-f 4", ">", um_name]
        # Do not output alignments where read and mate are unmapped
        # i.e. output all alignments where read and mate are mapped
        mm = [self.samtools, "view", bam_file, "-b", "-F12", ">", mm_name]

        self.extracted = file_path + ".mapped.bam"

        merge = [self.samtools, "merge", self.extracted, um_name, mu_name, mm_name]

        for cmd in (um, mu, mm, merge):
            self.run(cmd)

        return self.extracted

    def to_fastq(self, bam_file=None):

        """ Convert BAM to FASTQ """

        if bam_file is None:
            bam_file = self.extracted

        file_path, extension = os.path.splitext(bam_file)

        self.fastq = file_path + ".fastq"

        fastq_command = [self.bedtools, "bamtofastq", "-i", bam_file, "-fq", self.fastq]

        self.run(fastq_command)

        return self.fastq


class FLASH(Program):

    """ FLASh pre-assembles reads into longer contigs really fast. Might be worth seeing if it helps assembly. """

    def __init__(self, outdir, name="FLASh", version="0.1", flash="flash", log="flash.log", verbose=True, force=False):

        Program.__init__(self, outdir, name, version, log, verbose, force)

        self.flash = flash

    def run(self, pe1=None, pe2=None, prefix="out", *args):

        self.command = [self.flash, "-d", self.outdir, "-o", prefix, "--compress"] + list(args) + [pe1, pe2]

        self.run()

        return (os.path.join(self.outdir, prefix + ".extendedFrags.fastq.gz"),
                os.path.join(self.outdir, prefix + ".notCombined_1.fastq.gz"),
                os.path.join(self.outdir, prefix + ".notCombined_2.fastq.gz"))
