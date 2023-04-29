from Bio import SeqIO
from pathlib import Path




class GenomeAnalyser:
    def __init__(self, path_to_genome: Path = Path("covid_sequence.fasta")):
        self.genome = SeqIO.read(path_to_genome.as_posix(), "fasta").seq
        self._a_count = 0
        self._g_count = 0
        self._c_count = 0
        self._t_count = 0

    def _calculate_nucleotide_counts(self) -> None:
        self._a_count = 0
        self._g_count = 0
        self._c_count = 0
        self._t_count = 0
        for nucleotide in self.genome:
            if nucleotide == "A":
                self._a_count += 1
            elif nucleotide == "G":
                self._g_count += 1
            elif nucleotide == "C":
                self._c_count += 1
            elif nucleotide == "T":
                self._t_count += 1

    def get_and_print_nucleotide_counts(self) -> None:
        self._calculate_nucleotide_counts()
        print("Частота нуклеотида A: ", self._a_count / len(self.genome))
        print("Частота нуклеотида G: ", self._g_count / len(self.genome))
        print("Частота нуклеотида C: ", self._c_count / len(self.genome))
        print("Частота нуклеотида T: ", self._t_count / len(self.genome))


if __name__ == '__main__':
    genome_analyzer = GenomeAnalyser()
    genome_analyzer.get_and_print_nucleotide_counts()
