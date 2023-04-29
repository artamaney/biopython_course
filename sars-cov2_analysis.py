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

    def get_nucleotide_counts(self) -> None:
        self._calculate_nucleotide_counts()
        print("Частота нуклеотида A: ", self._a_count / len(self.genome))
        print("Частота нуклеотида G: ", self._g_count / len(self.genome))
        print("Частота нуклеотида C: ", self._c_count / len(self.genome))
        print("Частота нуклеотида T: ", self._t_count / len(self.genome))

    def get_gc_content(self):
        gc_content = 100 * sum([1 for base in self.genome if base in ["G", "C"]]) / len(self.genome)

        print(f"GC content of the genome is {gc_content:.2f}%")

    def get_rnk_and_protein(self) -> None:
        rna_seq = self.genome.transcribe()
        protein_seq = rna_seq.translate()
        with open(Path("rnk_seq.txt"), "w") as f:
            f.write(str(rna_seq))
        with open(Path("protein_seq.txt"), "w") as f:
            f.write(str(protein_seq))


if __name__ == '__main__':
    genome_analyzer = GenomeAnalyser()
    genome_analyzer.get_nucleotide_counts()
    print("==*******************************==")
    genome_analyzer.get_gc_content()
    print("==*******************************==")
    genome_analyzer.get_rnk_and_protein()
    print("Generated files with RNK sequence: rnk_seq.txt!")
    print("//0_0\\\\: Wow!")
    print("Generated files with Protein sequence: rnk_seq.txt!")
