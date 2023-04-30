from Bio import SeqIO
from Bio.SeqUtils import six_frame_translations
from Bio import PDB
import matplotlib.pyplot as plt
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
        with open(Path("rnk_seq.fasta"), "w") as f:
            f.write(">rna" + str(rna_seq))
        with open(Path("protein_seq.fasta"), "w") as f:
            f.write(">protein\n" +str(protein_seq))

    def get_codons(self) -> None:
        total_codons = len(self.genome) // 3
        codon_counts = {}
        for i in range(0, len(self.genome) - 2, 3):
            codon = str(self.genome[i:i + 3])
            if codon in codon_counts:
                codon_counts[codon] += 1
            else:
                codon_counts[codon] = 1
        codon_freqs = {codon: count / total_codons for codon, count in codon_counts.items()}
        with open(Path("codons.txt"), "w") as f:
            f.writelines(
                [f"{key}: {value}\n" for key, value in sorted(codon_freqs.items(), key=lambda a: a[1], reverse=True)])

    def get_orf(self) -> None:
        print(six_frame_translations(self.genome))

    def print_3d_covid19(self) -> None:
        parser = PDB.PDBParser()

        # загружаем структуру белка из файла covid.pdb
        structure = parser.get_structure('covid', 'covid.pdb')

        # выбираем первую модель структуры (обычно белки содержат только одну модель)
        model = structure[0]

        # создаем объекты для хранения координат атомов и их типов
        coords = []
        types = []

        for chain in model:
            for residue in chain:
                for atom in residue:
                    coords.append(atom.get_coord())
                    types.append(atom.get_name())

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for i in range(len(coords)):
            x, y, z = coords[i]
            if types[i] == 'CA':
                ax.scatter(x, y, z, c='blue', marker='o')
            elif types[i] == 'O':
                ax.scatter(x, y, z, c='red', marker='o')
            else:
                ax.scatter(x, y, z, c='gray', marker='.')

        # настраиваем параметры графика
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('COVID-19 spike protein')

        plt.show()


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
    print("^-^ ^_^ ^.^ ^o^ Getting codons!...")
    genome_analyzer.get_codons()
    print("O_o!! Look! Codons are in codons.txt file! ^-^")
    print("==*******************************==")
    print()
    print()
    print("ORF: ")
    genome_analyzer.get_orf()
    print()
    print("3D covid13!!!!")
    genome_analyzer.print_3d_covid19()
