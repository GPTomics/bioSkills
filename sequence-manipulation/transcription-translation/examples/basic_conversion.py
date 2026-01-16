'''Basic transcription and translation examples'''
from Bio.Seq import Seq

# Transcription: DNA to RNA
print('=== Transcription ===')
dna = Seq('ATGCGATCGATCG')
rna = dna.transcribe()
print(f'DNA: {dna}')
print(f'RNA: {rna}')

# Back transcription: RNA to DNA
print('\n=== Back Transcription ===')
back_to_dna = rna.back_transcribe()
print(f'RNA: {rna}')
print(f'DNA: {back_to_dna}')

# Translation: DNA to Protein
print('\n=== Translation ===')
coding_dna = Seq('ATGTTTGGTCATTAA')
protein = coding_dna.translate()
print(f'DNA: {coding_dna}')
print(f'Protein: {protein}')

# Translation stopping at stop codon
print('\n=== Translation to Stop ===')
protein_clean = coding_dna.translate(to_stop=True)
print(f'Protein (to_stop=True): {protein_clean}')

# Full pipeline
print('\n=== Full Pipeline ===')
dna = Seq('ATGTTTGGTCATTAA')
rna = dna.transcribe()
protein = rna.translate(to_stop=True)
print(f'DNA:     {dna}')
print(f'RNA:     {rna}')
print(f'Protein: {protein}')
