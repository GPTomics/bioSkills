'''SeqRecord slicing: what survives, what is silently dropped, and the GFF off-by-one'''
# Reference: biopython 1.83+ | Verify API if version differs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

record = SeqRecord(
    Seq('ATGCGATCGATCGATCGTAACCGGTTAACCGGTTAA'),
    id='chr_demo', name='demo', description='annotated demo record',
    annotations={'organism': 'Escherichia coli', 'taxonomy': ['Bacteria', 'Proteobacteria'], 'molecule_type': 'DNA'},
    dbxrefs=['BioProject:PRJNA000000'])
record.features = [
    SeqFeature(SimpleLocation(0, 36), type='source'),       # spans whole record, straddles any inner cut
    SeqFeature(SimpleLocation(4, 18), type='gene'),         # fully inside [4:24] -> kept, location recalculated
    SeqFeature(SimpleLocation(20, 30), type='CDS')]         # straddles the [4:24] cut -> dropped whole

print('=== Original record ===')
print(f'id={record.id}  len={len(record)}  organism={record.annotations.get("organism")}')
print(f'dbxrefs={record.dbxrefs}  features={[f.type for f in record.features]}')

sub = record[4:24]
print('\n=== After record[4:24] ===')
print(f'id preserved: {sub.id}')
print(f'molecule_type preserved: {sub.annotations.get("molecule_type")}')
print(f'annotations dict dropped (organism now {sub.annotations.get("organism")!r})')
print(f'dbxrefs dropped: {sub.dbxrefs}')
print(f'kept features: {[f.type for f in sub.features]}  (source + CDS straddled the cut, gone)')
gene = sub.features[0]
print(f'gene location recalculated to slice-relative coords: {int(gene.location.start)}-{int(gene.location.end)}')

print('\n=== Carrying annotations across deliberately ===')
sub.annotations = record.annotations.copy()
print(f'organism restored: {sub.annotations.get("organism")}')

print('\n=== Per-letter quality auto-slices to match ===')
read = SeqRecord(Seq('ACGTACGTAC'), id='read1', letter_annotations={'phred_quality': list(range(30, 40))})
window = read[2:6]
print(f'read qualities : {read.letter_annotations["phred_quality"]}')
print(f'window[2:6]    : {window.letter_annotations["phred_quality"]}  (sliced for free)')

print('\n=== GFF 1-based off-by-one ===')
genome = Seq('AAAATGCGATCGTAAGGGG')
gff_start, gff_end = 4, 15  # ATG..TAA ORF, 1-based inclusive, as written on a GFF line
wrong = genome[gff_start:gff_end]       # forgot to subtract 1 from start -> shifted left
right = genome[gff_start - 1:gff_end]   # subtract 1 from START only
print(f'naive seq[4:15]      : {wrong}  (off by one, misses the start codon)')
print(f'correct seq[4-1:15]  : {right}  (clean ATG..TAA ORF)')
