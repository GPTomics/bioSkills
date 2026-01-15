'''Reading GenBank files and extracting features'''
from Bio import SeqIO

print('=== Parsing GenBank File ===')
for record in SeqIO.parse('sample.gb', 'genbank'):
    print(f'ID: {record.id}')
    print(f'Description: {record.description}')
    print(f'Sequence length: {len(record.seq)} bp')
    print(f'Number of features: {len(record.features)}')

    print('\n=== Features ===')
    for feature in record.features:
        print(f'  Type: {feature.type}')
        print(f'  Location: {feature.location}')
        if feature.type == 'CDS':
            product = feature.qualifiers.get('product', ['Unknown'])[0]
            print(f'  Product: {product}')
        print()

    print('=== Annotations ===')
    for key, value in record.annotations.items():
        print(f'  {key}: {value}')
