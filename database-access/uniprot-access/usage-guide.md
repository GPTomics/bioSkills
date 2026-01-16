# UniProt Access Usage Guide

UniProt is the most comprehensive protein sequence and annotation database. Access it programmatically using the REST API.

## API Base URL

```
https://rest.uniprot.org/
```

## Quick Start

```python
import requests

# Fetch protein sequence
response = requests.get('https://rest.uniprot.org/uniprotkb/P04637.fasta')
print(response.text)

# Search for proteins
response = requests.get('https://rest.uniprot.org/uniprotkb/search',
                       params={'query': 'gene:BRCA1 AND organism_id:9606', 'format': 'json'})
results = response.json()
```

## Common Use Cases

### Get All Human Kinases

```python
query = 'organism_id:9606 AND keyword:kinase AND reviewed:true'
params = {'query': query, 'format': 'tsv', 'fields': 'accession,gene_names,protein_name'}
response = requests.get('https://rest.uniprot.org/uniprotkb/search', params=params)
```

### Map Gene Names to UniProt

```python
# Use ID mapping endpoint
response = requests.post('https://rest.uniprot.org/idmapping/run',
                        data={'ids': 'TP53,BRCA1,EGFR', 'from': 'Gene_Name', 'to': 'UniProtKB'})
```

### Download All Proteins for an Organism

```python
# Use stream endpoint for large results
response = requests.get('https://rest.uniprot.org/uniprotkb/stream',
                       params={'query': 'organism_id:9606 AND reviewed:true', 'format': 'fasta'})
```

## Output Formats

| Format | Use Case |
|--------|----------|
| `fasta` | Sequences for alignment/BLAST |
| `json` | Programmatic parsing |
| `tsv` | Tabular analysis, pandas |
| `xml` | Full structured data |
| `txt` | Human-readable flat file |
| `gff` | Feature coordinates |

## Tips

1. **Use reviewed:true** for Swiss-Prot (curated) entries only
2. **Use stream endpoint** for >500 results
3. **Specify fields** in TSV format to reduce data transfer
4. **Check rate limits** - be polite with requests
5. **Cache results** - UniProt updates frequently but not constantly
