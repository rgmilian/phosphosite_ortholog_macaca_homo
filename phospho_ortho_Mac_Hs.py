import pandas as pd
import requests
from time import sleep
from Bio import Align
import re
import os

# Load your data
df = pd.read_csv('limma_aged_vs_young.csv')

# Drop the unnamed column if it exists
if 'Unnamed: 0' in df.columns:
    df = df.drop('Unnamed: 0', axis=1)
if '' in df.columns:
    df = df.drop('', axis=1)

def parse_protein_phosphosite(phosphosite_string):
    """
    Parse Protein_Phosphosite column
    e.g., 'AHNAK_F6S2X3_S217' -> gene='AHNAK', uniprot='F6S2X3', aa='S', position=217
    """
    parts = phosphosite_string.split('_')
    if len(parts) >= 3:
        gene_name = parts[0]
        uniprot_id = parts[1]
        phospho_part = parts[2]
        
        # Extract amino acid and position (e.g., 'S217' -> 'S', 217)
        match = re.match(r'([STY])(\d+)', phospho_part)
        if match:
            aa_type = match.group(1)
            position = int(match.group(2))
            return gene_name, uniprot_id, aa_type, position
    return None, None, None, None

# Check if there's a temp file to resume from
temp_file = 'human_orthologs_mapped_temp.csv'
if os.path.exists(temp_file):
    print(f"Found temp file. Resuming from previous run...")
    existing_results = pd.read_csv(temp_file)
    processed_indices = set(existing_results.index.tolist())
    results = existing_results.to_dict('records')
    print(f"Resuming from row {len(results) + 1}/{len(df)}")
else:
    results = []
    processed_indices = set()

def make_request_with_retry(url, max_retries=5, initial_delay=1):
    """
    Make HTTP request with exponential backoff retry
    """
    delay = initial_delay
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                return response
            elif response.status_code == 429:  # Rate limit
                print(f"Rate limited. Waiting {delay * 2} seconds...")
                sleep(delay * 2)
                delay *= 2
            else:
                print(f"HTTP {response.status_code} for {url}")
                return None
        except requests.exceptions.ConnectionError:
            print(f"Connection error (attempt {attempt + 1}/{max_retries}). Retrying in {delay} seconds...")
            sleep(delay)
            delay *= 2
        except requests.exceptions.Timeout:
            print(f"Timeout error (attempt {attempt + 1}/{max_retries}). Retrying in {delay} seconds...")
            sleep(delay)
            delay *= 2
        except Exception as e:
            print(f"Unexpected error: {e}")
            sleep(delay)
            delay *= 2
    
    print(f"Failed after {max_retries} attempts")
    return None

def get_human_ortholog_and_gene(macaca_uniprot_id):
    """
    Get human ortholog UniProt ID and gene symbol from Macaca mulatta UniProt ID
    """
    try:
        # Query UniProt API for protein information
        url = f"https://rest.uniprot.org/uniprotkb/{macaca_uniprot_id}.json"
        response = make_request_with_retry(url)
        
        if response and response.status_code == 200:
            data = response.json()
            
            # Get gene name from Macaca protein
            gene_name = None
            if 'genes' in data and len(data['genes']) > 0:
                gene_name = data['genes'][0].get('geneName', {}).get('value')
            
            if gene_name:
                return search_human_by_gene(gene_name)
        
        sleep(0.5)  # Rate limiting
        return None, None
        
    except Exception as e:
        print(f"Error processing {macaca_uniprot_id}: {e}")
        return None, None

def search_human_by_gene(gene_name):
    """
    Search for human protein by gene name, return UniProt ID and gene symbol
    """
    try:
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_name}+AND+organism_id:9606+AND+reviewed:true&format=json&size=1"
        response = make_request_with_retry(url)
        
        if response and response.status_code == 200:
            data = response.json()
            if 'results' in data and len(data['results']) > 0:
                result = data['results'][0]
                uniprot_id = result['primaryAccession']
                human_gene = None
                if 'genes' in result and len(result['genes']) > 0:
                    human_gene = result['genes'][0].get('geneName', {}).get('value')
                return uniprot_id, human_gene
        
        sleep(0.5)
        return None, None
        
    except Exception as e:
        print(f"Error searching gene {gene_name}: {e}")
        return None, None

def get_protein_sequence(uniprot_id):
    """
    Get full protein sequence from UniProt
    """
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        response = make_request_with_retry(url)
        
        if response and response.status_code == 200:
            fasta = response.text
            sequence = ''.join(fasta.split('\n')[1:])
            return sequence
        
        sleep(0.5)
        return None
        
    except Exception as e:
        print(f"Error getting sequence for {uniprot_id}: {e}")
        return None

def map_position_with_alignment(macaca_seq, macaca_pos_num, human_seq):
    """
    Map Macaca position to human position using global alignment with Bio.Align
    """
    try:
        # Create PairwiseAligner
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        
        # Perform alignment
        alignments = aligner.align(macaca_seq, human_seq)
        
        if not alignments:
            return None
        
        # Use the best alignment (first one)
        best_alignment = alignments[0]
        
        # Map position from Macaca to human
        macaca_index = 0  # Current position in original Macaca sequence
        
        for i in range(len(best_alignment[0])):
            # Track position in original sequence (skip gaps)
            if best_alignment[0][i] != '-':
                macaca_index += 1
            
            # When we reach the target position in Macaca
            if macaca_index == macaca_pos_num:
                # Count position in human sequence (skip gaps)
                human_index = 0
                for j in range(i + 1):
                    if best_alignment[1][j] != '-':
                        human_index += 1
                
                return human_index
        
        return None
        
    except Exception as e:
        print(f"Error in alignment: {e}")
        return None

# Process each row
total = len(df)

for idx, row in df.iterrows():
    # Skip if already processed
    if idx in processed_indices:
        continue
    
    protein_phosphosite = row['Protein_Phosphosite']
    peptide_sequence = row['Sequence']
    
    # Parse the Protein_Phosphosite column
    macaca_gene, macaca_acc, aa_type, pos_num = parse_protein_phosphosite(protein_phosphosite)
    
    print(f"\nProcessing {len(results) + 1}/{total}: {protein_phosphosite}...")
    
    if not macaca_acc or not aa_type or not pos_num:
        results.append({
            'Macaca_Gene': macaca_gene,
            'Macaca_UniProtID': macaca_acc,
            'Macaca_Phosphosite': f"{aa_type}{pos_num}" if aa_type and pos_num else 'N/A',
            'Sequence': peptide_sequence,
            'Human_Gene': 'Error',
            'Human_UniProtID': 'Error',
            'Human_Phosphosite': 'Invalid format'
        })
        temp_df = pd.DataFrame(results)
        temp_df.to_csv(temp_file, index=False)
        continue
    
    # Get human ortholog and gene symbol
    human_acc, human_gene = get_human_ortholog_and_gene(macaca_acc)
    
    if not human_acc:
        results.append({
            'Macaca_Gene': macaca_gene,
            'Macaca_UniProtID': macaca_acc,
            'Macaca_Phosphosite': f"{aa_type}{pos_num}",
            'Sequence': peptide_sequence,
            'Human_Gene': 'Not found',
            'Human_UniProtID': 'Not found',
            'Human_Phosphosite': 'N/A'
        })
        temp_df = pd.DataFrame(results)
        temp_df.to_csv(temp_file, index=False)
        continue
    
    # Get full sequences
    macaca_full_seq = get_protein_sequence(macaca_acc)
    human_full_seq = get_protein_sequence(human_acc)
    
    if not macaca_full_seq or not human_full_seq:
        results.append({
            'Macaca_Gene': macaca_gene,
            'Macaca_UniProtID': macaca_acc,
            'Macaca_Phosphosite': f"{aa_type}{pos_num}",
            'Sequence': peptide_sequence,
            'Human_Gene': human_gene if human_gene else 'N/A',
            'Human_UniProtID': human_acc,
            'Human_Phosphosite': 'Sequence not available'
        })
        temp_df = pd.DataFrame(results)
        temp_df.to_csv(temp_file, index=False)
        continue
    
    # The pos_num from the Protein_Phosphosite is the absolute position in the full protein
    macaca_abs_pos = pos_num
    
    # Verify the position is correct in the full sequence
    position_verified = False
    if macaca_abs_pos <= len(macaca_full_seq):
        actual_aa = macaca_full_seq[macaca_abs_pos - 1]
        if actual_aa == aa_type:
            print(f"  ✓ Verified: {aa_type} at position {macaca_abs_pos} in Macaca sequence")
            position_verified = True
        else:
            print(f"  ⚠ Warning: Expected {aa_type} at position {macaca_abs_pos}, found {actual_aa}")
            # Position doesn't match - flag this
            results.append({
                'Macaca_Gene': macaca_gene,
                'Macaca_UniProtID': macaca_acc,
                'Macaca_Phosphosite': f"{aa_type}{pos_num}",
                'Sequence': peptide_sequence,
                'Human_Gene': human_gene if human_gene else 'N/A',
                'Human_UniProtID': human_acc,
                'Human_Phosphosite': f'Position mismatch: expected {aa_type}, found {actual_aa}'
            })
            temp_df = pd.DataFrame(results)
            temp_df.to_csv(temp_file, index=False)
            continue
    else:
        print(f"  ⚠ Position {macaca_abs_pos} out of range (protein length: {len(macaca_full_seq)})")
        results.append({
            'Macaca_Gene': macaca_gene,
            'Macaca_UniProtID': macaca_acc,
            'Macaca_Phosphosite': f"{aa_type}{pos_num}",
            'Sequence': peptide_sequence,
            'Human_Gene': human_gene if human_gene else 'N/A',
            'Human_UniProtID': human_acc,
            'Human_Phosphosite': 'Position out of range'
        })
        temp_df = pd.DataFrame(results)
        temp_df.to_csv(temp_file, index=False)
        continue
    
    # Use alignment to map to human sequence
    if position_verified:
        human_pos_num = map_position_with_alignment(
            macaca_seq=macaca_full_seq,
            macaca_pos_num=macaca_abs_pos,
            human_seq=human_full_seq
        )
        
        if human_pos_num:
            # Verify the amino acid at this position
            if human_pos_num <= len(human_full_seq):
                human_aa = human_full_seq[human_pos_num - 1]
                human_phosphosite = f"{human_aa}{human_pos_num}"
                print(f"  → Mapped to Human: {human_aa}{human_pos_num}")
            else:
                human_phosphosite = "Position out of range"
        else:
            human_phosphosite = "Alignment failed"
    else:
        human_phosphosite = "Macaca position not verified"
    
    results.append({
        'Macaca_Gene': macaca_gene,
        'Macaca_UniProtID': macaca_acc,
        'Macaca_Phosphosite': f"{aa_type}{pos_num}",
        'Sequence': peptide_sequence,
        'Human_Gene': human_gene if human_gene else 'N/A',
        'Human_UniProtID': human_acc,
        'Human_Phosphosite': human_phosphosite
    })
    
    # Save progress after EVERY row (safer for connection issues)
    temp_df = pd.DataFrame(results)
    temp_df.to_csv(temp_file, index=False)
    
    # Progress update every 10 rows
    if (len(results)) % 10 == 0:
        print(f"Progress saved: {len(results)}/{total}")

# Save final results
results_df = pd.DataFrame(results)
results_df.to_csv('human_orthologs_mapped.csv', index=False)

# Remove temp file after successful completion
if os.path.exists(temp_file):
    os.remove(temp_file)

print(f"\n✓ Processing complete! Results saved to 'human_orthologs_mapped.csv'")
print(f"Total processed: {len(results)}")
print(f"Human orthologs found: {sum(1 for r in results if r['Human_UniProtID'] not in ['Not found', 'Error'])}")
