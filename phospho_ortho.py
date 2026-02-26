import pandas as pd
import requests
from time import sleep
from Bio import Align
import re

# Load your data - don't use index_col since the unnamed column should be ignored
df = pd.read_csv('phosphosites_pos_aa_sequences.csv')

# Drop the unnamed column if it exists
if 'Unnamed: 0' in df.columns:
    df = df.drop('Unnamed: 0', axis=1)

def get_human_ortholog(macaca_uniprot_id):
    """
    Get human ortholog UniProt ID from Macaca mulatta UniProt ID
    """
    try:
        # Query UniProt API for protein information
        url = f"https://rest.uniprot.org/uniprotkb/{macaca_uniprot_id}.json"
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()
            
            # Get gene name from Macaca protein
            gene_name = None
            if 'genes' in data and len(data['genes']) > 0:
                gene_name = data['genes'][0].get('geneName', {}).get('value')
            
            if gene_name:
                return search_human_by_gene(gene_name)
        
        sleep(0.5)  # Rate limiting
        return None
        
    except Exception as e:
        print(f"Error processing {macaca_uniprot_id}: {e}")
        return None

def search_human_by_gene(gene_name):
    """
    Search for human protein by gene name
    """
    try:
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_name}+AND+organism_id:9606+AND+reviewed:true&format=json&size=1"
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()
            if 'results' in data and len(data['results']) > 0:
                return data['results'][0]['primaryAccession']
        
        sleep(0.5)
        return None
        
    except Exception as e:
        print(f"Error searching gene {gene_name}: {e}")
        return None

def get_protein_sequence(uniprot_id):
    """
    Get full protein sequence from UniProt
    """
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        response = requests.get(url)
        
        if response.status_code == 200:
            fasta = response.text
            sequence = ''.join(fasta.split('\n')[1:])
            return sequence
        
        sleep(0.5)
        return None
        
    except Exception as e:
        print(f"Error getting sequence for {uniprot_id}: {e}")
        return None

def extract_aa_position_number(aa_position):
    """
    Extract amino acid type and position number
    e.g., 'S1361' -> ('S', 1361)
    """
    match = re.match(r'([A-Z])(\d+)', aa_position)
    if match:
        aa_type = match.group(1)
        position = int(match.group(2))
        return aa_type, position
    return None, None

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

def find_position_in_full_sequence(peptide_seq, peptide_pos_num, full_seq):
    """
    Find where the peptide maps in the full sequence and calculate absolute position
    """
    try:
        # Try exact match first
        if peptide_seq in full_seq:
            peptide_start = full_seq.find(peptide_seq)
            absolute_position = peptide_start + peptide_pos_num
            return absolute_position
        
        # If exact match fails, use alignment
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        alignments = aligner.align(peptide_seq, full_seq)
        
        if not alignments:
            return None
        
        best_alignment = alignments[0]
        
        # Find where peptide starts in full sequence
        peptide_index = 0
        for i in range(len(best_alignment[0])):
            if best_alignment[0][i] != '-':
                peptide_index += 1
            
            if peptide_index == peptide_pos_num:
                # Count position in full sequence
                full_pos = 0
                for j in range(i + 1):
                    if best_alignment[1][j] != '-':
                        full_pos += 1
                return full_pos + best_alignment.aligned[1][0][0]
        
        return None
        
    except Exception as e:
        print(f"Error finding position: {e}")
        return None

# Process each row
results = []
total = len(df)

for idx, row in df.iterrows():
    macaca_acc = row['Accession']
    macaca_pos = row['AA_position']
    peptide_sequence = row['Sequence']
    
    print(f"Processing {idx + 1}/{total}: {macaca_acc} {macaca_pos}...")
    
    # Extract amino acid and position
    aa_type, pos_num = extract_aa_position_number(macaca_pos)
    
    if not aa_type or not pos_num:
        results.append({
            'Macaca_Accession': macaca_acc,
            'Macaca_Position': macaca_pos,
            'Sequence': peptide_sequence,
            'Human_Accession': 'Error',
            'Human_Position': 'Invalid position format'
        })
        continue
    
    # Get human ortholog
    human_acc = get_human_ortholog(macaca_acc)
    
    if not human_acc:
        results.append({
            'Macaca_Accession': macaca_acc,
            'Macaca_Position': macaca_pos,
            'Sequence': peptide_sequence,
            'Human_Accession': 'Not found',
            'Human_Position': 'N/A'
        })
        continue
    
    # Get full sequences
    macaca_full_seq = get_protein_sequence(macaca_acc)
    human_full_seq = get_protein_sequence(human_acc)
    
    if not macaca_full_seq or not human_full_seq:
        results.append({
            'Macaca_Accession': macaca_acc,
            'Macaca_Position': macaca_pos,
            'Sequence': peptide_sequence,
            'Human_Accession': human_acc,
            'Human_Position': 'Sequence not available'
        })
        continue
    
    # Find absolute position in Macaca full sequence
    macaca_abs_pos = find_position_in_full_sequence(
        peptide_seq=peptide_sequence, 
        peptide_pos_num=peptide_sequence.find(aa_type) + 1 if aa_type in peptide_sequence else pos_num,
        full_seq=macaca_full_seq
    )
    
    if not macaca_abs_pos:
        # Try using the provided position directly
        macaca_abs_pos = pos_num
    
    # Use alignment to map to human sequence
    human_pos_num = map_position_with_alignment(
        macaca_seq=macaca_full_seq,
        macaca_pos_num=macaca_abs_pos,
        human_seq=human_full_seq
    )
    
    if human_pos_num:
        # Verify the amino acid at this position
        if human_pos_num <= len(human_full_seq):
            human_aa = human_full_seq[human_pos_num - 1]
            human_position = f"{human_aa}{human_pos_num}"
        else:
            human_position = "Position out of range"
    else:
        human_position = "Alignment failed"
    
    results.append({
        'Macaca_Accession': macaca_acc,
        'Macaca_Position': macaca_pos,
        'Sequence': peptide_sequence,
        'Human_Accession': human_acc,
        'Human_Position': human_position
    })
    
    # Save intermediate results every 10 rows
    if (idx + 1) % 10 == 0:
        temp_df = pd.DataFrame(results)
        temp_df.to_csv('human_orthologs_mapped_temp.csv', index=False)
        print(f"Saved temporary results ({idx + 1}/{total})")

# Save final results
results_df = pd.DataFrame(results)
results_df.to_csv('human_orthologs_mapped.csv', index=False)
print(f"\nProcessing complete! Results saved to 'human_orthologs_mapped.csv'")
print(f"Total processed: {len(results)}")
print(f"Human orthologs found: {sum(1 for r in results if r['Human_Accession'] not in ['Not found', 'Error'])}")
