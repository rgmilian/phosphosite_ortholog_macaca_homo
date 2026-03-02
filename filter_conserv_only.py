import pandas as pd

# Load the results
results_df = pd.read_csv('human_orthologs_mapped.csv')

# Calculate success statistics
total_sites = len(results_df)
orthologs_found = len(results_df[~results_df['Human_UniProtID'].isin(['Not found', 'Error'])])
success_rate = (orthologs_found / total_sites) * 100

# Calculate amino acid conservation
conserved_aa = 0
for _, row in results_df.iterrows():
    if row['Human_UniProtID'] not in ['Not found', 'Error']:
        # Extract AA from positions
        macaca_aa = row['Macaca_Phosphosite'][0]  # First character (e.g., 'S' from 'S217')
        if row['Human_Phosphosite'] not in ['N/A', 'Alignment failed', 'Position out of range', 
                                              'Sequence not available', 'Macaca position not verified'] \
           and not row['Human_Phosphosite'].startswith('Position mismatch'):
            human_aa = row['Human_Phosphosite'][0]  # First character from human position
            if macaca_aa == human_aa:
                conserved_aa += 1

conservation_rate = (conserved_aa / orthologs_found) * 100 if orthologs_found > 0 else 0

print(f"\n=== SUCCESS STATISTICS ===")
print(f"Total macaque phosphosites: {total_sites}")
print(f"Human orthologs found: {orthologs_found} ({success_rate:.1f}%)")
print(f"Sites with conserved amino acid: {conserved_aa} ({conservation_rate:.1f}%)")
print(f"Sites with alignment issues: {orthologs_found - conserved_aa}")
print(f"\n=== FOR MANUSCRIPT ===")
print(f"This approach successfully mapped {success_rate:.1f}% of macaque "
      f"phosphosites to human orthologs, with {conservation_rate:.1f}% showing amino acid "
      f"conservation at the mapped position.")

# Optional: Create filtered version with conserved sites only
filtered_results = []
for _, row in results_df.iterrows():
    if row['Human_UniProtID'] not in ['Not found', 'Error']:
        macaca_aa = row['Macaca_Phosphosite'][0]
        if row['Human_Phosphosite'] not in ['N/A', 'Alignment failed', 'Position out of range', 
                                              'Sequence not available', 'Macaca position not verified'] \
           and not row['Human_Phosphosite'].startswith('Position mismatch'):
            human_aa = row['Human_Phosphosite'][0]
            if macaca_aa == human_aa:
                filtered_results.append(row)

if filtered_results:
    filtered_df = pd.DataFrame(filtered_results)
    filtered_df.to_csv('human_orthologs_mapped_conserved_only.csv', index=False)
    print(f"\nFiltered dataset saved: {len(filtered_df)} conserved phosphosites")
    print(f"File: human_orthologs_mapped_conserved_only.csv")
else:
    print(f"\nNo conserved phosphosites found for filtering.")

# Additional statistics breakdown
print(f"\n=== DETAILED BREAKDOWN ===")
error_types = results_df['Human_Phosphosite'].value_counts()
print("\nError/Status types:")
for status, count in error_types.items():
    if status not in ['N/A'] and (status.startswith('Position') or status.startswith('Alignment') 
                                   or status.startswith('Sequence') or status.startswith('Macaca')):
        print(f"  {status}: {count}")

# Show gene-level statistics
print(f"\nUnique Macaca genes processed: {results_df['Macaca_Gene'].nunique()}")
print(f"Unique Human genes mapped: {results_df[results_df['Human_Gene'] != 'Not found']['Human_Gene'].nunique()}")
