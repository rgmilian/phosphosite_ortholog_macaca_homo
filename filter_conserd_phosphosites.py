import pandas as pd

# Load the results
results_df = pd.read_csv('human_orthologs_mapped.csv')

# Calculate success statistics
total_sites = len(results_df)
orthologs_found = len(results_df[~results_df['Human_Accession'].isin(['Not found', 'Error'])])
success_rate = (orthologs_found / total_sites) * 100

# Calculate amino acid conservation
conserved_aa = 0
for _, row in results_df.iterrows():
    if row['Human_Accession'] not in ['Not found', 'Error']:
        # Extract AA from positions
        macaca_aa = row['Macaca_Position'][0]  # First character (e.g., 'S' from 'S1361')
        if row['Human_Position'] not in ['N/A', 'Alignment failed', 'Position out of range', 'Sequence not available']:
            human_aa = row['Human_Position'][0]  # First character from human position
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
    if row['Human_Accession'] not in ['Not found', 'Error']:
        macaca_aa = row['Macaca_Position'][0]
        if row['Human_Position'] not in ['N/A', 'Alignment failed', 'Position out of range', 'Sequence not available']:
            human_aa = row['Human_Position'][0]
            if macaca_aa == human_aa:
                filtered_results.append(row)

if filtered_results:
    filtered_df = pd.DataFrame(filtered_results)
    filtered_df.to_csv('human_orthologs_mapped_conserved_only.csv', index=False)
    print(f"\nFiltered dataset saved: {len(filtered_df)} conserved phosphosites")
    print(f"File: human_orthologs_mapped_conserved_only.csv")
