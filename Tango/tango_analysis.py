import pandas as pd
import subprocess

# Load your dataframe
df = pd.read_csv('codon_variants.csv', usecols=['barcode','library','aa_substitutions','n_aa_substitutions'])

# Filter for 'LibB' variants with at least one amino acid substitution
df_filtered = df[(df['library'] == 'LibB') & (df['n_aa_substitutions'] > 0)]

# Load the WT sequence
with open('WT.fa', 'r') as f:
    wt_sequence = ''.join(f.readlines()[1:]).strip()

# Function to generate variant sequence from amino acid substitutions
def generate_variant_sequence(wt_sequence, aa_substitutions):
    variant_sequence = list(wt_sequence)
    for sub in aa_substitutions.split():
        pos, new_aa = sub[:-1], sub[-1]  
        pos = int(pos[1:]) - 1 
        variant_sequence[pos] = new_aa
    return ''.join(variant_sequence)

# Function to call Tango and get the aggregation score (with error handling)
def get_tango_score(sequence):
    try:
        tango_process = subprocess.Popen(["Tango", "-", "nt=N", "ct=N", "ph=7.2", "te=310.15", "io=0.05"],
                                        stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        output, _ = tango_process.communicate(input=sequence)
        score_line = output.split('\n')[0] 
        agg_score = float(score_line.split()[-1]) 
        return agg_score
    except Exception as e:
        print(f"Error running Tango for sequence: {sequence}")
        print(e)
        return None 

# Generate sequences, calculate Tango scores, and store results
results = []
for index, row in df_filtered.iterrows():
    variant_sequence = generate_variant_sequence(wt_sequence, row['aa_substitutions'])
    agg_score = get_tango_score(variant_sequence)
    results.append({'barcode': row['barcode'], 'agg_score': agg_score})

# Create DataFrame from results
results_df = pd.DataFrame(results)

# Filter out rows with NaN scores
results_df = results_df.dropna(subset=['agg_score'])

# Save results to CSV
results_df.to_csv("tango_results.csv", index=False)

# Display confirmation message
print("Tango results have been saved to tango_results.csv")


