import pandas as pd
import subprocess
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

# Check Tango installation
if shutil.which("Tango") is None:
    raise FileNotFoundError("Tango executable not found in PATH. Please make sure Tango is installed and accessible.")

# Load your dataframe
df = pd.read_csv('codon_variants.csv', usecols=['barcode', 'library', 'aa_substitutions', 'n_aa_substitutions'])

# Filter for 'LibB' variants with at least one amino acid substitution
df_filtered = df[(df['library'] == 'LibB') & (df['n_aa_substitutions'] > 0)]

# Load the WT sequence
with open('WT.fa', 'r') as f:
    wt_sequence = ''.join(f.readlines()[1:]).strip()

def generate_variant_sequence(wt_sequence, aa_substitutions):
    variant_sequence = list(wt_sequence)
    for sub in aa_substitutions.split():
        pos, new_aa = sub[:-1], sub[-1]
        pos = int(pos[1:]) - 1
        variant_sequence[pos] = new_aa
    return ''.join(variant_sequence)

def process_sequence(args):
    index, row = args
    variant_sequence = generate_variant_sequence(wt_sequence, row['aa_substitutions'])
    scores = get_tango_scores(variant_sequence)
    return scores, row['barcode'], index

def get_tango_scores(sequence):
    cleaned_sequence = sequence.replace('\n', '')
    command = f"Tango P05100 ct=\"N\" nt=\"N\" ph=\"7.2\" te=\"310.15\" io=\"0.05\" seq=\"{cleaned_sequence}\""
    try:
        tango_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
        output, error = tango_process.communicate()
        if tango_process.returncode != 0:
            return None
        scores = {}
        output_items = output.split()
        for i in range(0, len(output_items) - 1, 2):
            key = output_items[i]
            value = output_items[i+1]
            scores[key] = float(value)
        return scores
    except Exception as e:
        return None

# Use ProcessPoolExecutor to run tasks in parallel
results = []
total_sequences = len(df_filtered)
with ProcessPoolExecutor() as executor:
    # Map all tasks to the executor
    futures = {executor.submit(process_sequence, (index, row)): index for index, row in df_filtered.iterrows()}
    completed = 0
    for future in as_completed(futures):
        scores, barcode, index = future.result()
        completed += 1
        print(f"Processing sequence {completed}/{total_sequences} ({(completed/total_sequences)*100:.2f}% completed)")
        if scores is not None:
            results.append({'barcode': barcode, **scores})

# Create DataFrame from results
results_df = pd.DataFrame(results)

# Filter out rows with NaN scores in any column
results_df = results_df.dropna()

# Save results to CSV
results_df.to_csv("tango_results.csv", index=False)

print("Tango results have been saved to tango_results.csv")
