#import collections.Counter
import matplotlib.pyplot as plt

#1.Read the data
#2.Remove the newline characters from the data
#3.strip the 'N' at the start and end of the DNA sequence 
#4.return a string as the format
def read_file_sequence(file_path):
    sequence = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip() #remove the newline characters from each line
                if line:
                    sequence.append(line)
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None
    
    return "".join(sequence).strip('N')

# analyze_sequence call and pass sequeces' character emerging times(counts) to visualize_coutns to
# draw bar chart
def visualize_counts(title, counts_dict, total_bases):

    # sorted() to make sure the final base displayed in the x-axis in ascending form
    bases_to_plot = sorted([b for b in counts_dict.keys() if b != 'N'])
    
    frequencies = [(counts_dict[base] / total_bases) * 100 for base in bases_to_plot]
    
    # background size, w = 20, h = 6
    plt.figure(figsize=(10, 6))
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    plt.bar(bases_to_plot, frequencies, color=colors[:len(bases_to_plot)])
    
    plt.title(f'{title} Nucleotide proportion chart', fontsize=14)
    plt.xlabel('Nucleotide Symbol', fontsize=12)
    plt.ylabel('Frequency (%)', fontsize=12)
    
    # mark the specific value at the top of each bar
    for i, freq in enumerate(frequencies):
        plt.text(i, freq + 0.5, f'{freq:.2f}%', ha='center', fontsize=10) 

    plt.ylim(0, max(frequencies) * 1.2)
    plt.grid(axis='y', linestyle='--')
    

    filename = f'{title.replace(" ", "_")}_distribution.png'
    plt.savefig(filename)
    plt.close()
    print(f"Graph saved as: {filename}")


def analyze_sequence(sequence, label):
    if not sequence:
        print(f"Warning: {label} sequence is empty.")
        return
    #counts extract character in the DNA sequence and caculate corresponding occurence counts in dict form
    #counts = collections.Counter(sequence)    
    counts = {}
    
    for base in sequence:
        counts[base] = counts.get(base, 0) + 1
    
    total_len = len(sequence)
    
    print(f"\n=== {label} quantitative outcome ===")
    print(f"Total Length: {total_len}")
    
    for base in sorted(counts.keys()):
        count = counts[base]
        print(f"  {base}: {count} ({(count/total_len)*100:.2f}%)")
    
    visualize_counts(label, counts, total_len)


if __name__ == "__main__":
    file_paths = {
        'Gene A (DNA)': 'D:/bioinformatics assignment/gene_a.fa',
        'RNA A (mRNA)': 'D:/bioinformatics assignment/rna_a.fa'
    }
    
    
    for label, path in file_paths.items():
        
        seq_data = read_file_sequence(path)
    
        if seq_data:
            analyze_sequence(seq_data, label)
            