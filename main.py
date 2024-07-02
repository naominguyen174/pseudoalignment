from collections import defaultdict
import matplotlib.pyplot as plt

class Transcript:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

def parse_fasta(filename):
    sequences = []
    try:
        with open(filename, 'r') as f:
            name = None
            sequence = ""
            for line in f:
                if line.startswith(">"):
                    if name is not None:
                        sequences.append(Transcript(name, sequence))
                    name = line.strip()[1:]
                    sequence = ""
                else:
                    sequence += line.strip()
            if name is not None:
                sequences.append(Transcript(name, sequence))
    except FileNotFoundError:
        print(f"Error: The file {filename} was not found.")
    except Exception as e:
        print(f"An error occurred while parsing the file: {e}")
    return sequences

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_sequence = reversed(sequence)
    reverse_complement_seq = []
    for base in reversed_sequence:
        if base in complement:
            reverse_complement_seq.append(complement[base])
        else:
            reverse_complement_seq.append(base)  # handle non-standard bases
    return ''.join(reverse_complement_seq)

def extract_reads_from_fasta(filename):
    reads = []
    try:
        with open(filename, 'r') as f:
            sequence = ""
            for line in f:
                if line.startswith(">"):
                    if sequence:
                        reads.append(sequence)
                        sequence = ""
                else:
                    sequence += line.strip()
            if sequence:
                reads.append(sequence)
    except FileNotFoundError:
        print(f"Error: The file {filename} was not found.")
    except Exception as e:
        print(f"An error occurred while parsing the file: {e}")
    return reads

def build_index(transcripts, k):
    index = defaultdict(list)
    for transcript in transcripts:
        sequence = transcript.sequence
        rev_sequence = reverse_complement(sequence)
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            rev_kmer = rev_sequence[i:i+k]
            index[kmer].append(transcript.name)
            index[rev_kmer].append(transcript.name)  # Index reverse complement k-mer as well
    return index

def pseudoalignment(index, reads, k):
    equivalence_classes = defaultdict(set)
    for read in reads:
        kmers_seen = set()
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            rev_kmer = reverse_complement(kmer)  # Calculate reverse complement k-mer
            if kmer in index and kmer not in kmers_seen:
                for transcript in index[kmer]:
                    equivalence_classes[transcript].add(kmer)
                kmers_seen.add(kmer)
            elif rev_kmer in index and rev_kmer not in kmers_seen:
                for transcript in index[rev_kmer]:
                    equivalence_classes[transcript].add(rev_kmer)
                kmers_seen.add(rev_kmer)
    return equivalence_classes

def format_equivalence_class_counts(equivalence_classes, all_transcript_names):
    counts = defaultdict(list)
    matched_transcripts = set(equivalence_classes.keys())

    for transcript, kmers in equivalence_classes.items():
        count = len(kmers)
        counts[count].append(transcript)

    formatted_counts = []

    # Add formatted counts for matched transcripts
    for count in sorted(counts):
        transcripts = counts[count]
        formatted_counts.append(
            f"count: {count}\n"
            f"number of items in equivalence class: {len(transcripts)}\n"
            f"isoforms in equivalence class: {', '.join(transcripts)}\n"
        )

    # Add entry for missing isoforms
    #it just wont print anything unaligned. It should only show aligned transcripts/isoforms
    unmatched_transcripts = set(all_transcript_names) - matched_transcripts
    if unmatched_transcripts:
        formatted_counts.append(
            f"count: {count}\n"
            f"number of items in equivalence class: 0\n"
            f"isoforms in equivalence class: NA\n"
        )

    return formatted_counts


def main():
    transcripts = parse_fasta("chr11.fasta")
    print(f"Total number of transcripts: {len(transcripts)}")

    reads = extract_reads_from_fasta("reads.fasta")
    print(f"Total number of reads extracted: {len(reads)}")

    k = 25
    index = build_index(transcripts, k)
    print(f"Index size: {len(index)}")

    equivalence_classes = pseudoalignment(index, reads, k)

    num_classes = len(equivalence_classes)
    avg_size = sum(len(kmers) for kmers in equivalence_classes.values()) / num_classes if num_classes > 0 else 0

    all_transcript_names = [transcript.name for transcript in transcripts]
    formatted_counts = format_equivalence_class_counts(equivalence_classes, all_transcript_names)

    print("Total Mapped Reads:", len(reads))
    print("Number of Equivalence Classes:", num_classes)
    print("Average Equivalence Class Size:", avg_size, "\n")

    # Print all equivalence class counts
    print('\n'.join(formatted_counts))

    sizes = []
    for formatted_count in formatted_counts:
        count_line = formatted_count.split('\n')[0]
        size = int(count_line.split(':')[1])
        sizes.append(size)

    if sizes:
        plt.hist(sizes, bins=range(0, max(sizes) + 10, 10), color='skyblue', edgecolor='black', linewidth=1.2)
        plt.xlabel('Size of Equivalence Classes')
        plt.ylabel('Frequency')
        plt.title('Equivalence Class Size Distribution')
        plt.grid(True)
        plt.show()
    else:
        print("No equivalence class sizes found.")

    # Plot histogram for mapped reads per transcript
    mapped_reads_counts = [len(kmers) for kmers in equivalence_classes.values()]
    plt.hist(mapped_reads_counts, bins=range(0, max(mapped_reads_counts) + 10, 10), color='salmon', edgecolor='black', linewidth=1.2)
    plt.xlabel('Number of Mapped Reads')
    plt.ylabel('Frequency')
    plt.title('Mapped Reads per Transcript Distribution')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
