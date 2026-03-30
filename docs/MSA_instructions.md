# How to Contruct an MSA

## Description

This git provides a walk through for how to build an MSA starting from compiling a dataset to assessing the completed MSA. 

**Overview of MSA Building**

* Collect data either using a PFAM Family or a PSIBlast search of the non-redundant protein
database
* Clean data by removing duplicate sequence, non-natural constructs (if using PSIBlast),
sequence fragments/sequences that are too short, and those that are too long as determined
by plotting a distribution of sequence lengths.
  * For large alignments (> 20,000 sequences), use FAMSA to construct an MSA.
  *  MAFFT can be used for smaller alignments (< 20,000 sequences) or after using MMSeq
to cluster and sub-sample the data. A typical clustering identity parameter is 0.8 sequence
similarity.
  * For very small datasets (in the 1000s of sequences), other, well-reviewed alignment
algorithms, like Promals3D, can be used.
* Check distribution and length of gaps within the MSA. Gaps should be evenly distributed
and not in overly long runs. One can imagine a segregated dataset (ie. with two distinct
homologies) separating within the MSA. This could create long gaps clustering at the start
alignment for one set of sequences where it fails to align with the other set, and at the end of
the alignment for that second set, where it is not aligning with the first.
* MSAs can be processed using the SCA pre-processing parameters
* Processed MSAs can be used to build a Potts Model using a Stochastic Boltzmann Machine.
The resulting Jij coupling matrix can be compared to the known contact graph for the protein
of interest. Alignment between Jij and the contact map suggests a high-quality MSA.

## Data Collection 

### PFAM 
The most straight-forward method of data collection is finding the associated PFAM family for your sequence or protein family of interest. If the PFAM family is a sufficient representation for your protein family, you can download the PFAM generated alignment as a Stockholm file 
```
alignment (on the left) > choose full > download 
```
It is worth noting that a PFAM family may only contain a segment of a larger protein or one domain of a multi-domain protein. These cases take some careful thought as to how to best collect your data of interest (the intersection between two PFAM families, looking a one PFAM but obtaining full sequences from Uniprot etc.). Any PFAM data can be realigned by removing gaps from the provided alignment and realigned with one of the following MSA Building Tools.  

### PSIBlast
Homologous sequences can be compiled thought an iterative search of the Non-Redundant Protein Database (see Additional Notes for download of the nr database). PSIBlast uses a positional-scoring matrix to search through the database for sequences similar to your query sequence. PSIBlast can be downloaded by obtaining the ncbi package for your computer type and following install instructions: 
```
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```
The following command is used to run a psiblast search 
```
psiblast 
-query <query_sequence_file> 
-db <path/to/nr_database> 
-num_iterations 6 
-evalue 0.0001 
-out <output_file>  
-save_each_pssm  #optional save at each iteration 
-outfmt "format_of_choice, see psiblast options" 
-max_target_seqs 50000 
-num_treads 8
```
### Data Cleaning
Before assembling an MSA, you should remove any sequences that might throw off your alignment such as artificial constructs, sequence fragments, or sequences that are too long. 
Artificial constructs can be an issue when scraping the nr database, to remove use a variation for the following code: 
```python
excluded_keywords = [
    'synthetic construct',
    'vector',
    'synthetic',
    'artificial',
    'cloning',
    'binary',
    'plasmid',
    'transformation',
    'recombinant',
    'engineered',
    'chimeric',
    'fusion'
]

filtered_file = <filename>

kept = []
removed = []

with open(<initial_file>) as f:
    for line in SeqIO.parse(f, "fasta"):
        header = line.description.lower()
        if any(kw in header for kw in excluded_keywords):
            removed.append(line.id)
        else:
            kept.append(line)

SeqIO.write(kept, filtered_file, "fasta")

print(f"Number of Artificial Sequences: {len(removed)}")
print(f"Number of Remaining Sequences: {len(kept)}")
```
To remove sequences based off of length, plot a histogram of sequence length and trim the dataset accordingly. It can be helpful to create an additional histogram after trimming your data.  
```python
#create histogram 
file = <filename>
lengths = [] 

for line in SeqIO.parse(file, "fasta"):
    lengths.append(len(line))
        
plt.hist(lengths, bins=50)
print(len(lengths))
plt.title("Length of Sequences before MSA and Trimming")
plt.xlabel("Length")
plt.ylabel("Number of Sequences")
```
```python
#remove sequences based on length 
kept_len = []
removed_len = []
p = <lower_limit>
n = <upper_limit> 
with open(file) as f:
    for line in SeqIO.parse(f, "fasta"):
        seq = line.seq
        if len(seq) < p or len(seq) > n:
            removed_len.append(line.id)
        else:
            kept_len.append(line)

SeqIO.write(kept_len, <output_file>, "fasta")

print(f"Number of Removed Sequences: {len(removed_len)}")
print(f"Number of Remaining Sequences: {len(kept_len)}")
```

## MSA Building Tools 

### MAFFT 
MAFFT downloading instructions can be found here 
```
https://mafft.cbrc.jp/alignment/software/macosx.html
``` 

To run MAFFT, use the following line
To run a FAMSA from terminal, use the following line 
```
mafft [arguments] --auto input.txt > output.txt
```
The addition of --auto allows mafft to pick with alignment algorithm it thinks is best based off of your data size and shape, but it not necessary and can be replaced with a specific option if desired. 

### FAMSA 
FAMSA, or pyFAMSA, installation instructions can be found here 
```
https://github.com/refresh-bio/FAMSA
https://pyfamsa.readthedocs.io/en/stable/ 
```

To run a FAMSA from terminal, use the following line 
```
famsa [options] <input_file> <output_file>
```
To run a pyFAMSA in python or Jupyter notebook, use a variation of the following block 
```python
for line in SeqIO.parse(<input_file>, "fasta"):
    seq_id = line.id.encode("utf-8")
    seq_str = str(line.seq).encode("utf-8")
    seqs.append(Sequence(seq_id, seq_str))

aligner = Aligner() 
alignment = aligner.align(seqs)

with open(<output_file>, "w") as f:
    for seq in alignment:
        f.write(f">{seq.id.decode('utf-8')}\n")
        f.write(f"{seq.sequence.decode('utf-8')}\n")
```

## MSA Assessment 

### Gapping Distribution 
To determine the location and length of gaps across your MSA, plot "runs" of gaps as a heatmap over starting position. Ideally, no runs should be too long or clustered at the beginning and end of your alignment. 

```python
#count the runs of gaps
gap_runs = []

for line in SeqIO.parse(<file>, "fasta"):
    seq_str = str(line.seq)
    gap_start = None

    gap_count = 0
    for i, char in enumerate(seq_str):
        if char == "-":
            gap_count += 1
            if gap_start == None:
                gap_start = i
        elif gap_count > 0:
            gap_runs.append((gap_start, gap_count))
            gap_count = 0
            gap_start = None
    if gap_count > 0: 
        gap_runs.append((gap_start, gap_count))

#plot gapping distribution 

pos_bin_size = 100  

run_bins = [0, 1, 5, 10, 20, 50, 100, 250, 500, 1000, 10000] #bin size/labels can be changes depending on MSA 
run_labels = [
    "1",
    "2-5",
    "6-10",
    "11–20",
    "21–50",
    "51–100",
    "101–250",
    "251–500",
    "501–1000",
    "1001+"
]

df = pd.DataFrame(gap_runs, columns=["pos", "run_length"])
df["pos_bin"] = (df["pos"] // pos_bin_size) * pos_bin_size

df["run_bin"] = pd.cut(
    df["run_length"],
    bins=run_bins,
    labels=run_labels,
    right=True,
    include_lowest=True
)

heatmap = (df.groupby(["run_bin", "pos_bin"]).size().unstack(fill_value=0))

plt.figure(figsize=(18, 6))
plt.imshow(heatmap.values, aspect="auto", origin="lower", cmap="viridis")
plt.colorbar(label="Frequency")
plt.yticks(ticks=np.arange(len(heatmap.index)),labels=heatmap.index)
ticks=np.arange(len(heatmap.columns))
labels = [col if i % 2 == 0 else '' for i, col in enumerate(heatmap.columns)]
plt.xticks(ticks=ticks, labels=labels, rotation=90)
plt.xlabel(f"Gap Run Start Position")
plt.ylabel("Gap run length")
plt.title("Gap Run Start Position vs Run Length")

plt.tight_layout()
plt.show()
```

### Potts Model Contact Map
Run a Potts Model on your MSA and compare the predicted contacts to the experimentally determined contact map for your sequence or family of interest. 

To assemble a contact map (for residues with alpha-carbons < 8 A apart) from a PDB file: 
```python
from Bio.PDB import *
parser = PDBParser()
pdb_file = "path/to/pdb/file"
structure = parser.get_structure("", pdb_file)

threshold = 8.0
residues = [res for res in lucf_structure[0]["A"]
                if res.has_id("CA")]
original_seq = "".join(
    seq1(residue.resname) 
    for residue in lucf_structure[0]["A"].get_residues() 
    if residue.id[0] == " ")

coords = np.array([res["CA"].get_vector().get_array() for res in residues])
    
diff = coords[:, None] - coords[None, :]
dist_matrix = np.sqrt((diff ** 2).sum(-1))

rows, cols = dist_matrix.shape
i, j = np.indices((rows, cols))

contact_map = (np.abs(i - j) >= 5) & (dist_matrix < threshold)
L = len(original_seq)

plt.figure(figsize=(8, 8))
plt.imshow(contact_map, cmap="Blues", origin="upper")
plt.title("Contact Map for Firefly Luciferase")
plt.xlabel("Residue index")
plt.ylabel("Residue index")
plt.tight_layout()
plt.show()
```

From your Potts Model Results, align and compare predicted contacts to your contact map: 
```python 
potts_file = np.load("path/to/potts/output", allow_pickle=True)

item = potts_file.item()

J = item['J']
h = item['h']

print("J shape:", J.shape)   # (q, q, L, L)
print("h shape:", h.shape)   # (q, L)

J_fnorm = np.sqrt(np.sum(J**2, axis=(2,3))) #calculate the frobenius norm for J values  
print("J_fnorm shape:", J_fnorm.shape)
```
**Further Instructions About Mapping to Come**

 