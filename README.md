COULDNT UPLOAD READS.FASTA BC IT WAS TOO BIG

2 Pseudoalignment implementation
[10 points extra credit simply for attempting this one]
Given some RNA-seq data in FASTA format, find the vector of equivalence class counts.
Your implementation should in the least do the following:
• Given a gene annotation, an index should be produced which can be read in by your
pseudoalignment procedure.
• Your pseudoalignment procedure should take in (1) the previously mentioned index,
(2) RNA-seq data and output equivalence class counts, and (3) a k-mer length.
Please provide equivalence class counts in the following format:
counts number of items in equivalence class isoforms in equivalence class
30                             0                  NA
5                              1             ENST000003679
10379                          2             ENST000003679,ENST000009216
You may implement either the naive version presented in lecture using the hash table or
the skipping (colored de Bruijn graph) version[1]. If you are feeling extra weird, you can also
implement it using a suffix array.
Note: there are some annoyances we didn’t talk about in class. The most immediate is
dealing with the reverse complement. Basically, you have to map both the forward and the
reverse strand. The other is that you have to keep track of new equivalence classes that arise
when taking the intersection in novel ways.
Finally, provide some basic statistics about the size of the equivalence classes and how
much data is mapping to them. This sort of summary can be provided in a figure/plot.
An additional tip: you can take the aligned data and compute equivalence classes directly
from it to check your answer.

a file reads.fasta is missing because it is larger than 25MB
