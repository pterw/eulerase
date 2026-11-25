# Changelog

## 0.1.0

The first release. Given sequences and a k, Eulerase builds the de Bruijn graph, finds a minimum set
of Eulerian trails that uses every edge once, and assembles contigs from the graph; it is not a
production assembler. It runs from the command line, which prints the trails and writes figure, CSV,
GFA, FASTA and GenBank files on request, and from a window that opens when `eulerase` is given no
arguments.
