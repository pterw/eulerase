# Eulerase

Eulerase is a desktop and command-line tool for learning and exploring de Bruijn graphs. Given
sequences and a k, it builds the graph (nodes are (k-1)-mers, edges are k-mers), finds a minimum set
of Eulerian trails that uses every edge once, and assembles contigs from the graph, with every
simplification visible and chosen by the user. It is for students, teachers and bench biologists
working with short sequences. It is not a production assembler.

This README is minimal: the package is being built in steps, so it describes only what exists. At
this step the repository holds the package layout, the licence and the gates; the tool itself is not
there yet. The README is completed with the tool's features, limits and examples when they ship.

## What it claims

Three claims are kept apart, in the code, the labels and this file:

1. A trail is correct for the graph: it uses each edge once and spells a string whose k-mers are
   those edges.
2. A contig is a valid path in the assembly graph built from the k-mers given.
3. A sequence exists in an organism.

Eulerase makes the first two and never the third. A contig is not a validated assembly.

## What it does not do

Eulerase shows what is in the sequences it is given. It never reads genotype files, never looks up
variant identifiers, and never interprets a result for a person.

## Install and check

Python 3.13 or later is required.

```
python -m pip install -e .[dev]
python -m pytest -q -p no:cacheprovider
python -m ruff check .
python tools/check_fresh_export.py
```

The third command runs the tests again on a copy of the tracked files only, with nothing installed
and nothing untracked in reach, so a test that passes only because of the local setup fails there.
Stage a new file before running it.

## Licence

MIT. See `LICENSE`.
