# Assets

## The Eulerase logo

The mark is the de Bruijn graph of the word EULERASE at k = 2. Its nodes
are the word's letters, and its edges are the word's seven 2-mers: EU, UL,
LE, ER, RA, AS, SE. The left loop is E-U-L-E and the right loop is
E-R-A-S-E.

An Eulerian circuit of this graph, read from E, spells the word. From E
there are two: EULERASE and ERASEULE. Eulerase prints the circuit read
from A, `ASEULERA`, which is the same loop.

To see it, run:

```text
eulerase --kmers "EU UL LE ER RA AS SE"
```

It prints:

```text
Eulerian circuit · 7 edges · 6 nodes
Note: GC and ORFs are not available because the input holds letters other than A, C, G, T and N: E, L, R, S, U.
trail 1 (7 edges): ASEULERA
```

## Files

- `eulerase-logo-light.svg` is for a light page.
- `eulerase-logo-dark.svg` is for a dark page.

The letterforms are Noto Sans, under the SIL Open Font License 1.1,
converted to outlines. The files hold no text elements.
