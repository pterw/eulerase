---
name: Feature request
about: Suggest a change, starting from the problem you have
title: "feature: "
---

## The problem

<!-- What are you trying to do, and what stops you? Describe the problem
before any solution. -->

## Before you file

Read the table of kept, parked and dropped features in section 3 of
the design, `docs/superpowers/specs/2026-11-09-eulerase-v0.1-design.md`,
and tick what applies.

- [ ] I read section 3 of the design. My request is not a feature it parks
      or drops.
- [ ] My request is a parked feature (the automatic error model, contig
      binning, or comparison of contigs with a reference). I name the truth
      test it would pass, from section 8.3 of the design (the error model
      and binning each has its own), or from section 10 (the comparison
      design) for the comparison.

A parked feature returns only when it passes its truth tests and the
design document is changed first; being easy to build is not a reason.
Dropped features (molecular weight, the outlier flag and the
principal-component view) are not coming back.

## The boundary

Eulerase shows what is in the sequences it is given. It never reads
genotype files, never looks up variant identifiers, and never interprets
a result for a person (section 1.1 of the design). A request that needs
any of these is out of scope.

## Your idea

<!-- The change you propose, and how you would know it is right: what
independent answer it could be checked against. -->
