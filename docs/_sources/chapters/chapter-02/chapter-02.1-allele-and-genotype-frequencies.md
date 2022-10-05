# Allele And Genotype Frequencies

```{margin}
A *locus* (plural: *loci*) is a specific spot in the genome. The term allele was coined by Edith Rebecca Saunders and William Bateson in 1902 in their paper “The facts of heredity in the light of Mendel’s discovery” .
```

In this chapter we will work through how the basics of Mendelian genetics play out at the population level in sexually reproducing organisms.

Loci and alleles are the basic currency of population genetics--and indeed of genetics. A locus may be an entire gene, or a single nucleotide base pair such as A-T. At each locus, there may be multiple genetic variants segregating in the population---these different genetic variants are known as *alleles*. If all individuals in the population carry the same allele, we say that the locus is *monomorphic*; at this locus there is no genetic variability in the population. If there are multiple alleles in the population at a locus, we say that this locus is *polymorphic* (this is sometimes referred to as a segregating site).

:::{margin}
```{figure} ../../illustration_images/alleles_genotypes/Drosophila_mel/Drosophila_mel_mutants.jpg
---
name: figure-2.1
align: left
---
\- *Drosophila melanogaster* holds a special place in the history of genetics and population genetics. From Morgan's fly room discovering the principals of genetics to Dobzhansky's early work on natural genetic variation. <span style="font-size: smaller;">Contributions to the genetics of *Drosophila melanogaster* (1919). Morgan T.H., Bridges C.B., Sturtevant A. H. Image from the [Biodiversity Heritage Library](https://www.biodiversitylibrary.org/page/805594#page/147/mode/1up). Contributed by MBLWHOI Library. Not in copyright.</span>
:::

{numref}`table-2.1` shows a small stretch of orthologous sequence for the ADH locus from samples from *Drosophila melanogaster*, *D. simulans*, and *D. yakuba*. *D. melanogaster* and *D. simulans* are sister species and *D. yakuba* is a close outgroup to the two. Each column represents a single haplotype from an individual (the individuals are diploid but were inbred so they're homozygous for their haplotype). Only sites that differ among individuals of the three species are shown. Site 834 is an example of a polymorphism; some *D. simulans* individuals carry a C allele while others have a T. *Fixed differences* are sites that differ between the species but are monomorphic within the species. Site 781 is an example of a fixed difference between *D. melanogaster* and the other two species.

We can also annotate the alleles and loci in various ways. For example, position 781 is a non-synonymous fixed difference. We call the less common allele at a polymorphism the *minor allele* and the common allele the *major allele*, e.g. at site 1068 the T allele is the minor allele in *D. melanogaster*. We call the more evolutionarily recent of the two alleles the *derived allele* and the older of the two the *ancestral allele*. We infer that the T allele at site 1068 is the derived allele because the C is found in both other species, suggesting that the T allele arose via a C &rarr; T mutation.

:::{list-table} - Variable sites in exons 2 and 3 of the ADH gene in *Drosophila* {cite:p}`mcdonald:91`. The first column (pos.) gives the position in the gene; exon 2 begins at position 778 and we've truncated the dataset at site 1175. The second column gives the consensus nucleotide (con.), i.e. the most common base at that position; individuals with nucleotides that match the consensus are marked with a dash.  The first columns of sequence (a-l) are from *D. melanogaster*; the next columns (a-f) give sequences from *D. simulans*, and the final set of columns (a-l ) from *D. yakuba*. The last column shows whether the difference is a non-synonymous (N) or synonymous (S) change. 
:name: table-2.1

* -
:::

<div style="overflow-x: scroll;">

:::{table}
| pos.   | con.   | a   | b   | c   | d   | e   | f   | g   | h   | i   | j   | k   | l   |    | a   | b   | c   | d   | e   | f   |    | a   | b   | c   | d   | e   | f   | g   | h   | i   | j   | k   | l   | NS/S   |
| :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
| 781   | G   | T   | T   | T   | T   | T   | T   | T   | T   | T   | T   | T   | T   |    | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | NS   |
| 789   | T   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    | C   | C   | C   | C   | C   | C   | C   | C   | C   | C   | C   | C   | S   |
| 808   | A   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    | G   | G   | G   | G   | G   | G   | G   | G   | G   | G   | G   | G   | NS   |
| 816   | G   | T   | T   | T   | T   | -   | -   | -   | -   | -   | -   | -   | T   |    | T   | T   | T   | T   | T   | T   |    | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | S   |
| 834   | T   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | C   | C   | -   | -   | -   | C   |    | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | S   |
| 859   | C   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    | G   | G   | G   | G   | G   | G   | G   | G   | G   | G   | G   | G   | NS   |
| 867   | C   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    | G   | G   | G   | G   | G   | A   | G   | G   | G   | G   | G   | G   | S   |
| 870   | C   | T   | T   | T   | T   | T   | T   | T   | T   | T   | T   | T   | T   |    | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | S   |
| 950   | G   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | A   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | S   |
| 974   | G   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | T   | -   | T   | T   | T   | T   |    | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | S   |
| 983   | T   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    | C   | C   | C   | C   | C   | C   | C   | C   | C   | C   | C   | C   | S   |
| 1019   | C   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | A   | -   | -   | -   | -   | -   | -   | -   | S   |
| 1031   | C   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   | -   | -   | A   | -   | -   | -   | S   |
| 1034   | T   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    |    | C   | C   | C   | C   | C   | -   | -   | C   | -   | C   | C   | S   |
| 1043   | C   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | A   | -   | -   | -   | -   | -   | -   | -   | S   |
| 1068   | C   | T   | T   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | S   |
| 1089   | C   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | A   | A   | A   | A   | A   | A   |    | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | NS   |
| 1101   | G   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    | A   | A   | A   | A   | A   | A   | A   | A   | A   | A   | A   | A   | NS   |
| 1127   | T   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    | C   | C   | C   | C   | C   | C   | C   | C   | C   | C   | C   | C   | S   |
| 1131   | C   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | T   | -   | -   | -   | -   | -   | -   | -   | S   |
| 1160   | T   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | -   |    | -   | -   | -   | -   | -   | -   |    | C   | C   | C   | C   | C   | C   | C   | C   | C   | C   | C   | C   | S   |
:::

</div>


:::{admonition} Question 1
:name: question-2.1
**A)** How many segregating sites does the sample from *D. simulans* have in the ADH gene?

**B)** How many fixed differences are there between *D. melanogaster* and *D. yakuba*?
:::