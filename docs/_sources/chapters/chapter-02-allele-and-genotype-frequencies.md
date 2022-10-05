# Allele and Genotype Frequencies

```{margin}
A *locus* (plural: *loci*) is a specific spot in the genome. The term allele was coined by Edith Rebecca Saunders and William Bateson in 1902 in their paper “The facts of heredity in the light of Mendel’s discovery” .
```

In this chapter we will work through how the basics of Mendelian genetics play out at the population level in sexually reproducing organisms.

Loci and alleles are the basic currency of population genetics--and indeed of genetics. A locus may be an entire gene, or a single nucleotide base pair such as A-T. At each locus, there may be multiple genetic variants segregating in the population---these different genetic variants are known as *alleles*. If all individuals in the population carry the same allele, we say that the locus is *monomorphic*; at this locus there is no genetic variability in the population. If there are multiple alleles in the population at a locus, we say that this locus is *polymorphic* (this is sometimes referred to as a segregating site).

:::{margin}
```{figure} ../illustration_images/alleles_genotypes/Drosophila_mel/Drosophila_mel_mutants.jpg
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

## Allele frequencies

Allele frequencies are a central unit of population genetics analysis, but from diploid individuals we only get to observe genotype counts. Our first task then is to calculate allele frequencies from genotype counts. Consider a diploid autosomal locus segregating for two alleles ($A_1$ and $A_2$). We'll use these arbitrary labels for our alleles, merely to keep this general. Let $N_{11}$ and $N_{12}$ be the number of $A_1A_1$ homozygotes and $A_1A_2$ heterozygotes, respectively. Moreover, let $N$ be the total number of diploid individuals in the population. We can then define the relative frequencies of $A_1A_1$ and $A_1A_2$ genotypes as $f_{11} = N_{11}/N$ and $f_{12} = N_{12}/N$, respectively. The frequency of allele $A_1$ in the population is then given by

:::{math}
:label: eq-2.1
p = \frac{2 N_{11} + N_{12}}{2N} = f_{11} + \frac{1}{2} f_{12}.
:::

Note that this follows directly from how we count alleles given individuals' genotypes, and holds independently of Hardy--Weinberg proportions and equilibrium (discussed below). The frequency of the alternate allele ($A_2$) is then just $q=1-p$.

### Measures of genetic variability

##### Nucleotide diversity ($\pi$)

One common measure of genetic diversity is the average number of single nucleotide differences between haplotypes chosen at random from a sample. This is called *nucleotide diversity* and is often denoted by $\pi$. For example, we can calculate $\pi$ for our ADH locus from {numref}`table-2.1` above: we have 6 sequences from *D. simulans* (a-f), there's a total of 15 ways of pairing these sequences, and

:::{math}
:label: eq-2.2
\pi=\frac{1}{15} \big( (2 + 1 + 1 + 1 + 0 ) + (3 + 3 + 3 + 2 ) +(0 + 0 + 1) + (0 + 1) + (1)  \big)=1.2\overline{6}
:::

where the first bracketed term gives the pairwise differences between a and b-f, the second bracketed term the differences between b and c-f and so on.

````{margin}
Technically we would need to divide by the total number of possible point mutations that would result in a synonymous change; this is because some mutational changes at a particular nucleotide will result in a non-synonymous or synonymous change depending on the base-pair change.
````
Our $\pi$ measure will depend on the length of sequence it is calculated for. Therefore, $\pi$ is usually normalized by the length of sequence, to be a per site (or per base) measure. For example, our ADH sequence covers 397bp of DNA and so $\pi = \frac{1.2\overline{6}}{397}=0.0032$ per site in *D. simulans* for this region. Note that we could also calculate $\pi$ per synonymous site (or non-synonymous). For synonymous site $\pi$, we would count up number of synonymous differences between our pairs of sequences, and then divide by the total number of sites where a synonymous change could have occurred.

##### Number of segregating sites

Another measure of genetic variability is the total number of sites that are polymorphic (segregating) in our sample. One issue is that the number of segregating sites will grow as we sequence more individuals (unlike $\pi$). Later in the course, we'll talk about how to standardize the
number of segregating sites for the number of individuals sequenced (see **\eqn \eqref{watterson_theta}**).

##### The frequency spectrum

We also often want to compile information about the frequency of alleles across sites.  We call alleles that are found once in a sample *singletons*, alleles that are found twice in a sample *doubletons*, and so on. We count up the number of loci where an allele is found $i$ times out of $n$, e.g. how many singletons are there in the sample, and this is called the *frequency spectrum*. We'll want to do this in some consistent manner, such as calculating the frequency spectrum of the minor allele or the derived allele.

:::{admonition} Question 2
:name: question-2.2
How many minor-allele singletons are there in *D. simulans* in the ADH region? [Defining minor allele just within *D. simulans*.]
:::

##### Levels of genetic variability across species

Two observations have puzzled population geneticists since the inception of molecular population genetics. The first is the relatively high level of genetic variation observed in most obligately sexual species. This first observation, in part, drove the development of the Neutral theory of molecular evolution, the idea that much of this molecular polymorphism may simply reflect a balance between genetic drift and mutation.

The second observation is the relatively narrow range of polymorphism across species with vastly different census sizes. This observation represented a puzzle as the Neutral theory predicts that levels of genetic diversity should scale with population size. Much effort in theoretical and empirical population genetics has been devoted to trying to reconcile models with these various observations. We'll return to discuss these ideas throughout our course.

:::{margin}
```{figure} ../illustration_images/alleles_genotypes/Ciona_intestinalis/21016139168_2a8a57ded3_z.jpg
---
name: figure-2.2
align: left
---
\- Sea Squirt (*Ciona intestinalis*). <span style="font-size: smaller;">Einleitung in die vergleichende gehirnphysi- ologie und Vergleichende psychologie. Loeb, J. 1899. Image from the Biodiversity Heritage Library. Contributed by MBLWHOI Library. No known copyright restrictions.</span>
```
:::

The first observations of molecular genetic diversity within natural populations were made from surveys of allozyme data, but we can revisit these general patterns with modern data. For example, {cite:p}`leffler:12` compiled data on levels of within-population, autosomal nucleotide diversity ($\pi$) for 167 species across 14 phyla from non-coding and synonymous sites ({numref}`figure-2.3`). The species with the lowest levels of $\pi$ in their survey was Lynx, with $\pi = 0.01\%$, i.e. only $\frac{1}{10000}$ bases differed between two sequences. In contrast, some of the highest levels of diversity were found in *Ciona savignyi*, Sea Squirts, where a remarkable $\frac{1}{12}$ bases differ between pairs of sequences. This 800-fold range of diversity seems impressive, but census population sizes have a much larger range.

```{figure} ../Journal_figs/alleles_genotypes/Leffer_riddle/Leffer_riddle_diversity.png
---
name: figure-2.3
align: left
---
\- Levels of autosomal nucleotide diversity for 167 species across 14 phyla. Figure 1 from {cite:p}`leffler:12`, licensed under CC BY 4.0. Points are ranked by their $\pi$, and coloured by their phylum. Note the log-scale.
```

:::{margin}
```{figure} ../illustration_images/alleles_genotypes/Lynx/20731949565_8a065700af_z.jpg
---
name: figure-2.4
align: left
---
\- Eurasian Lynx (*Lynx lynx*). <span style="font-size: smaller;">An introduction to the study of mammals living and extinct. Flower, W.H. and Lydekker, R. 1891. Image from the [Biodiversity Heritage Library](https://www.flickr.com/photos/internetarchivebookimages/20731949565/in/photolist-x5Jzv2-x6QVyp-xir9rH-wYHrQD-wPn1sP-w9PsqY-xDcqri-sMcQoB-trrkVd-x6Nx1H-wPea7N-sM28N9-tJ3zsp-xneVdx-wGJRtQ-xnfHZ8-wPfga7-xCUPrN-x7kXDV-xmAb9E-xm3x4k-xBoSKb-wGTgyB-xBoSbf-wGGvzA-xmzYTJ-oeKJcH-xA1Ffr-xA1Eji-xqWTQZ-xF4Lru-oxJfrH-x7ojSn-xra8zP-wGGibY-xgb21y-xY1jH9-xY1iyf-wGHxTS-wGQEoR-xmtPQh-x8uFKK-xdTkoU-wPggQf-wPfvHN-wPfc27-w9YnGF-wPeauS-wPiVxK-w6aiSi). Contributed by Cornell University Library. No known copyright restrictions.</span>
```
:::

### Hardy--Weinberg proportions

Imagine a population mating at random with respect to genotypes, i.e. no inbreeding, no assortative mating, no population structure, and no sex differences in allele frequencies. The frequency of allele $A_1$ in the population at the time of reproduction is $p$. An $A_1A_1$ genotype is made by reaching out into our population and independently drawing two $A_1$ allele gametes to form a zygote. Therefore, the probability that an individual is an $A_1A_1$ homozygote is $p^2$. This probability is also the expected frequencies of the $A_1A_1$ homozygote in the population. The expected frequency of the three possible genotypes are

<center>
<div style="width: 50%;">

:::{table}
| $f_{11}$ | $f_{12}$ | $f_{22}$ |
| :------: | :------: | :------: |
| $p^2$    | $2pq$    | $q^2$    |
:::

</div>
</center>

i.e. their Hardy-Weinberg expectations {cite:p}`hardy1908mendelian,weinberg1908ber`. Note that we only need to assume random mating with respect to our focal allele in order for these expected frequencies to hold in the zygotes forming the next generation. Evolutionary forces, such as selection, change allele frequencies within generations, but do not change this expectation for new zygotes, as long as $p$ is the frequency of the $A_1$ allele in the population at the time when gametes fuse. We only need the assumptions of no migration, selection, and mutation in order for these Hardy-Weinberg expectations of genotypes to represent a long term equilibrium.

:::{margin}
```{figure} ../illustration_images/alleles_genotypes/Kermode_bear/EIiK2dsWoAAmUOf.jpeg
---
name: figure-2.5
align: left
---

\- Kermode's bear (*Ursus americanus kermodei*). It's possible that this morph is favoured as the salmon these bears eat have a harder time seeing the light morph {cite:p}`klinka2009adaptive`. The adaptive value of tasting like cinnamon is unknown. <span style="font-size: smaller;">Field book of North American mammals; descriptions of every mammal known north of the Rio Grande. Anthony, (1928) H. E. Image from the [Biodiversity Heritage Library](https://www.biodiversitylibrary.org/item/38166#page/115/mode/1up). Contributed by MBLWHOI Library. No known copyright restrictions.</span>
```
:::

:::{admonition} Question 3
:name: question-2.3
On the coastal islands of British Columbia there is a subspecies of black bear (*Ursus americanus kermodei*, Kermode's bear). Many members of this black bear subspecies are white; they're sometimes called spirit bears. These bears aren't hybrids with polar bears, nor are they albinos. They are homozygotes for a recessive change at the MC1R gene. Individuals who are $GG$ at this SNP are white, while $AA$ and $AG$ individuals are black.

Below are the genotype counts for the MC1R polymorphism in a sample of bears from British Columbia's island populations from {cite:t}`RITLAND:01`.

<center>
<div style="width: 50%;">

```{table}
| $AA$ | $AG$ | $GG$ |
| :--: | :--: | :--: |
| 42 | 24 | 21 |
```

</div>
</center>

What are the expected frequencies of the three genotypes under HW?
:::

See {numref}`figure-2.6` for a nice empirical demonstration of Hardy--Weinberg proportions. The mean frequency of each genotype closely matches its HW expectations, and much of the scatter of the dots around the expected line is due to our small sample size (~60 individuals). While HW often seems like a silly model, it often holds remarkably well within populations. This is because individuals don't mate at random, but they do mate at random with respect to their genotype at most of the loci in the genome.

:::{admonition} Question 4
:name: question-2.4
You are investigating a locus with three alleles, A, B, and C, with allele frequencies $p_A$, $p_B$, and $p_C$. What fraction of the population is expected to be homozygotes under Hardy--Weinberg?
:::

:::{margin}
CODIS: Combined DNA Index System
:::

Microsatellites are regions of the genome where individuals vary for the number of copies of some short DNA repeat that they carry. These regions are often highly variable across individuals, making them a suitable way to identify individuals from a DNA sample. This so-called DNA fingerprinting has a range of applications from establishing paternity and identifying human remains to matching individuals to DNA samples from a crime scene. The FBI make use of the CODIS database. The CODIS database contains the genotypes of over 13 million people, most of whom have been convicted of a crime. Most of the profiles record genotypes at 13 microsatellite loci that are tetranucleotide repeats (since 2017, 20 sites have been genotyped).

The allele counts for two loci (D16S539 and TH01) are shown in table {numref}`table-2.2` and {numref}`table-2.3` for a sample of 155 people of European ancestry. You can assume these two loci are on different chromosomes.

:::{list-table} - Data for 155 Europeans at the D16S539 microsatellite from CODIS from {cite:t}`algee:16`. The top row gives the number of tetranucleotide repeats for each allele, the bottom row gives the sample counts.
:name: table-2.2
:header-rows: 0
* - allele name
  - 80
  - 90
  - 100
  - 110
  - 120
  - 121
  - 130
  - 140
  - 150
* - allele count
  - 3
  - 34
  - 13
  - 102
  - 97
  - 1
  - 44
  - 13
  - 3
:::

:::{list-table} - Same as {numref}`table-2.2` but for the TH01 microsatellite.
:name: table-2.3
:header-rows: 0

* - allele name
  - 60
  - 70
  - 80
  - 90
  - 93
  - 100
  - 110
* - allele count
  - 84
  - 42
  - 37
  - 67
  - 77
  - 1
  - 2
:::

:::{admonition} Question 5
:name: question-2.5
You extract a DNA sample from a crime scene. The genotype is 100/80 at the D16S539 locus and 70/93 at TH01.

**A)** You have a suspect in custody. Assuming this suspect is innocent and of European ancestry, what is the probability that their genotype would match this profile by chance (a false-match probability)?

**B)** The FBI uses $\geq$ 13 markers. Why is this higher number necessary to make the match statement convincing evidence in court?

**C)** An early case that triggered debate among forensic geneticists was a crime among the Abenaki, a Native American community in Vermont (see {cite:p}`lewontin:94` for discussion). There was a DNA sample from the crime scene, and the perpetrator was thought likely to be a member of the Abenaki community. Given that allele frequencies vary among populations, why would people be concerned about using data from a non-Abenaki population to compute a false match probability?
:::

```{figure} ../figures/CEU_YRI_separately_HWE.png
:name: figure-2.6
:align: left

\- Demonstrating Hardy--Weinberg proportions using 10,000 SNPs from the HapMap European (CEU)  and African (YRI) populations. Within each of these populations the allele frequency against the frequency of the 3 genotypes; each SNP is represented by 3 different coloured points. The solid lines show the mean genotype frequency. The dashed lines show the predicted genotype frequency from Hardy--Weinberg equilibrium. Code [here](https://github.com/cooplab/popgen-notes/blob/master/Rcode/HWE_exercise/HWE_HAPMAP.R). Blog post on figure [here](http://gcbias.org/2011/10/13/population-genetics-course-resources-Hardy--Weinberg-eq/). 
```

### Assortative mating

One major violation of the assumptions of Hardy Weinberg is non-random mating with respect to the genotype at a locus. One way that individuals can mate non-randomly is if individuals choose to mate based on a phenotype determined by (in part) the genotype at a locus. This non-random mating can be between: 1) individuals with similar phenotype, so called positive assortative mating or 2) individuals with dissimilar phenotypes, negative assortative mating or disassortative mating. Here we'll briefly discuss a couple of real examples of assortative mating to make sure we're all on the same page. We'll encounter other forms of non-random mating, due to inbreeding and population structure, in the next few chapters.

:::{figure} ../Journal_figs/alleles_genotypes/Heliconius_Merrill_assort_mating/trimmed_Heliconius_Merrill_assort_mating.png
:name: figure-2.7
:align: left

\- Wing pattern phenotypes of top, *H. cydno chioneus* (left), *H. melpomene rosina* (right), their nonmimetic first-generation hybrid (center); and bottom, their sympatric comimics *H. sapho sapho* (left) and *H. erato demophoon* (right). Figure and caption modified from {cite:p}`merrill2019genetic`, licensed under CC BY 4.0.
:::

Positive assortative mating on the basis of a phenotype can create an excess of homozygotes. Heliconius butterflies are famous for their mimicry, where poisonous pairs of distantly related species mimic each others' bright colour patterns and so share the benefits of being avoided by visual predators (Müllerian mimics). *H. melpomene rosina* and *H. cydno chioneus* are closely related species that co-occur in central Panama, but mimic different other co-occuring species ({numref}`figure-2.7`). These differences in colouration pattern are due to a few loci with large phenotypic effects. The two species can hybridize and produce viable F1 hybrids. These F1 hybrids are heterozygotes at the colour loci, and their intermediate appearance means that they're poor mimics and so are quickly eaten by predators. However, these heterzygote (F1) hybrids are very rare in nature $<\frac{1}{1000}$, as the parental species show strong positive assortatively mating based on colour pattern, based on genetic differences in mate preference {cite:p}`merrill2019genetic`.

:::{margin}
```{figure} ../illustration_images/alleles_genotypes/White_throated_sparrow/white_throated_sparrow.jpg
:name: figure-2.8
:align: left

\- White-throated sparrows (*Zonotrichia albicollis*) with a white morph (bottom, male) and tan morph (top, female). The difference between the morphs wasn't fully appreciated until the 1960's {cite:p}`lowther1961polymorphism`, previously birders thought the tan morphs were just young or females individuals (so Audubon's male and female labels may well by wrong). There are also a number of behavioural differences, with both sexes of the white-striped morph invest more in territorial defense and the tan-striped morphs more parental care. <span style="font-size: smaller;">From John James Audubon's Birds of America (1827). Image from [Audubon.org](https://www.audubon.org/birds-of-america/white-throated-sparrow), public domain.
```
:::

Disassortative mating, mating of unlike individuals, can lead to an excess of heterozygotes and a deficit of homozygotes. One example of very strong disassortative mating is offered by white-throated sparrows (*Zonotrichia albicollis*). In white-throated sparrows, there is a white-striped and a tan-striped morph, with female and male white-striped morphs have a much brighter white stripe and throat. There is very strong disassortative mating in this system, with 1099 out of 1116 nesting pairs consisting of one tan- and one white-striped morph and only 17 of these nesting pairs being different morphs {cite:p}`tuttle2016divergence`. The difference between these morphs has a simple inheritance pattern, with white being due to a single dominant allele (called 2m) and tan colour from a recessive allele called 2. Thus strong disassortative mating has a strong effect on the genotype frequencies:

<center>
<div style="width: 50%;">

:::{table}
| Tan | White | (Super)White |
| :-: | :---: | :----------: |
| 2/2 | 2/2m  | 2m/2m        |
| 978 | 1011  | 3            |
:::

</div>
</center>

:::{margin}
```{figure} ../illustration_images/alleles_genotypes/Sparassis_crispa/Sparassis_crispa_cauliflower_mushroom.pdfSparassis_crispa_cauliflower_mushroom.pdf
:name: figure-2.9
:align: left

\- Cauliflower mushrooms (*Sparassis crispa*) parasitize tree roots and form these amazing, edible fruiting bodies, which can weigh in at up to 30lb and apparently taste like noodles. In a collection of 18 fruiting bodies from a *Sparassis* population, all individuals were heterozygotes for mating type and 17 different mating types were genetically identified {cite:p}`james2015mushrooms,martin1978decay`. <span style="font-size: smaller;">Atlas champignons comestibles et vénéneux (1891 ). Dufour, L Image from the Biodiversity Heritage Library. Contributed by New York Botanical Garden. Not in copyright.</span>
```
:::

There are almost no 2m homozygotes (so called Super white individuals) despite the 2m allele being common in the population (data from {cite:p}`tuttle2016divergence`, table S1).

Another important example of disassortative mating are mating type systems, which are present in many fungi, algae, and protozoa. Gametes of the same species can only fuse to form a zygote if they differ in mating type. The mating type of gametes is genetically controled by a mating type locus, and so individuals are nearly always heterozygous at this locus. In some groups of organisms, there are just two different alleles, in other clades these loci have tens or hundreds of alleles.

## Allele sharing among related individuals and Identity by Descent

All of the individuals in a population are related to each other by a giant pedigree (family tree). For most pairs of individuals in a population these relationships are very distant (e.g. distant cousins), while some individuals will be more closely related (e.g. sibling/first cousins). All individuals are related to one another by varying levels of relatedness, or *kinship*. Related individuals can share alleles that have both descended from the shared common ancestor. To be shared, these alleles must be inherited through all meioses connecting the two individuals (e.g. surviving the $\frac{1}{2}$ probability of segregation each meiosis). As closer relatives are separated by fewer meioses, closer relatives share more alleles. In {numref}`figure-2.10`, we show the sharing of chromosomal regions between two cousins. As we'll see, many population and quantitative genetic concepts rely on how closely related individuals are, and thus we need some way to quantify the degree of kinship among individuals.

```{figure} ../figures/Cousins_IBD_chromo_cartoon.png
:name: figure-2.10
:align: left

\- First cousins sharing a stretch of chromosome identical by descent. The different grandparental diploid chromosomes are coloured so we can track them and recombinations between them across the generations. Notice that the identity by descent between the cousins persists for a long stretch of chromosome due to the limited number of generations for recombination. The squares represent males and circles females.
```

:::{margin}
{cite:p}`cotterman:40,malecot:48`
:::

We will define two alleles to be identical by descent (IBD) if they are identical due to transmission from a common ancestor in the past few generations. For the moment, we ignore mutation, and we will be more precise about what we mean by 'past few generations' later on. For example, parent and child share exactly one allele identical by descent at a locus, assuming that the two parents of the child are randomly mated individuals from the population. In {numref}`figure-2.16`, I show a pedigree demonstrating
some configurations of IBD.

:::{margin}
```{figure} ../figures/sharing_relatives/IBD_0_1_2.pdf
:name: figure-2.11
:align: left

\- A pair of diploid individ- uals (i and j) sharing 0, 1, or 2 alleles IBD where lines show the sharing of alleles by descent (e.g. from a shared ancestor). Here we’ll focus on IBD of outbred individuals. Dealing with sharing between inbred individuals requires 6 more identity-by-descent r coefficients, which honestly makes my head spin.
```
:::

One summary of how related two individuals (let's call them $i$ and $j$) are is the probability that our pair of individuals share 0, 1, or 2 alleles identical by descent (see {numref}`figure-2.11`). We denote these identity-by-descent probabilities by $r_0$, $r_1$, and $r_2$ respectively. See {numref}`table-2.4` for some examples. We can also interpret these probabilities as genome-wide averages. For example, on average, at a quarter of all their autosomal loci full-sibs share zero alleles identical by descent.

One summary of relatedness that will be important is the probability that two alleles (I & J) picked at random, one from each of the two different individuals $i$ and $j$, are identical by descent ($P(\text{I&J IBD})$). We call this quantity the *coefficient of kinship* of individuals $i$ and $j$, and denote it by $F_{ij}$. It is calculated as


:::{math}
:label: eq-2.3
F_{ij} = P(\text{I&J IBD})
:::

:::{math}
:label: eq-2.4
F_{ij} = P(\text{I&J IBD | i&j 0 IBD})P(\text{i&j 0 IBD}) 

+ P(\text{I&J IBD | i&j 1 IBD})P(\text{i&j 1 IBD}) 

+ P(\text{I&J IBD | i&j 2 IBD})P(\text{i&j 2 IBD})
:::

:::{math}
:label: eq-2.5
F_{ij} = 0 \times r_0 + \frac{1}{4} r_1  + \frac{1}{2} r_2.
:::

In the above step, Equation {eq}`eq-2.4`, we're summing the conditional probability of alleles $I$ & $J$ being IBD over whether our individuals $i$ & $j$ share $0$, $1$, or $2$ alleles IBD, an example of using the Law of Total Probability (see Appendix Equation {eq}`eq-A.12`). We've then, in Equation {eq}`eq-2.5`, used the fact that we can calculate our condition probabilities of I & J being IBD using the rules of Mendelian transmision. Consider the probability $P(\text{I&J IBD} ~|~ \text{i&j  1 IBD})$, i.e. that our pair of alleles ($I$ & $J$) drawn from individuals $i$ and $j$ are IBD given that $i$ and $j$ share one allele IBD, this is a $\frac{1}{4}$ as we need to draw the allele that is IBD from both $i$ and $j$, i.e. drawing both black alleles in the middle panel of {numref}`figure-2.11`, which happens with probability $\frac{1}{2} \times \frac{1}{2}$. The coefficient of kinship will appear multiple times, in both our discussion of inbreeding and in the context of phenotypic resemblance between relatives.

:::{table} - Probability that two individuals ($i$ and $j$) of a given relationship share 0, 1, or 2 alleles identical by descent on the autosomes. *Assuming that our individuals are outbred and that this the only close relationship the pair shares.
:name: table-2.4

| Relationship (i,j)*    | P(i&j 0 IBD)  | P(i&j 1 IBD)  | P(i&j 2 IBD)  | P(I&J IBD)      |
| :-----------------:    | :----------:  | :----------:  | :----------:  | :--------:      |
| Relationship (i,j)*    | $r_{0}$       | $r_{1}$       | $r_{2}$       | $F_{ij}$        |
| parent-child           | 0             | 1             | 0             | $\frac{1}{4}$   |
| full siblings          | $\frac{1}{4}$ | $\frac{1}{2}$ | $\frac{1}{4}$ | $\frac{1}{4}$   |
| Monozygotic twins      | 0             | 0             | 1             | $\frac{1}{2}$   |
| 1<sup>st</sup> cousins | $\frac{3}{4}$ | $\frac{1}{4}$ | 0             | $\frac{1}{16}$  |
:::

:::{admonition} Question 6
:name: question-2.6
What are $r_0$, $r_1$, and $r_2$ for $\frac{1}{2}$ sibs? ($\frac{1}{2}$ sibs share one parent but not the other).
:::

:::{admonition} Question 7
:name: question-2.7
Explain in words why $P(\text{I&J IBD} ~|~ \text{i&j  2 IBD}) = \frac{1}{2}$.
:::

##### Genotypic sharing between pairs of individuals

Our $r$ coefficients are going to have various uses. For example, they allow us to calculate the probability of the genotypes of a pair of relatives. Consider a biallelic locus where allele $A_1$ is at frequency $p$, and two individuals who have IBD allele sharing probabilities $r_0$, $r_1$, $r_2$. What is the overall probability that these two individuals are both homozygous for allele 1? Well that's

:::{math}
:label: eq-2.6
P(\text{both }A_{1}A_{1}) = P(\text{both }A_{1}A_{1}\text{ | 0 alleles IBD})P(\text{0 alleles IBD})

+ P(\text{both }A_{1}A_{1}\text{ | 1 alleles IBD})P(\text{1 alleles IBD})

+ P(\text{both }A_{1}A_{1}\text{ | 2 alleles IBD})P(\text{2 alleles IBD})
:::

Or, in our $r_0$, $r_1$, $r_2$ notation:

:::{math}
:label: eq-2.7
P(\text{both }A_{1}A_{1}) = P(\text{both }A_{1}A_{1}\text{ | 0 alleles IBD})r_{0}

+ P(\text{both }A_{1}A_{1}\text{ | 1 alleles IBD})r_{1}

+ P(\text{both }A_{1}A_{1}\text{ | 2 alleles IBD})r_{2}
:::

If our pair of relatives share $0$ alleles IBD, then the probability that they are both homozygous is $P(\textrm{both } A_1 A_1 | \text{0 alleles IBD}) =p^2 \times p^2$, as all four alleles represent independent draws from the population. If they share $1$ allele IBD, then the shared allele is of type $A_1$ with probability $p$, and then the other non-IBD allele, in both relatives, also needs to be $A_1$ which happens with probability $p^2$, so $P(\textrm{both } A_1 A_1 | \text{1 alleles IBD})=p \times p^2$. Finally, our pair of relatives can share two alleles IBD, in which case $P(\textrm{both } A_1 A_1 | \text{2 alleles IBD}) = p^2$, because if one of our individuals is homozygous for the $A_1$ allele, both individuals will be. Putting this all together our Equation {eq}`eq-2.7` becomes

:::{math}
:label: eq-2.8
P(\text{both }A_{1}A_{1}) = p^{4}r_{0} + p^{3}r_{1} + p^{2}r_{2}
:::

Note that for specific cases we could also calculate this by summing over all the possible genotypes their shared ancestor(s) had; however, that would be much more involved and not as general as the form we have derived here.

We can write out terms like Equation {eq}`eq-2.8` for all of the possible configurations of genotype sharing/non-sharing between a pair of individuals. Based on this we can write down the expected number of polymorphic sites where our individuals are observed to share 0, 1, or 2 alleles.

:::{margin}
```{figure} ../illustration_images/alleles_genotypes/Florida_scrub_jay/8576533889_3a131ffc4c_z.jpg
:name: figure-2.12
:align: left

\- Florida Scrub-Jays (*Aphelocoma coerulescens*). <span style="font-size: smaller;">The birds of America : from drawings made in the United States and their territories. 1880. Audubon J.J. Image from the [Biodiversity Heritage Library](https://www.biodiversitylibrary.org/page/40447048\#page/169/mode/1up). Contributed by Smithsonian Libraries. Licensed under CC BY-2.0.</span>
```
:::

:::{admonition} Question 8 (Trickier question)
:name: question-2.8
The genotype of our suspect in [Question 5](question-2.5) turns out to be 100/80 for D16S539 and 70/80 at TH01. The suspect is not a match to the DNA from the crime scene; however, they could be a sibling.

Calculate the joint probability of observing the genotype from the crime
and our suspect:

- **A)** Assuming that they share no close relationship.
- **B)** Assuming that they are full sibs.
- **C)** Briefly explain your findings.
:::

There's a variety of ways to estimate the relationships among individuals using genetic data. An example of using allele sharing to identify relatives is offered by the work of Nancy Chen [in collaboration with Stepfanie Aguillon, see {cite:p}`chen:16,Aguillon:17`. {cite:t}`chen:16` has collected genotyping data from thousands of Florida Scrub Jays at over ten thousand loci. These Jays live at the Archbold field site, and have been carefully monitored for many decades allowing the pedigree of many of the birds to be known. Using these data, she estimates allele frequencies at each locus. Then by equating the observed number of times that a pair of individuals share $0$, $1$, or $2$ alleles to the theoretical expectation, she estimates the probability of $r_0$, $r_1$, and $r_2$ for each pair of birds. A plot of these are shown in {numref}`figure-2.13`, showing how well the estimates match those known from the pedigree.

:::{figure} ../figures/FSJ_IBD.jpg
:name: figure-2.13
:align: left

\- Estimated coefficient of kinship from Florida Scrub Jays. Each point is a pair of individuals, plotted by their estimated IBD ($r_1$ and $r_2$) from their genetic data. The points are coloured by their known pedigree relationships. Note that most pairs have low kinship, and no recent genealogical relationship, and so appear as black points in the lower left corner. Thanks to Nancy Chen for supplying the data. Code [here](https://github.com/cooplab/popgen-notes/blob/master/Rcode/FSJ_IBD/FSJ_IBD_plotting.R).
:::

##### Sharing of genomic blocks among relatives.

:::{figure} ../figures/sharing_relatives/First_cousin_overlap.pdf
:name: figure-2.14
:align: left

\- A simulation of sharing between first cousins. The regions of your grandmother's 22 autosomes that you inherited are coloured red, those that your cousins inherited are coloured blue. In the third panel we show the overlapping genomic regions in purple, these regions will be IBD in you and your cousin. If you are full first cousins, you will also have shared genomic regions from your shared grandfather, not shown here. Details about how we made these simulations [here](https://gcbias.org/2013/12/02/how-many-genomic-blocks-do-you-share-with-a-cousin/).
:::

We can more directly see the sharing of the genome among close relatives using high-density SNP genotyping arrays. In {numref}`figure-2.14` we show a simulation of first cousins' genomic sharing from their shared grandmother. Colored purple are regions where they have matching genomic material, due to having inherited it IBD from their shared grandmother.

First cousins will share at least one allele of your genotype at all of the polymorphic loci in these purple regions. There's a range of methods to detect such sharing. One way is to look for unusually long stretches of the genome where two individuals are never homozygous for different alleles. By identifying pairs of individuals who share an unusually large number of such putative IBD blocks, we can hope to identify unknown relatives in genotyping datasets. In fact, companies like 23&me and Ancestry.com use signals of IBD to help identify family ties.

As another example, consider the case of third cousins. You share one of eight sets of great-great-grandparents with each of your (likely many) third cousins. On average, you and each of your third cousins each inherit one-sixteenth of your genome from each of those two great-great-grandparents. This turns out to imply that on average, a little less than one percent of your and your third cousin's genomes ($2 \times (\frac{1}{16})^2 =0.78\%$) will be identical by virtue of descent from those shared ancestors. A simulated example where third cousins share blocks of their genome (on chromosome 16 and 2) due to their great-great-grandmother is shown in {numref}`figure-2.15`.

:::{figure} ../figures/sharing_relatives/Third_cousin_overlap_1.png
:name: figure-2.15
:align: left

\- A simulation of sharing between third cousins, the details are the same as in {numref}`figure-2.14`.
:::

:::{margin}
Indeed the suspect in case of the Golden State Killer was identified through identifying third cousins that genetically matched a DNA sample from an old crime scene (see [here](https://gcbias.org/2018/05/07/how-lucky-was-the-genetic-investigation-in-the-golden-state-killer-case/) for more details).
:::

Note how if you compare {numref}`figure-2.15` and {numref}`figure-2.14`, individuals inherit less IBD from a shared great-great-grandmother than from a shared grandmother, as they inherit from more total ancestors further back. Also notice how the sharing occurs in shorter genomic blocks, as it has passed through more generations of recombination during meiosis. These blocks are still detectable, and so third cousins can be detected using high-density genotyping chips, allowing more distant relatives to be identified than single marker methods alone. More distant relations than third cousins, e.g. fourth cousins, start to have a significant probability of sharing none of their genome IBD. But you have many fourth cousins, so you will share some of your genome IBD with some of them; however, it gets increasingly hard to identify the degree of relatedness from genetic data the deeper in the family tree this sharing goes.

### Inbreeding

We can define an inbred individual as an individual whose parents are more closely related to each other than two random individuals drawn from some reference population.

:::{margin}
```{figure} ../figures/Child_first_cousins_Homozy_BD.png
:name: figure-2.16
:align: left

\- Alleles being transmitted through an inbred pedigree. The two sisters (mum and aunt) share two alleles identical by descent (IBD). The cousins share one pink  allele IBD. The offspring of first cousins is homozygous by descent at this locus for the pink allele. You can simulate more examples for yourself [here](https://cooplab.shinyapps.io/Shiny_pedigree_sims/).
```
:::

When two related individuals produce an offspring, that individual can receive two alleles that are identical by descent, i.e. they can be homozygous by descent (sometimes termed autozygous), due to the fact that they have two copies of an allele through different paths through the pedigree. This increased likelihood of being homozygous relative to an outbred individual is the most obvious effect of inbreeding. It is also the one that will be of most interest to us, as it underlies a lot of our ideas about inbreeding depression and population structure. For example, in {numref}`figure-2.16` our offspring of first cousins is homozygous by descent having received the same IBD allele via two different routes around an inbreeding loop.

As the offspring receives a random allele from each parent ($i$ and $j$), the probability that those two alleles are identical by descent is equal to the kinship coefficient $F_{ij}$ of the two parents (Equation {eq}`eq-2.5`). This follows from the fact that the genotype of the offspring is made by sampling an allele at random from each of our parents.

<center>
<div style="width: 50%;">

:::{table} - Generalized Hardy– Weinberg. Dropped the indices i and j for simplicity.
:name: table-2.5

| $f_{11}$        | $f_{12}$   | $f_{22}$        |
| :------:        | :------:   | :------:        |
| $(1-F)p^{2}+Fp$ | $(1-F)2pq$ | $(1-F)q^{2}+Fq$ |
:::

</div>
</center>

The only way the offspring can be heterozygous ($A_1 A_2$) is if their two alleles at a locus are not IBD (otherwise they would necessarily be homozygous). Therefore, the probability that they are heterozygous is

:::{math}
:label: eq-2.9

P(A_{1}A_{2}) = P(A_{1}A_{2} \text{ | I&J not IBD})P(\text{I&J not IBD}) = 2pq(1-F_{ij}),
:::

The offspring can be homozygous for the $A_1$ allele in two different ways. They can have two non-IBD alleles that are not IBD but happen to be of the allelic type $A_1$, or their two alleles can be IBD, such that they inherited allele $A_1$ by two different routes from the same ancestor. Thus, the probability that an offspring is homozygous for $A_1$ is

:::{math}
:label: eq-2.10

P(A_{1}A_{2}) = P(A_{1}A_{2} \text{ | I&J not IBD})P(\text{I&J not IBD})

+ P(A_{1}A_{2} \text{ | I&J IBD})P(\text{I&J IBD})

= p^{2}(1-F_{ij}+pF_{ij})
:::

using the Law of Total Probability (see Appendix Equation {eq}`eq-A.12`). Therefore, the frequencies of the three possible genotypes can be written as given in {numref}`table-2.5`, which provides a generalization of the Hardy--Weinberg proportions.

:::{admonition} Question 9
:name: question-2.9
The frequency of the $A_1$ allele is $p$ at a biallelic locus. Assume that our population is randomly mating and that the genotype frequencies in the population follow from HW. We select two individuals at random to mate from this population. We then mate the children from this cross. What is the probability that the child from this full sib-mating is homozygous?
:::

##### Multiple inbreeding loops in a pedigree.

Up to this point we have assumed that there is at most one inbreeding loop in the recent family history of our individuals, i.e. the parents of our inbred individual have at most one recent genealogical connection. However, an individual who has multiple inbreeding loops in their pedigree can be homozygous by descent thanks to receiving IBD alleles via multiple different different loops. To calculate inbreeding in pedigrees of arbitrary complexity, we can extend beyond our original relatedness coefficients $r_0$, $r_1$, and $r_2$ to account for higher order sharing of alleles IBD among relatives. For example, we can ask, what is the probability that *both* of the alleles in the first individual are shared IBD with one allele in the second individual? There are nine possible relatedness coefficients in total to completely describe kinship between two diploid individuals, and we won't go in to them here as it's a lot to keep track of. However, we will show how we can calculate the inbreeding coefficient of an individual with multiple inbreeding loops more directly.

Let's say the parents of our inbred individual (B and C) have $K$ shared ancestors, i.e. individuals who appear in both B and C's recent family trees. We denote these shared ancestors by $A_1, \dots,A_K$, and we denote by $n$ the total number of individuals in the chain from B to C via ancestor $A_i$, including B, C, and $A_i$. For example, if B is C's aunt, then B and C share two ancestors, which are B's parents and, equivalently, C's grandparents. In this case, there are n=4 individuals from B to C through each of these two shared ancestor. In the general case, the kinship coefficient of B and C, i.e. the inbreeding coefficient of their child, is

:::{math}
:label: eq-2.11
F = \sum_{i=1}^K \frac{1}{2^{n_i}}(1+f_{A_i})
:::

:::{margin}
For example, in the case of our aunt-nephew case, assuming that the aunt’s two parents are their only recent shared ancestors, then $F = \frac{1}{24} + \frac{1}{24} = \frac{1}{8}$, in agreement with the answer we would obtain from Equation {eq}`eq-2.5`.
:::

where $f_{A_i}$ is the inbreeding coefficient of the ancestor $A_i$. What's happening here is that we sum over all the mutually-exclusive paths in the pedigree through which B and C can share an allele IBD. With probability $\frac{1}{2^{n_i}}$, a pair of alleles picked at random from B and C is descended from the same ancestral allele in individual $A_i$, in which case the alleles are IBD. However, even if B inherits the maternal allele and C inherits the paternal allele of shared ancestor $A_i$, if $A_i$ was themselves inbred, with probability $f_{A_i}$ those two alleles are themselves IBD. Thus a shared *inbred* ancestor further increases the kinship of B and C.

:::{margin}
```{figure} ../illustration_images/alleles_genotypes/Carlos_second/378px-Juan_de_Miranda_Carreno_002.jpg
:name: figure-2.18
:align: left

Charles II of Spain (by Juan Carreño de Miranda, 1685). [Public Domain](https://it.wikipedia.org/wiki/Carlo_II_di_Spagna\#/media/File:Juan_de_Miranda_Carreno_002.jpg).
```
:::

:::{figure} ../Journal_figs/alleles_genotypes/Charles_second_pedigree/Carlos_second_pedigree_2_trimmed.pdf
:name: figure-2.17
:align: left

The pedigree of King Charles II of Spain. Pedigree from [wikimedia](https://commons.wikimedia.org/wiki/File:Carlos_segundo80.png) drawn by [Lec CRP1](https://en.wikipedia.org/wiki/User:Lec_CRP1), public domain.
:::

Multiple inbreeding loops increase the probability that a child is homozygous by descent at a locus, which can be calculated simply by plugging in $F$, the child's inbreeding coefficient, into our generalized HW equation.

As one extreme example of the impact of multiple inbreeding loops in an individual's pedigree, let's consider king Charles II of Spain, the last of the Spanish Habsburgs. Charles was the son of Philip IV of Spain and Mariana of Austria, who were uncle and niece. If this were the only inbreeding loop, then Charles would have had an inbreeding coefficient of $\frac{1}{8}$. Unfortunately for Charles, the Spanish Habsburgs had long kept wealth and power within their family by arranging marriages between close kin. The pedigree of Charles II is shown in {numref}`figure-2.17`, and multiple inbreeding loops are apparent. For example, Phillip III, Charles II's grandfather and great-grandfather, was himself a child of an uncle-niece marriage.

:::{margin}
Pedro Gargantilla, who performed Charles's autopsy, stated that his body "did not contain a single drop of blood; his heart was the size of a peppercorn; his lungs corroded; his intestines rotten and gangrenous; he had a single testicle, black as coal, and his head was full of water." While some of this description may refer to actual medical conditions, some of these details seem a little unlikely. See [here](https://www.thevintagenews.com/2017/03/23/when-charles-ii-of-spain-died-the-autopsy-stated-that-his-body-did-not-contain-a-single-drop-of-blood-and-his-head-was-full-of-water/).
:::

{cite:t}`alvarez:09` calculated that Charles II had an inbreeding coefficient of $0.254$, equivalent to a full-sib mating, thanks to all of the inbreeding loops in his pedigree. Therefore, he is expected to have been homozygous by descent for a full quarter of his genome. As we'll talk about later in these notes, this means that Charles may have been homozygous for a number of recessive disease alleles, and indeed he was a very sickly man who left no descendants due to his infertility. Thus plausibly the end of one of the great European dynasties came about through inbreeding.

### Calculating inbreeding coefficients from genetic data

If the observed heterozygosity in a population is $H_O$, and we assume that the generalized Hardy--Weinberg proportions hold, we can set $H_O$ equal to $f_{12}$, and solve Equation {eq}`eq-2.9` for $F$ to obtain an estimate of the inbreeding coefficient as

:::{margin}
The hat on the F here denotes that this is an estimate from a sample.
:::

:::{math}
:label: eq-2.12

\hat{F} = 1-\frac{f_{12}}{2pq} = \frac{2pq - f_{12}}{2pq}.
:::

As before, $p$ is the frequency of allele $A_{1}$ in the population. This can be rewritten in terms of the observed heterozygosity ($H_O$) and the heterozygosity expected in the absence of inbreeding, $H_E=2pq$, as

:::{math}
:label: eq-2.13

\hat{F} = \frac{H_E-H_O}{H_E} = 1 - \frac{H_O}{H_E}.
:::

Hence, $\hat{F}$ quantifies the deviation due to inbreeding of the observed heterozygosity from the one expected under random mating, relative to the latter.

:::{admonition} Question 10
:name: question-2.10
Suppose the following genotype frequencies were observed for an esterase
locus in a population of *Drosophila* (A denotes the "fast" allele and
B denotes the "slow" allele):

<center>
<div style="width: 50%;">

| $AA$ | $AB$ | $BB$ |
| :--: | :--: | :--: |
| 0.6  | 0.2  | 0.2  |

</div>
</center>

What is the estimate of the inbreeding coefficient at the esterase locus?
:::

:::{margin}
```{figure} ../illustration_images/alleles_genotypes/grey_wolf/5988399184_0c36a8e51c_z.jpg
:name: figure-2.19
:align: left

\- Grey wolf (*Canis lupus*). <span style="font-size: smaller;">Dogs, jackals, wolves, and foxes: a monograph of the Canidae. 1890. y J.G. Keulemans. Image from the [Biodiversity Heritage Library](https://www.biodiversitylibrary.org/page/17002968\#page/58/mode/1up). Contributed by University of Toronto - Gerstein Science Information Centre. Not in copyright.</span>
```
:::

If we have multiple loci, we can replace $H_O$ and $H_E$ by their means over loci, $\bar{H}_O$ and $\bar{H}_E$, respectively. Note that, in principle, we could also calculate $F$ for each individual locus first, and then take the average across loci. However, this procedure is more prone to introducing a bias if sample sizes vary across loci, which is not unlikely when we are dealing with real data.

Genetic markers are commonly used to estimate inbreeding for wild and/or captive populations of conservation concern. As an example of this, consider the case of the Mexican wolf (*Canis lupus baileyi*), a sub-species of gray wolf.

They were extirpated in the wild during the mid-1900s due to hunting, and the remaining five Mexican wolves in the wild were captured to start a breeding program. @vonHoldt:11 estimated the current-day, average expected heterozygosity to be $0.18$, based on allele frequencies at over forty thousand SNPs. However, the average Mexican wolf individual was only observed to be heterozygous at $12\%$ of these SNPs. Therefore, the average inbreeding coefficient for the Mexican wolf is $\hat{F} = 1 -\frac{0.12}{0.18}$, i.e. $\sim 33 \%$ of a lobo's genome is homozygous due to recent inbreeding in their pedigree.

##### Genomic blocks of homozygosity due to inbreeding.

As we saw above, close relatives are expected to share alleles IBD in large genomic blocks. Thus, when related individuals mate and transmit alleles to an inbred offspring, they transmit these alleles in big blocks through meiosis. As an example, let's return to the case of our hypothetical first cousins from {numref}`figure-2.10`. If this pair of individuals had a child, one possible pattern of genetic transmission is shown in {numref}`figure-2.20`. The child has inherited the red stretch of chromosome via two different routes through their predigree from the grandparents. This is an example of an autozygous segment, where the child is homozygous by descent at all of the loci in this red region. The inbreeding coefficient of the child sets the proportion of their genome that will be in these autozygous segments. For example, a child of first full cousins is expected to have $1/16$ of their genome in these segments. The more distant the loop in the pedigree, the more meioses that chromosomes have been through and the shorter individual blocks will be. A child of first cousins will have longer blocks than a child of second cousins, for example.

:::{figure} ../figures/sharing_relatives/first_cousin_offspring.pdf
:name: figure-2.20
:align: left

\- A pedigree showing the offspring of first cousins. The chromosomes of their great-grandparents are coloured different colours so their transmission can be tracked. The child is homozygous by descent (HBD) for a section of the red chromosome.
:::

Individuals with multiple inbreeding loops in their family tree can have a high inbreeding coefficient due to the combined effect of many small blocks of autozygosity. For example, Charles II had an inbreeding coefficient that is equivalent to that of the child of full-sibs, with a quarter of his genome expected to homozygous by descent, but this would be made up of many shorter blocks.

We can hope to detect these blocks by looking for unusually long genomic runs of homozygosity (ROH) sites in an individual's genome. One way to estimate an individual's inbreeding coefficient is then to total up the proportion of an individual's genome that falls in such ROH regions. This estimate is called $F_{ROH}$.

An example of using $F_{ROH}$ to study inbreeding comes from the work of {cite:t}`sams2018fine`, who identified runs of homozygosity in 2,500 dogs, ranging from 500kb up to many megabases. {numref}`figure-2.21` shows the distribution of $F_{ROH}$ of individuals in each dog breed for the X and autosome. In {numref}`figure-2.23` this is broken down by the length of ROH segments.

:::{figure} ../figures/sharing_relatives/dogs_FROH.pdf
:name: figure-2.21
:align: left

\- The distribution of FROH of individuals from various dog breeds from {cite:t}`Sams:18`, licensed under CC BY 4.0.
:::

:::{margin}
```{figure} ../illustration_images/alleles_genotypes/english_bulldog/14752595581_4330377c97_z.jpg
:name: figure-2.22
:align: left

\- English bulldog. The dogs of Boytown. 1918. Dyer, W. A.
```
:::

:::{figure} ../figures/sharing_relatives/dog_FROH_dist.pdf
:name: figure-2.23
:align: left

\- Cumulative density of ROH length, measured in megabases (Mb) from {cite:t}`Sams:18` for various dog breeds (li- censed under CC BY 4.0). Note that longer lengths of ROH are on the left of the plot.
:::

Dog breeds have been subject to intense breeding that has resulted in high levels of inbreeding. Of the population samples examined, Doberman Pinschers have the highest levels of their genome in runs of homozygosity ($F_{ROH}$), somewhat higher than English bulldogs. In {numref}`figure-2.23`, we can see that English bulldogs have more short ROH than Doberman Pinschers, but that Doberman Pinschers have more of their genome in very large ROH ($>16 Mb$). This suggests that English bulldogs have had long history of inbreeding as they have many small blocks, but that Doberman Pinschers have a lot of recent inbreeding as their autozygosity is contained in long blocks relatively unbroken by recombination.

:::{admonition} Chapter Summary
:name: chapter-summary-2
- This chapter developed the relationship between allele frequencies and genotype frequencies within a generation and among relatives.
- Under random mating, we derived expectations of the genotype frequencies (Hardy-Weinberg), and we can identify deviations away from these expectations.
- Identity by descent (IBD) refers to the sharing of alleles due to a recent shared biological relationship.
- We can predict the probability and expected level of sharing of alleles IBD among pairs of relatives using mendelian transmission probabilities (as contained in coefficients $r_0$, $r_1$, and $r_2$). One useful summary of relatedness for a pair of individuals is the kinship coefficient $F_{i,j}$.
- We can also learn about genetic relationships from the sharing of genomic segments among relatives, with many long shared segments revealing a closer relationship.
- An inbred individual has parents who are more closely related than random draws from some reference population.
- Inbreeding results in decreased heterozygosity and a complementary increase in homozygosity. We can use the kinship coefficient of the parents to estimate the distortion away from Hardy-Weinberg and the expected level of heterozygosity.
- Inbreeding coefficients can be calculated from genetic data, either for multiple individuals at a single locus or for multiple loci for a single individual.
:::

:::{admonition} Question 11
:name: question-2.11

Calculate $r_{0}$, $r_{1}$, $r_{2}$ and the coefficient of kinship $F$
between:

- **A)** A grandparent and their grandchild
- **B)** A great grandparent and their great grandchild
- **C)** Full siblings
- **D)** A great aunt and her grand nephew (your great aunt = your parent's aunt)
:::

:::{admonition} Question 12
:name: question-2.12

You are studying a codominant flower colour polymorphism. Skipping
through a meadow of flowers you and compile the following data:

| red | pink | white |
| :-: | :--: | :---: |
| 200 | 100  | 200   |

- **A)** What frequencies would you expect at this locus under Hardy-Weinberg equilibrium?
- **B)** Calculate the inbreeding coefficient at this locus.
- **C)** Name two distinct processes that could lead to the deviation you see, and describe how they would result in a deficit of heterozygotes.
:::

:::{admonition} Question 13
:name: question-2.13

What are the relatedness coefficients of the X chromosome between:

- **A)** Two male full siblings?
- **B)** Two female full siblings?
- **C)** What is the probability that a female offspring of a full sib mating is homozygous by descent at a locus on her X chromosome?
:::

:::{admonition} Question 14
:name: question-2.14

You are studying the wing spot polymorphism in a butterfly species. From crosses in the lab you find that the presence of wing spots is determined by a dominant allele.

You collect 100 butterflies, 84 of them have the wing spots. What is the frequency of the wing-spot allele? What assumption did you have to make to come to your answer?
:::

:::{admonition} Question 15
:name: question-2.15

An allele has frequency of $0.001$ in the population. What is the probability that both you and your first (full) cousin are heterozygote for the allele?
:::

:::{admonition} Question 16
:name: question-2.16

The kinship coefficient of the parents is the inbreeding coefficient of the offspring. Explain, with reference to the weighting of relatedness coefficients in the inbreeding coefficient, why the inbreeding coefficient is the probability that a locus is homozygous by descent.
:::

:::{admonition} Question 17
:name: question-2.17

In terms of identity by descent, explain why multiple inbreeding loops in an individual's pedigree lead to higher levels of inbreeding.
:::