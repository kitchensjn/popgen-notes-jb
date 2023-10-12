# Allele Frequencies

Allele frequencies are a central unit of population genetics analysis, but from diploid individuals we only get to observe genotype counts. Our first task then is to calculate allele frequencies from genotype counts. Consider a diploid autosomal locus segregating for two alleles ($A_1$ and $A_2$). We'll use these arbitrary labels for our alleles, merely to keep this general. Let $N_{11}$ and $N_{12}$ be the number of $A_1A_1$ homozygotes and $A_1A_2$ heterozygotes, respectively. Moreover, let $N$ be the total number of diploid individuals in the population. We can then define the relative frequencies of $A_1A_1$ and $A_1A_2$ genotypes as $f_{11} = \frac{N_{11}}{N}$ and $f_{12} = \frac{N_{12}}{N}$, respectively. The frequency of allele $A_1$ in the population is then given by

:::{math}
:label: eq-2.1

p = \frac{2 N_{11} + N_{12}}{2N} = f_{11} + \frac{1}{2} f_{12}.
:::

Note that this follows directly from how we count alleles given individuals' genotypes, and holds independently of Hardy--Weinberg proportions and equilibrium (discussed below). The frequency of the alternate allele ($A_2$) is then just $q=1-p$.

## Measures of genetic variability

#### Nucleotide diversity ($\pi$)

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

#### Number of segregating sites

Another measure of genetic variability is the total number of sites that are polymorphic (segregating) in our sample. One issue is that the number of segregating sites will grow as we sequence more individuals (unlike $\pi$). Later in the course, we'll talk about how to standardize the
number of segregating sites for the number of individuals sequenced (see **\eqn \eqref{watterson_theta}**).

#### The frequency spectrum

We also often want to compile information about the frequency of alleles across sites.  We call alleles that are found once in a sample *singletons*, alleles that are found twice in a sample *doubletons*, and so on. We count up the number of loci where an allele is found $i$ times out of $n$, e.g. how many singletons are there in the sample, and this is called the *frequency spectrum*. We'll want to do this in some consistent manner, such as calculating the frequency spectrum of the minor allele or the derived allele.

:::{admonition} Question 2
:name: question-2.2
How many minor-allele singletons are there in *D. simulans* in the ADH region? [Defining minor allele just within *D. simulans*.]
:::

#### Levels of genetic variability across species

Two observations have puzzled population geneticists since the inception of molecular population genetics. The first is the relatively high level of genetic variation observed in most obligately sexual species. This first observation, in part, drove the development of the Neutral theory of molecular evolution, the idea that much of this molecular polymorphism may simply reflect a balance between genetic drift and mutation.

The second observation is the relatively narrow range of polymorphism across species with vastly different census sizes. This observation represented a puzzle as the Neutral theory predicts that levels of genetic diversity should scale with population size. Much effort in theoretical and empirical population genetics has been devoted to trying to reconcile models with these various observations. We'll return to discuss these ideas throughout our course.

:::{margin}
```{figure} ../../illustration_images/alleles_genotypes/Ciona_intestinalis/21016139168_2a8a57ded3_z.jpg
---
name: figure-2.2
align: left
---
\- Sea Squirt (*Ciona intestinalis*). <span style="font-size: smaller;">Einleitung in die vergleichende gehirnphysi- ologie und Vergleichende psychologie. Loeb, J. 1899. Image from the Biodiversity Heritage Library. Contributed by MBLWHOI Library. No known copyright restrictions.</span>
```
:::

The first observations of molecular genetic diversity within natural populations were made from surveys of allozyme data, but we can revisit these general patterns with modern data. For example, {cite:p}`leffler:12` compiled data on levels of within-population, autosomal nucleotide diversity ($\pi$) for 167 species across 14 phyla from non-coding and synonymous sites ({numref}`figure-2.3`). The species with the lowest levels of $\pi$ in their survey was Lynx, with $\pi = 0.01\%$, i.e. only $\frac{1}{10000}$ bases differed between two sequences. In contrast, some of the highest levels of diversity were found in *Ciona savignyi*, Sea Squirts, where a remarkable $\frac{1}{12}$ bases differ between pairs of sequences. This 800-fold range of diversity seems impressive, but census population sizes have a much larger range.

```{figure} ../../Journal_figs/alleles_genotypes/Leffer_riddle/Leffer_riddle_diversity.png
---
name: figure-2.3
align: left
---
\- Levels of autosomal nucleotide diversity for 167 species across 14 phyla. Figure 1 from {cite:p}`leffler:12`, licensed under CC BY 4.0. Points are ranked by their $\pi$, and coloured by their phylum. Note the log-scale.
```

:::{margin}
```{figure} ../../illustration_images/alleles_genotypes/Lynx/20731949565_8a065700af_z.jpg
---
name: figure-2.4
align: left
---
\- Eurasian Lynx (*Lynx lynx*). <span style="font-size: smaller;">An introduction to the study of mammals living and extinct. Flower, W.H. and Lydekker, R. 1891. Image from the [Biodiversity Heritage Library](https://www.flickr.com/photos/internetarchivebookimages/20731949565/in/photolist-x5Jzv2-x6QVyp-xir9rH-wYHrQD-wPn1sP-w9PsqY-xDcqri-sMcQoB-trrkVd-x6Nx1H-wPea7N-sM28N9-tJ3zsp-xneVdx-wGJRtQ-xnfHZ8-wPfga7-xCUPrN-x7kXDV-xmAb9E-xm3x4k-xBoSKb-wGTgyB-xBoSbf-wGGvzA-xmzYTJ-oeKJcH-xA1Ffr-xA1Eji-xqWTQZ-xF4Lru-oxJfrH-x7ojSn-xra8zP-wGGibY-xgb21y-xY1jH9-xY1iyf-wGHxTS-wGQEoR-xmtPQh-x8uFKK-xdTkoU-wPggQf-wPfvHN-wPfc27-w9YnGF-wPeauS-wPiVxK-w6aiSi). Contributed by Cornell University Library. No known copyright restrictions.</span>
```
:::

## Hardy--Weinberg proportions

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
```{figure} ../../illustration_images/alleles_genotypes/Kermode_bear/EIiK2dsWoAAmUOf.jpeg
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

```{figure} ../../figures/CEU_YRI_separately_HWE.png
:name: figure-2.6
:align: left

\- Demonstrating Hardy--Weinberg proportions using 10,000 SNPs from the HapMap European (CEU)  and African (YRI) populations. Within each of these populations the allele frequency against the frequency of the 3 genotypes; each SNP is represented by 3 different coloured points. The solid lines show the mean genotype frequency. The dashed lines show the predicted genotype frequency from Hardy--Weinberg equilibrium. Code [here](https://github.com/cooplab/popgen-notes/blob/master/Rcode/HWE_exercise/HWE_HAPMAP.R). Blog post on figure [here](http://gcbias.org/2011/10/13/population-genetics-course-resources-Hardy--Weinberg-eq/). 
```

## Assortative mating

One major violation of the assumptions of Hardy Weinberg is non-random mating with respect to the genotype at a locus. One way that individuals can mate non-randomly is if individuals choose to mate based on a phenotype determined by (in part) the genotype at a locus. This non-random mating can be between: 1) individuals with similar phenotype, so called positive assortative mating or 2) individuals with dissimilar phenotypes, negative assortative mating or disassortative mating. Here we'll briefly discuss a couple of real examples of assortative mating to make sure we're all on the same page. We'll encounter other forms of non-random mating, due to inbreeding and population structure, in the next few chapters.

:::{figure} ../../Journal_figs/alleles_genotypes/Heliconius_Merrill_assort_mating/trimmed_Heliconius_Merrill_assort_mating.png
:name: figure-2.7
:align: left

\- Wing pattern phenotypes of top, *H. cydno chioneus* (left), *H. melpomene rosina* (right), their nonmimetic first-generation hybrid (center); and bottom, their sympatric comimics *H. sapho sapho* (left) and *H. erato demophoon* (right). Figure and caption modified from {cite:p}`merrill2019genetic`, licensed under CC BY 4.0.
:::

Positive assortative mating on the basis of a phenotype can create an excess of homozygotes. Heliconius butterflies are famous for their mimicry, where poisonous pairs of distantly related species mimic each others' bright colour patterns and so share the benefits of being avoided by visual predators (Müllerian mimics). *H. melpomene rosina* and *H. cydno chioneus* are closely related species that co-occur in central Panama, but mimic different other co-occuring species ({numref}`figure-2.7`). These differences in colouration pattern are due to a few loci with large phenotypic effects. The two species can hybridize and produce viable F1 hybrids. These F1 hybrids are heterozygotes at the colour loci, and their intermediate appearance means that they're poor mimics and so are quickly eaten by predators. However, these heterzygote (F1) hybrids are very rare in nature $<\frac{1}{1000}$, as the parental species show strong positive assortatively mating based on colour pattern, based on genetic differences in mate preference {cite:p}`merrill2019genetic`.

:::{margin}
```{figure} ../../illustration_images/alleles_genotypes/White_throated_sparrow/white_throated_sparrow.jpg
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
```{figure} ../../illustration_images/alleles_genotypes/Sparassis_crispa/Sparassis_crispa_cauliflower_mushroom.pdf
:name: figure-2.9
:align: left

\- Cauliflower mushrooms (*Sparassis crispa*) parasitize tree roots and form these amazing, edible fruiting bodies, which can weigh in at up to 30lb and apparently taste like noodles. In a collection of 18 fruiting bodies from a *Sparassis* population, all individuals were heterozygotes for mating type and 17 different mating types were genetically identified {cite:p}`james2015mushrooms,martin1978decay`. <span style="font-size: smaller;">Atlas champignons comestibles et vénéneux (1891 ). Dufour, L Image from the Biodiversity Heritage Library. Contributed by New York Botanical Garden. Not in copyright.</span>
```
:::

There are almost no 2m homozygotes (so called Super white individuals) despite the 2m allele being common in the population (data from {cite:p}`tuttle2016divergence`, table S1).

Another important example of disassortative mating are mating type systems, which are present in many fungi, algae, and protozoa. Gametes of the same species can only fuse to form a zygote if they differ in mating type. The mating type of gametes is genetically controled by a mating type locus, and so individuals are nearly always heterozygous at this locus. In some groups of organisms, there are just two different alleles, in other clades these loci have tens or hundreds of alleles.