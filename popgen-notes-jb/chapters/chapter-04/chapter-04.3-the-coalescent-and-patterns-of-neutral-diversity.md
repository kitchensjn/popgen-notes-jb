# The Coalescent and patterns of neutral diversity

*"Life can only be understood backwards; but it must be lived forwards"* - Kierkegaard

### Pairwise Coalescent time distribution and the number of pairwise differences

Thinking back to our calculations we made about the loss of neutral
heterozygosity and equilibrium levels of diversity (in Sections {ref}`sec-4.1` and {ref}`sec-4.1.1`), you'll note that we could first
specify which generation a pair of sequences coalesce in, and then
calculate some properties of heterozygosity based on that. That's
because neutral mutations do not affect the probability that an
individual transmits an allele, and so don't affect the way in which we
can trace ancestral lineages back through the generations.

As such, it will often be helpful to consider the time to the common
ancestor of a pair of sequences ($T_2$), and then think of the impact of
that time to coalescence on patterns of diversity. See {numref}`figure-4.15` for an example of this.

```{figure} ../../figures/Coalescent.png
:name: figure-4.15
:align: left

\- A simple demonstration of the coalescent process. The simulation
consists of a diploid population of 10 individuals (20 alleles). In each
generation, each individual is equally likely to be the parent of an
offspring (and the allele transmitted is indicated by a light grey
line). We track a pair of alleles, chosen in the present day, back 14
generations until they find a common ancestor. Deeper in time than 14
generations those two alleles have the same ancestral lineage and
completely share their history, e.g. the mutations that occur on that
lineage.
```

The probability that a pair of alleles have failed to coalesce in $t$
generations and then coalesce in the $t+1$ generation back is

:::{math}
:label: eq-4.21
    P(T_2=t+1) = \frac{1}{2N} \left(1- \frac{1}{2N} \right)^{t}
:::

For example, the probability that a pair of sequences coalesce three
generations back is the probability that they fail to coalesce in
generation 1 and 2, which is $\left(1- \frac{1}{2N} \right) \times \left(1- \frac{1}{2N} \right)$,
multipled by the probability that they find a common ancestor, i.e.
coalesce, in the third generation, which happens with probability
$\frac{1}{2N}$.

From the form of equation {eq}`eq-4.21`, we can see that the coalescent time of
our pair of alleles is a Geometrically distributed random variable,
where the probability of success is $p=\frac{1}{2N}$. The waiting
time for a pair of lineages to coalesce is like the number of tails
thrown while waiting for a head on a coin with the probability of a head
is $\frac{1}{2N}$, i.e. if the population is large we might be
waiting for a long time for our pair to coalesce. We'll denote this
geometric distribution by $T_2 \sim  \text{Geo}(1/(2N))$. The expected
(i.e. the mean over many replicates) coalescent time of a pair of
alleles is then


:::{math}
:label: eq-4.22
    E(T_2) = 2N
:::

generations. This form to the
expectation follows from the fact that the mean of an geometric random
variable is $\frac{1}{p}$.\
Conditional on a pair of alleles coalescing $t$ generations ago, there
are $2t$ generations in which a mutation could occur. See {numref}`figure-4.16` for an example. If the per generation
mutation rate is $\mu$, then the expected number of mutations between a
pair of alleles coalescing $t$ generations ago is $2 t\mu$ (the alleles
have gone through a total of $2t$ meioses since they last shared a
common ancestor).

So we can write the expected number of mutations ($S_2$) separating two
alleles drawn at random from the population as

:::{math}
:label: eq-4.23
    \begin{aligned}
E(S_2) &= \sum_{t=0}^{\infty} E(S_2 | T_2=t) P(T_2=t) \nonumber\\
& =\sum_{t=0}^{\infty} 2 \mu t P(T_2=t) \nonumber\\
& =2\mu E(T_2)  \nonumber\\
& = 4 \mu N \end{aligned}
:::

this makes use of the law of total
expectation (see Appendix equation {eq}`eq-A.27`) to average which generation our pair
of sequences coalesce in. We'll assume that mutation is rare enough that
it never happens at the same basepair twice, i.e. no multiple hits, such
that we get to see all of the mutation events that separate our pair of
sequences. This is assumption that repeat mutation is vanishingly rare
at a basepair is called the *i*nfinitely-many-sites assumption, which
should hold if $N\mu_{BP} \ll 1$, where $\mu_{BP}$ is the mutation rate
per basepair. Thus the number of mutations between a pair of sites is
the observed number of differences between a pair of sequences. In the
previous chapter we denote the observed number of pairwise differences
at putatively neutral sites separating a pair of sequences as $\pi$ (we
usually average this over a number of pairs of sequences for a region).
Therefore, under our simple, neutral, constant population-size model we
expect

:::{math}
:label: eq-4.24
    E(\pi) = 4 N \mu = \theta
:::

So we
can get an empirical estimate of $\theta$ from $\pi$, let's call this
$\widehat{\theta}_{\pi}$, by setting $\widehat{\theta}_{\pi}=\pi$, i.e.
our observed level of pairwise genetic diversity. If we have an
independent estimate of $\mu$, then from setting
$\pi =\widehat{\theta}_{\pi} = 4N\mu$ we can furthermore obtain an
estimate of the population size $N$ that is consistent with our levels
of neutral polymorphism. If we estimate the population size this way, we
should call it the effective coalescent population size ($N_e$). It's
best to think about $N_{e}$ estimated from neutral diversity as a
long-term effective population size for the species, but there are many
caveats that come along with that assumption. For example, past
bottlenecks and population expansions are all subsumed into a single
number and so this estimated $N_{e}$ may not be very representative of
the population size at any time. That said, it's not a bad place to
start when thinking about the rate of genetic drift for neutral
diversity in our population over long time-periods.

Let's take a moment to distinguish our expected heterozygosity (equation {eq}`eq-4.13`) from our expected number of pairwise
differences ($\pi$). Our expected heterozygosity is the probability that
two alleles at a locus, sampled from a population at random, are
different from each other. If one or more mutations have occurred since
a pair of alleles last shared a common ancestor, then our sequences will
be different from each other. On the other hand, our $\pi$ measure keeps
track of the average total number of differences between our loci. As
such, $\pi$ is often a more useful measure, as it records the number of
differences between the sequences, not just whether they are different
from each other (however, for certain types of loci, e.g.
microsatellites, heterozygosity is often used as we cannot usually count
up the minimum number of mutations in a sensible way). In the case where
our locus is a single basepair, the two measures will usually be close
to one another, as $H \approx \theta$ for small values of $\theta$. For
example, comparing two sequences at random in humans,
$\pi \approx 1/1000$ per basepair, and the probability that a specific
base pair differs between two sequences is $\approx 1/1000$. However,
these two quantities start to differ from each other when we consider
regions with higher mutation rates. For example, if we consider a 10kb
region, our mutation rate will 10,000 times larger than a single base
pair. For this length of sequence the probability that two randomly
chosen haplotypes differ is quite different from the number of
mutational differences between them. (Try a mutation rate of $10^{-8}$
per base and a population size of $10,000$ in our calculations of
$E[\pi]$ and H to see this.)

:::{admonition} Question 7
:name: question-4.7

{cite:t}`robinson:16` found that the endangered Californian Channel Island fox on
San Nicolas had very low levels of diversity
($\pi =0.000014 \text{bp}^{-1}$) compared to its close relative the
California mainland gray fox ($0.0012\text{bp}^{-1}$).\
**A)** Assuming a mutation rate of $2\times 10^{-8}$ per bp, what
effective population sizes do you estimate for these two populations?\
**B)** Why is the effective population size of the Channel Island fox so
low? \[Hint: quickly google Channel island foxes to read up on their
history, also to see how ridiculously cute they are.\]
:::

:::{admonition} Question 8
:name: question-4.8

In your own words describe why the coalescent time of a pair of lineages
scales linearly with the (effective) population size.
:::

### More details on the pairwise coalescent and the randomness of mutation

We found that our pairwise coalescent times followed a Geometric
distribution, equation {eq}`eq-4.21`. However, that assumes discrete
generations, and we'll often was to think about populations that lack
discrete generations (i.e. individuals reproducing at random times with
some mean generation time). Using our exponential approximation, we can
see that is 

:::{math}
:label: eq-4.25
    \approx \frac{1}{2N} e^{-t/(2N)}
:::

and so think of a
continuous random variable, i.e. we could say that the coalescent time
of a pair of sequences ($T_2$) is approximately exponentially
distributed with a rate $1/(2N)$, i.e.
$T_2 \sim \text{Exp}\left( 1/(2N) \right)$. Formally we can do this by
taking the limit of the discrete process more carefully. See Appendix
equation {eq}`eq-A.34` for more on exponential random variables.

We've derived the expected number of differences between a pair of
sequences and talked about the variability of the coalescent time for a
pair of sequences. The mutation process is also very variable; even if
two sequences coalesce in the very distant past by chance, they may
still be identical in the present if there was no mutation during that
time.

Conditional on the coalescent time $t$, the probability that our pair of
alleles are separated by $S_2$ mutations since they last shared a common
ancestor is binomially distributed


:::{math}
:label: eq-4.25
    P(S_2 | T_2 = t ) = {2t \choose j} \mu^{j} (1-\mu)^{2t-j}
:::

i.e.
mutations happen in $j$ generations and do not happen in $2t-j$
generations (with ${2t \choose j}$ ways this combination of events can
possibly happen). See Appendix equation {eq}`eq-A.28` for discussion of the binomial
distribution. Assuming that $\mu \ll 1$ and that $2t-j \approx 2t$, then
we can approximate the probability that we have $S_2$ mutations as a
Poisson distribution:

:::{math}
:label: eq-4.25
    P(S_2 | T_2 = t ) = \frac{ (2 \mu t )^{j} e^{-2\mu t}}{j!}
:::

i.e. a
Poisson with mean $2\mu t$. This is an example of taking the binomial
distribution to its Poisson distribution limit, see Appendix equation {eq}`eq-A.32` for more details. We'll not make much
use of this result, but it is very useful in thinking about how to
simulate the process of mutation.