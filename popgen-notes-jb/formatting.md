# Formatting

This file contains the formatting to be used within the Jupyter Notebook.

## Citations

- Within text citation, year in parentheses - {cite:t}`janicke:16`
- Parentheses citation - {cite:p}`hardy1908mendelian`
    - If multiple citations, separate with a comma - {cite:p}`hardy1908mendelian,weinberg1908ber`

## References

- Equations - equation {eq}`eq-4.13`
- Figures - {numref}`figure-4.9`

Note that the figures reference automatically adds "Fig." to the reference.

## Math

:::{math}
:label: eq-4.14
    \theta = 4N\mu
:::

## Figures

```{figure} ../../figures/Cousins_IBD_chromo_cartoon.png
:name: figure-2.10
:align: left

\- First cousins sharing a stretch of chromosome identical by descent. The different grandparental diploid chromosomes are coloured so we can track them and recombinations between them across the generations. Notice that the identity by descent between the cousins persists for a long stretch of chromosome due to the limited number of generations for recombination. The squares represent males and circles females.
```

## Margins
:::{margin}
{cite:p}`cotterman:40,malecot:48`
:::

:::{margin}
```{figure} ../../figures/sharing_relatives/IBD_0_1_2.pdf
:name: figure-2.11
:align: left

\- A pair of diploid individ- uals (i and j) sharing 0, 1, or 2 alleles IBD where lines show the sharing of alleles by descent (e.g. from a shared ancestor). Here we’ll focus on IBD of outbred individuals. Dealing with sharing between inbred individuals requires 6 more identity-by-descent r coefficients, which honestly makes my head spin.
```
:::