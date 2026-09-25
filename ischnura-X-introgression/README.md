# Restricted X chromosome introgression in hybridizing damselflies

Analysis code for:

> Swaegers J, Sánchez-Guillén RA, Chauhan P, Wellenreuther M, Hansson B (2022).
> Restricted X chromosome introgression and support for Haldane's rule in hybridizing damselflies.
> *Proceedings of the Royal Society B* 289: 20220968. https://doi.org/10.1098/rspb.2022.0968

## About the study

*Ischnura elegans* has been expanding its range into northern and western Spain, where it
hybridizes with its sister species *I. graellsii*. Using genome-wide SNPs from RAD sequencing
of 253 damselflies, we tested whether the X chromosome introgresses less than the autosomes
between the two species, and whether hybrid males (X0) are underrepresented, as predicted
by Haldane's rule. Both predictions were supported: X-linked SNPs showed steeper genomic
clines and less introgression in ABBA-BABA analyses than autosomal SNPs, and males were
underrepresented among admixed individuals.

## Workflow

The pipeline runs from demultiplexed reads to the introgression tests. Steps 1–8 are Bash
scripts for a SLURM cluster (the analyses were run on UPPMAX); the R scripts handle SNP
filtering, sexing and the statistical tests.

| Step | Script | What it does |
|---|---|---|
| 1 | `scripts/01_clone_filter.slurm` | Removes PCR duplicates (Stacks `clone_filter`) |
| 2 | `scripts/02_align_bowtie2.slurm` | Aligns reads to the *I. elegans* genome (Bowtie2) and sorts BAMs |
| 3 | `scripts/03_ref_map.slurm` | Builds loci and calls SNPs (Stacks `ref_map`), with and without outgroup species |
| 4 | `scripts/04_populations.slurm` | Filters SNPs into the full and outgroup SNP sets (Stacks `populations`) |
| 5 | `scripts/05_diagnostic_snps.slurm` | Finds fixed differences between allopatric populations and genotypes them in all samples |
| 6 | `R/01_full_snp_set.R` | Genotypic sexing from X homozygosity; removes X-linked SNPs heterozygous in males; one random SNP per RAD-tag |
| 7 | `R/02_diagnostic_snp_set.R` | Same filtering for the diagnostic SNP set |
| 8 | `scripts/06_split_X_autosomes.sh` | Splits VCFs into X-linked and autosomal SNPs |
| 9 | `scripts/07_pca_admixture.sh` | PCA (PLINK) and supervised ADMIXTURE, with males treated as haploid on the X |
| 10 | `scripts/08_dsuite.sh` | ABBA-BABA statistics in sliding windows (Dsuite `Dinvestigate`) |
| 11 | `R/03_bgc.R` | Bayesian genomic clines: X-linked vs autosomal alpha and beta (figure 2) |
| 12 | `R/04_abba_baba.R` | X-linked vs autosomal f_dM, D and f_d; permutation and Wilcoxon tests (figure 3) |

Shared R functions (reading VCFs, the male-homozygosity filter, SNP thinning and the
permutation test) are in `R/functions.R`.

## Key methods

- **X-linked SNPs** were identified using the X-linked scaffolds of the *I. elegans* genome
  assembly (Chauhan et al. 2021, *Genomics*).
- **Sexing.** Ischnura males carry a single X chromosome (X0), so they are homozygous at all
  X-linked SNPs. Individuals were classified as male when homozygosity at X-linked SNPs was
  above 0.96, the valley of the bimodal distribution.
- **Filtering X-linked SNPs.** X-linked SNPs where any male was heterozygous were removed,
  as these calls cannot be correct in a hemizygous sex.
- **Permutation tests.** To compare X-linked and autosomal SNPs with equal sample sizes,
  10,000 random autosomal subsets the size of the X-linked set were drawn and their medians
  compared with the observed X-linked median.

## Data

Raw reads are available from the NCBI Sequence Read Archive (PRJNA850104). VCF files and
sample metadata are on Dryad: https://doi.org/10.5061/dryad.gqnk98sp8. The genome assembly
is described in Chauhan et al. (2021), https://doi.org/10.1016/j.ygeno.2021.04.003.

Data files are not included in this repository. The scripts expect:

- `config.sh`: paths for the cluster scripts (edit before running)
- `data/` for the R scripts: the Stacks VCFs, `X_chrom_ids` and `Autosome_ids`
  (one scaffold name per line), a sample information table (`ID`, `PopID`, sex), the bgc
  output and the Dsuite output folders

## Software

Stacks 2.2, Bowtie2 2.3, SAMtools, VCFtools, PLINK 1.9, ADMIXTURE 1.3.0, Dsuite, bgc,
and R with adegenet, vcfR, dplyr, ggplot2 and ClineHelpR.

## Running

Cluster steps are submitted from the repository root, for example:

```bash
sbatch -A <account> scripts/01_clone_filter.slurm
```

R scripts are run from the repository root as well:

```bash
Rscript R/01_full_snp_set.R
```

## Notes

The scripts were reorganised from my working notes and R notebooks for this repository:
hard-coded paths were replaced by settings, repeated code was turned into functions, and
comments were added. Parameters and filtering steps are those used for the paper.
The one-SNP-per-RAD-tag subset is drawn at random; the published analysis did not fix a
random seed, so a rerun gives a slightly different subset.

Not included here: demultiplexing (`process_radtags`), preparation and running of bgc,
and the plots and tests for figures 1 and 4 (admixture proportions and Haldane's rule).

## Author

Janne Swaegers (ORCID [0000-0003-1952-3170](https://orcid.org/0000-0003-1952-3170))

## License

MIT; see `LICENSE`.
