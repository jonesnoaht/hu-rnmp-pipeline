# Preferential rNMP incorporation and failed excision by gene

**v1 scope.** The pipeline v1 is wild-type versus RNase H2 (`rnh201`) rNMP mapping.
This note is an rnh201 HydEn-seq gene table, not that comparison. Hydroxyurea is
Future work / v2 (Meroni et al. 2019, PMID 30847483). v1 does not assign
replicative polymerases.

SGD GFF gene spans intersected with Zhou 2019 HydEn-seq on sacCer3.
Libraries: SRR8503261 (rnh201 NT, no RNase HII pre-cleavage) and
SRR8503270 (rnh201 RNase HII). Same rnh201 background. Koh ribose-seq
is not in this table: rnh201 vs rnh1+rnh201 is underpowered and is a
different assay.

## What the ratio is

RNase HII / NT is the within-assay contrast for rNMPs that remained
embedded in rnh201 cells. The HII library records ends after in vitro
RNase HII cleavage. The NT library records alkali / hydrolytic ends
without that cleavage. These are separate libraries, not a before/after
subtraction of one DNA sample. Pipeline coordinates already apply the
1-nt Alk vs RHII offset (Ribose-Map coordinate.sh); this intersection
does not shift them again.

The ranked ratio is unique HII positions / unique NT positions on the
SGD gene span. Read-end counts are a depth column. Genome-wide, HII
has fewer ends (ratio 0.8370) and more distinct positions,
so an end-count ratio mostly tracks library depth. A call requires at
least 100 ends and 10 distinct positions on each side,
and HII/NT ≥ 2.0 on positions. NT sites = 0 is class `HII_only`
and has no ratio. The pseudocount (HII+1)/(NT+1) is not used.

Gene models: `sacCer3.sgd.gff.gz` (7154 gene-like features). SGD
assembly header R64-4-1; chromosome lengths match sacCer3 / R64-1-1,
including chrM = chrmt = 85,779 bp. Gene spans are the February 2024
transcript-extended parent features, not CDS-only intervals. SGD
`chrmt` is matched to pipeline `chrM`.

## Library totals

- NT ends: 9,640,282; of which no SGD gene hit: 2,243,992
- HII ends: 8,068,542; of which no SGD gene hit: 2,855,609
- Genome-wide HII/NT end ratio: 0.8370
- NT hotspot-window mass: 4,086,695
- HII hotspot-window mass: 4,568,073

Gene end sums double-count a position that falls in two gene spans.
Intergenic is the complement: ends whose position hits no gene.

## Where HII ends sit

| Compartment | Genes | HII ends | Share of library |
|---|---:|---:|---:|
| Mitochondrial | 55 | 1,101,434 | 13.7% |
| rRNA / RDN | 22 | 318,939 | 4.0% |
| Transposable element | 91 | 118,478 | 1.5% |
| Dubious ORF | 639 | 339,041 | 4.2% |
| End-proximal nuclear ORF | 252 | 295,571 | 3.7% |
| Nuclear interior ORF | 5694 | 4,428,314 | 54.9% |

Mitochondrial genes, rDNA, and subtelomeres can dominate raw end totals.
Incorporation rank below is nuclear interior ORFs (biotype `gene`,
not dubious, not within 20 kb of a chromosome end).

## Preferential incorporation (nuclear interior ORFs)

Primary rank: distinct HydEn HII positions per kb. That is where embedded
rNMPs are spread, not where one site was sequenced deeply. A second table
ranks HII ends per kb so pileups stay visible. Neither table is a human
disease rank.

### Distinct positions per kb (top 25)

| Symbol | Systematic | Length | HII sites/kb | HII sites | NT sites | Site HII/NT | End HII/NT |
|---|---|---:|---:|---:|---:|---:|---:|
| YLR154C-G | YLR154C-G | 150 | 986.6667 | 148 | 148 | 1.0 | 0.6413 |
| TAR1 | YLR154W-C | 375 | 984.0 | 369 | 360 | 1.025 | 0.6587 |
| FLO11 | YIR019C | 4104 | 845.5166 | 3470 | 2744 | 1.2646 | 1.6221 |
| BSC1 | YDL037C | 987 | 712.2594 | 703 | 520 | 1.3519 | 0.864 |
| YLR154C-H | YLR154C-H | 132 | 689.3939 | 91 | 89 | 1.0225 | 0.8348 |
| YBL081W | YBL081W | 1793 | 624.6514 | 1120 | 592 | 1.8919 | 1.5271 |
| YHR054C-B | YHR054C-B | 174 | 620.6897 | 108 | 96 | 1.125 | 0.6208 |
| URA3 | YEL021W | 1051 | 619.4101 | 651 | 329 | 1.9787 | 0.7758 |
| CUP1-2 | YHR055C | 186 | 612.9032 | 114 | 101 | 1.1287 | 0.5845 |
| DAN4 | YJR151C | 3486 | 575.7315 | 2007 | 1525 | 1.3161 | 0.9029 |
| AGA1 | YNR044W | 2403 | 575.1144 | 1382 | 809 | 1.7083 | 0.9982 |
| WSC3 | YOL105C | 2118 | 572.7101 | 1213 | 770 | 1.5753 | 0.9341 |
| WSC4 | YHL028W | 1984 | 566.0282 | 1123 | 688 | 1.6323 | 0.7713 |
| ZDS1 | YMR273C | 3062 | 559.7649 | 1714 | 650 | 2.6369 | 1.3214 |
| AIP5 | YFR016C | 3847 | 554.718 | 2134 | 1382 | 1.5441 | 0.7378 |
| STP22 | YCL008C | 1377 | 550.472 | 758 | 339 | 2.236 | 1.2418 |
| YVH1 | YIR026C | 1415 | 549.1166 | 777 | 373 | 2.0831 | 1.1646 |
| SPB1 | YCL054W | 2687 | 548.195 | 1473 | 641 | 2.298 | 1.4184 |
| RPS1B | YML063W | 1269 | 545.3113 | 692 | 428 | 1.6168 | 0.583 |
| MND2 | YIR025W | 1422 | 544.3038 | 774 | 418 | 1.8517 | 1.0837 |
| IST3 | YIR005W | 762 | 543.3071 | 414 | 233 | 1.7768 | 0.6933 |
| GLC8 | YMR311C | 1051 | 538.5347 | 566 | 239 | 2.3682 | 1.0298 |
| PXR1 | YGR280C | 1660 | 536.747 | 891 | 546 | 1.6319 | 0.6714 |
| YPT1 | YFL038C | 1022 | 536.2035 | 548 | 368 | 1.4891 | 0.4951 |
| DPM1 | YPR183W | 1423 | 534.0829 | 760 | 407 | 1.8673 | 0.7935 |

### End depth per kb (top 15)

High ends/kb with a site ratio near 1 is a pileup, not more incorporation sites.

| Symbol | Systematic | HII ends/kb | HII ends | End HII/NT | Site HII/NT |
|---|---|---:|---:|---:|---:|
| FLO11 | YIR019C | 12766.3255 | 52,393 | 1.6221 | 1.2646 |
| YLR154C-G | YLR154C-G | 12060.0 | 1,809 | 0.6413 | 1.0 |
| TAR1 | YLR154W-C | 11608.0 | 4,353 | 0.6587 | 1.025 |
| DAN4 | YJR151C | 7942.9145 | 27,689 | 0.9029 | 1.3161 |
| YLR154C-H | YLR154C-H | 5090.9091 | 672 | 0.8348 | 1.0225 |
| BSC1 | YDL037C | 4854.1033 | 4,791 | 0.864 | 1.3519 |
| NOP6 | YDL213C | 4579.3539 | 6,521 | 0.2925 | 1.1256 |
| NUP49 | YGL172W | 3529.4864 | 5,566 | 0.1088 | 0.9692 |
| WSC3 | YOL105C | 3447.1199 | 7,301 | 0.9341 | 1.5753 |
| YBL081W | YBL081W | 3333.5192 | 5,977 | 1.5271 | 1.8919 |
| URA3 | YEL021W | 3332.0647 | 3,502 | 0.7758 | 1.9787 |
| MTL1 | YGR023W | 2924.2783 | 6,179 | 1.1792 | 1.6765 |
| YPT1 | YFL038C | 2703.5225 | 2,763 | 0.4951 | 1.4891 |
| AGA1 | YNR044W | 2481.4815 | 5,963 | 0.9982 | 1.7083 |
| AHT1 | YHR093W | 2369.7632 | 1,301 | 0.6689 | 1.7554 |

## Failed excision classes

Class uses the unique-position ratio. End HII/NT is stored and is not the call.

| Class | Genes | Meaning |
|---|---:|---|
| HII_enriched | 1884 | both sides clear the floor and site HII/NT ≥ 2.0 |
| shared | 4599 | both sides clear the floor and 0.5 < site ratio < 2.0 |
| NT_enriched | 0 | both sides clear the floor and site HII/NT ≤ 0.5 |
| sparse | 664 | both sides nonzero, floor not met; ratio stored, not called |
| HII_only | 3 | NT sites or NT ends = 0; not a ratio |
| NT_only | 0 | HII sites or HII ends = 0; not a ratio |
| absent | 4 | no HydEn ends on the gene span |

Both sides nonzero: 7147 genes. Supported calls: 6483.
Nuclear interior supported calls: 5554, of which 1696 are HII_enriched.

### Nuclear HII-enriched (top 25 by site ratio)

More distinct positions after RNase HII cleavage than in the alkali-only
library, on a nuclear interior ORF. In rnh201 cells that is the
within-assay reading of embedded rNMPs the alkali library under-called.
It is not an excision-rate measurement.

| Symbol | Systematic | Length | Site HII/NT | HII sites | NT sites | End HII/NT |
|---|---|---:|---:|---:|---:|---:|
| YIP4 | YGL198W | 1118 | 3.6 | 198 | 55 | 1.0514 |
| TSR4 | YOL022C | 1432 | 3.5185 | 380 | 108 | 1.2896 |
| ILV6 | YCL009C | 1153 | 3.4907 | 377 | 108 | 1.343 |
| MEH1 | YKR007W | 701 | 3.4746 | 205 | 59 | 1.1003 |
| HEM12 | YDR047W | 1710 | 3.3056 | 357 | 108 | 1.2903 |
| KIN3 | YAR018C | 1559 | 3.2671 | 526 | 161 | 1.3609 |
| MIM1 | YOL026C | 667 | 3.2373 | 191 | 59 | 1.3346 |
| ATG19 | YOL082W | 1700 | 3.237 | 437 | 135 | 1.1437 |
| AIM45 | YPR004C | 1229 | 3.2368 | 246 | 76 | 1.0349 |
| YEL009C-A | YEL009C-A | 408 | 3.2188 | 103 | 32 | 1.7459 |
| ATG34 | YOL083W | 1555 | 3.1825 | 436 | 137 | 1.2232 |
| QCR10 | YHR001W-A | 1087 | 3.1688 | 244 | 77 | 1.0129 |
| VTC4 | YJL012C | 2372 | 3.157 | 543 | 172 | 1.0618 |
| TIR1 | YER011W | 1368 | 3.1513 | 375 | 119 | 1.2645 |
| XKS1 | YGR194C | 1803 | 3.1429 | 418 | 133 | 1.0767 |
| ERP1 | YAR002C-A | 986 | 3.141 | 245 | 78 | 1.0693 |
| RPE1 | YJL121C | 1239 | 3.1389 | 226 | 72 | 0.9095 |
| SFB3 | YHR098C | 3049 | 3.1381 | 568 | 181 | 1.069 |
| YOL085C | YOL085C | 342 | 3.1212 | 103 | 33 | 1.0602 |
| YLR111W | YLR111W | 621 | 3.1207 | 181 | 58 | 1.257 |
| ROG3 | YFR022W | 2632 | 3.1075 | 867 | 279 | 1.3068 |
| PDX1 | YGR193C | 2592 | 3.1027 | 695 | 224 | 1.3909 |
| HBN1 | YCL026C-B | 770 | 3.0938 | 297 | 96 | 1.2055 |
| UTP23 | YOR004W | 1034 | 3.0778 | 277 | 90 | 1.02 |
| SFC1 | YJR095W | 1259 | 3.069 | 356 | 116 | 0.9678 |

### Nuclear NT-enriched on positions

Site HII/NT ≤ 0.5 would mean the HII library lost positions relative to
alkali. That class is empty or short here because the HII library
expands unique positions genome-wide. End-count NT bias (NUP49, LYS1,
MEP2) is a depth pileup: those genes keep a site ratio near 1. Do not
read either pattern as successful excision.

No nuclear interior ORF clears the floor with site HII/NT ≤ 0.5.

## What this does not say

- No FTD, MAPT, or tau conclusion. Yeast has no MAPT.
- Hydroxyurea is Future work / v2 (Meroni et al. 2019, PMID 30847483), not a result ranked here.
- No Koh genotype contrast. Those libraries are not inputs.
- HII_only genes, and any NT=0 pseudocount ratio, are not failed-excision calls.
- An end-count ratio above 2 on mtDNA intron ORFs is depth, not the position rank.

Tables: `docs/gene-excision/gene_rnmp_sgd.tsv`,
`docs/gene-excision/incorporation_rank.tsv`,
`docs/gene-excision/failed_excision_supported.tsv`.
