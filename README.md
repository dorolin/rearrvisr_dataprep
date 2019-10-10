# rearrvisr_dataprep

*`[under construction]`*

The scripts and BASH code snippets in this repository were used to prepare example input files for the R package [`rearrvisr`](https://github.com/dorolin/rearrvisr). Some of them are more sophisticated Perl scripts that have multiple command line arguments and can be adjusted to work on a variety of input files. Others are rather unsophisticated ("quick-and-dirty") BASH code snippets that might easily fail on non-Linux operation systems, and where a lot of improvement in efficiency could be done. There is absolutely no warranty that the workflow presented below will produce accurate results when repeated on fasta or gff3 files that are formatted differently than the ones for which the scripts were written.

## Genome assemblies

Peptide sequences and annotation information (pep.all.fa and gff3 files) for genes from 12 *Drosophila* species were downloaded on Dec 23 2017 from Ensemble Release 91 ([http://dec2017.archive.ensembl.org](http://dec2017.archive.ensembl.org); *D. melanogaster*) or Ensemble Metazoa Release 37 ([http://oct2017-metazoa.ensembl.org](http://oct2017-metazoa.ensembl.org); *D. ananassae*, *D. erecta*, *D. grimshawi*, *D. mojavensis*, *D. persimilis*, *D. pseudoobscura*, *D. sechellia*, *D. simulans*, *D. virilis*, *D. willistoni*, and *D. yakuba*).

Species | ID | Genome version | Source | Date
:-------|:---|:---------------|:-------|:----
*D. ananassae*     | ANA | dana_caf1 | Ensemble Metazoa 37 | 15/08/17 
*D. erecta*        | ERE | dere_caf1 | Ensemble Metazoa 37 | 15/08/17 
*D. grimshawi*     | GRI | dgri_caf1 | Ensemble Metazoa 37 | 15/08/17 
*D. melanogaster*  | MEL | BDGP6     | Ensemble Release 91 | 22/11/17 
*D. mojavensis*    | MOJ | dmoj_caf1 | Ensemble Metazoa 37 | 15/08/17 
*D. persimilis*    | PER | dper_caf1 | Ensemble Metazoa 37 | 15/08/17 
*D. pseudoobscura* | PSE | Dpse_3.0  | Ensemble Metazoa 37 | 15/08/17 
*D. sechellia*     | SEC | dsec_caf1 | Ensemble Metazoa 37 | 15/08/17 
*D. simulans*      | SIM | dsim_caf1 | Ensemble Metazoa 37 | 15/08/17 
*D. virilis*       | VIR | dvir_caf1 | Ensemble Metazoa 37 | 15/08/17 
*D. willistoni*    | WIL | dwil_caf1 | Ensemble Metazoa 37 | 15/08/17 
*D. yakuba*        | YAK | dyak_caf1 | Ensemble Metazoa 37 | 15/08/17 

## Preparation of FASTA files

Fasta files for all 12 species were prepared as shown below for *D. melanogaster*.

* Modify header lines in the peptide fasta files to contain only species_ID and peptide_ID, separated by `|` (e.g., `>MEL|FBpp0075257`), and remove trailing `*` (if any) from sequences, as they are not allowed in OMA standalone v2.2.0 (Altenhoff et al. 2015) input files.

    ```bash
    scripts/fastaAdjust.pl Drosophila_melanogaster.BDGP6.pep.all.fa MEL " " 1
    ```
    This script is based on code from the orthoMCL `orthomclAdjustFasta` tool from the OrthoMCL software (Li et al. 2003) [https://orthomcl.org/common/downloads/software/v2.0/].

    &rarr; *Returns file MEL.fasta.*

* Exclude sequences with intermediate `*` character (i.e., premature stop codons) or that are shorter than 50 amino acids.

    ```bash
    scripts/fastaFilter.pl MEL.fasta 50 0
    ```
    This script is based on code from the orthoMCL `orthomclFilterFasta` tool from the OrthoMCL software (Li et al. 2003) [https://orthomcl.org/common/downloads/software/v2.0/].
    
    &rarr; *Returns files MEL_good.fasta and MEL_poor.fasta.*
    
* Compute the lengths of retained sequences.

    ```bash
    scripts/fastaInfo.pl MEL_good.fasta 2
    ```
    
    &rarr; *Returns file MEL_good_info.txt.*
    
## Preparation of GFF3  files

Gff3 files for all 12 species were prepared as shown below for *D. melanogaster*.

* For retained sequences, extract their genome position from gff3 files, and generate a list with alternative transcripts per gene.

    ```bash
    scripts/gff3Info3.pl -l MEL_good_info.txt -c 2 -g Drosophila_melanogaster.BDGP6.91.gff3 \
                         -n protein_id -m gene_id -o MEL_info.txt
    ```
    
    &rarr; *Returns files MEL_info.txt and MEL_info_transcripts.txt.*

* Verify that all genes have gene info.

    ```bash
    grep 'NA' MEL_info.txt
    ```
    
    &rarr; *Nothing should be returned. Requires trouble-shooting otherwise.*

* Temporarily remove alternative transcripts. *Note:* the bash commands below are not very efficient and take a long time to complete with many alternative transcripts.

    ```bash
    cut -d ' ' -f3- MEL_info_transcripts.txt | sed '/^\s*$/d' | tr ' ' '\n' >MEL_toRmTest.txt
    cat MEL_info.txt | grep -v -w -f MEL_toRmTest.txt >MEL_info_test.txt
    ```
    
    &rarr; *Make temporary files MEL_toRmTest.txt and MEL_info_test.txt.*
    
* Test whether gene positions overlap with genome positions of other genes when alternative transcripts are excluded, using min/max positions and taking strand into account.

    ```bash
    Rscript scripts/testGeneOverlaps.R MEL_info_test.txt MEL_overlapping_test_any.txt
    ```
    
    &rarr; *Returns file MEL_overlapping_test_any.txt (contains many overlapping genes, particularly for MEL and PSE).*

* Test whether gene positions overlap with genome positions of other genes when alternative transcripts are excluded, using midpoint positions and not taking strand into account.

    ```bash
    Rscript scripts/makeGeneListsOC.R MEL_info_test.txt MEL
    ```
    
    &rarr; *Returns files MEL_genome.txt and MEL_overlapping.txt (contains one to two overlapping genes for ANA, ERE, GRI, PER, PSE, and YAK).*
    
* Figure out more about the overlapping genes, if any (not applicable to MEL).

    ```bash
    grep -f <(cat ANA_overlapping.txt | tr ' ' '\n') ANA_genome.txt
    ```
    
    &rarr; *Overlapping genes are all on opposite strands except for PSE, where overlapping sequences seem to be duplicates (FBpp0302555 and FBpp0284283).*
    
* Remove all test files that were created above.

    ```bash
    rm -f *toRmTest.txt *info_test.txt *overlapping_test_any.txt *overlapping.txt *genome.txt
    ```
    
## Identification of orthologs


## Phylogeny


## Preparation of genome maps


## Ancestral genome reconstruction


## Identification of rearrangements


# References

Altenhoff,A.M. et al. (2015) The OMA orthology database in 2015: function predictions, better plant support, synteny view and other improvements. *Nucleic Acids Research*, **43**, D240–D249.

Li L, Stoeckert Jr. CJ, Roos DS (2003) OrthoMCL: identification of ortholog groups for eukaryotic genomes. *Genome Research*, **13**, 2178–2189.
