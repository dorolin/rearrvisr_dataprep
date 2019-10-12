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
    
    &rarr; *Makes temporary files MEL_toRmTest.txt and MEL_info_test.txt.*
    
* Test whether gene positions overlap with genome positions of other genes when alternative transcripts are excluded, using min/max positions and taking strand into account.

    ```bash
    Rscript scripts/testGeneOverlaps.R MEL_info_test.txt MEL_overlapping_test_any.txt
    ```
    
    &rarr; *Returns file MEL_overlapping_test_any.txt (contains many overlapping genes, particularly for MEL and PSE).*

* Test whether gene positions overlap with genome positions of other genes when alternative transcripts are excluded, using midpoint positions and not taking strand into account.

    ```bash
    Rscript scripts/makeGeneListsOC.R MEL_info_test.txt MEL
    ```
    
    &rarr; *Returns files MEL_genome.txt and MEL_overlapping.txt (\*_overlapping.txt contains one to two overlapping genes for ANA, ERE, GRI, PER, PSE, and YAK; the file is not created when there are no overlapping genes).*
    
* Figure out more about the overlapping genes, if any (not applicable to MEL).

    ```bash
    grep -f <(cat ANA_overlapping.txt | tr ' ' '\n') ANA_genome.txt
    ```
    
    &rarr; *Overlapping genes are all on opposite strands except for PSE, where overlapping sequences seem to be duplicates (FBpp0302555 and FBpp0284283). **The duplicate FBpp0302555 was removed from all PSE\* files manually prior to further analysis.***
    
* Remove all test files that were created above.

    ```bash
    rm -f *toRmTest.txt *info_test.txt *overlapping_test_any.txt *overlapping.txt *genome.txt
    ```

* Prepare alternative splice variant files from MEL_info_transcripts.txt (remove gene names in first column, separate transcripts by `; `, and only keep rows with >1 transcript).

    ```bash
    cat MEL_info_transcripts.txt | sed 's/.*: \(.*$\)/\1/g' | sed 's/ /; /g' | grep ";"  >MEL.splice
    ```
    
    &rarr; *Makes file MEL.splice.*


## Identification of orthologs

Orthologs were identified for retained sequences of all 12 species with OMA standalone v2.2.0 (Altenhoff et al. 2015), using as guidance tree the phylogeny published in Drosophila 12 Genomes Consortium (2007), and default settings otherwise.

* Species tree in Newick format as guidance tree for OMA, which was added to the parameters/parameters.drw file.
    
    `SpeciesTree := '((((((MEL,(SIM,SEC)),(YAK,ERE)),ANA),(PSE,PER)),WIL),((MOJ,VIR),GRI));';`
    
* Set up OMA database in the directory where OMA will be run (path/to/OMA/).

    ```bash
    mkdir DB
    cp *_good.fasta *.splice DB/
    cd DB
    ## change file name extension and replace '|' delimiter in sequence header by '_'
    rename _good.fasta .fa ./*fasta
    sed -i 's/|/_/g' *.fa
    ## add species_ID to splice variant files
    for file in *.splice; do
        id="${file%.*}"
        sed -i "s/^/${id}_/g" ${id}.splice
        sed -i "s/; /; ${id}_/g" ${id}.splice
    done
    cd ..
    ```

* Place the parameters/parameters.drw file in the directory where OMA will be run (path/to/OMA/) and run OMA in three steps. More information on running OMA and usage examples (e.g., for parallelization) are in the [OMA standalone documentation](https://omabrowser.org/standalone/#downloads).
    
    ```bash
    OMA -c 
    OMA -s
    OMA
    ```

    &rarr; *Creates the directory Output/ that contains files that will be used in the steps below.*

## Phylogeny

A species tree with branch lengths is required as input file for the ancestral genome reconstruction software ANGES v1.01 (Jones et al. 2012). To generate such a tree, sequences of OMA orthologous groups that only included one-to-one orthologous genes present in all 12 *Drosophila* species were extracted and aligned separately for each orthologous group with MAFFT v7.407 (Katoh et al. 2002; Katoh and Standley 2013). Alignments with not more than 20% missing data were concatenated, and a phylogenetic tree was computed with RAxML v8.2.12 (Stamatakis 2014).

* Identify orthologous groups in OMA output that have genes for all species. This is based on data in the path/to/OMA/Output/PhyleticProfileOMAGroups.txt file generated by OMA.

    ```bash
    Rscript scripts/getOGs.R path/to/OMA/Output
    ```
    
    &rarr; *Creates file path/to/OMA/Output/fullOGs.txt.*

* For the orthologous groups defined in fullOGs.txt, extract orthologs with one-to-one relationship in all species pairs. One-to-one orthologs are identified based on data in the path/to/OMA/Output/PairwiseOrthologs/ directory generated by OMA. *Note:* the bash script below is not very efficient and takes a long time to complete.

    ```bash
    scripts/getOneOne.sh path/to/OMA/Output
    ```
    
    &rarr; *Creates file path/to/OMA/Output/fullOneOneOGfiles.txt.*
    
* Copy the fasta files listed in path/to/OMA/Output/fullOneOneOGfiles.txt and that are located in the path/to/OMA/Output/OrthologousGroupsFasta/ directory to the directory where MAFFT will be run (path/to/MAFFT/).

    ```bash
    while read file; do
        cp path/to/OMA/Output/OrthologousGroupsFasta/${file} path/to/MAFFT/
    done < path/to/OMA/Output/fullOneOneOGfiles.txt
    ```
    
* Align sequences for each orthologous group with MAFFT, using the iterative refinement method incorporating local pairwise alignment information, as shown for one example below. More information on running MAFFT are in the [MAFFT manual](https://mafft.cbrc.jp/alignment/software/manual/manual.html).

    ```bash
    mafft --localpair --maxiterate 1000 path/to/MAFFT/OG10003.fa >path/to/MAFFT/OG10003.msa
    ```
    
    &rarr; *Returns file path/to/MAFFT/OG10003.msa.*
    
* Convert all alignments to phylip format using the bash script [convertFasta2Phylip.sh](https://github.com/stamatak/standard-RAxML/blob/master/usefulScripts/convertFasta2Phylip.sh) distributed with RAxML, as shown for one example below.

    ```bash
    convertFasta2Phylip.sh path/to/MAFFT/OG10003.msa >path/to/MAFFT/OG10003.phy
    ```

    &rarr; *Returns file path/to/MAFFT/OG10003.phy.*
    
* Obtain information for all alignments in phylip format that are in the path/to/MAFFT/ directory.

    ```bash
    scripts/getAlignmentInfo.sh path/to/MAFFT
    ```
    
    &rarr; *Creates file path/to/MAFFT/alignInfo.txt.*
    
* Determine alignments with <=20% missing data (where missing data are `-` or `X` characters in the alignment).

    ```bash
    Rscript scripts/filterAlignments.R path/to/MAFFT
    ```
    
    &rarr; *Creates file path/to/MAFFT/filteredAlignments.txt.*
    
* Concatenate alignments for kept orthologous groups listed in path/to/MAFFT/filteredAlignments.txt.

    ```bash
    scripts/concatAlignments.sh path/to/MAFFT
    ```
    
    &rarr; *Returns file path/to/MAFFT/concatOGAlignment.phy.*
    
* 



## Preparation of genome maps


## Ancestral genome reconstruction


## Identification of rearrangements


# References

Altenhoff AM, Skunca N, Glover N, et al. (2015) The OMA orthology database in 2015: function predictions, better plant support, synteny view and other improvements. *Nucleic Acids Research*, **43**, D240–D249.

Drosophila 12 Genomes Consortium (2007) Evolution of genes and genomes on the *Drosophila* phylogeny. *Nature*, **450**, 203–218.

Jones BR, Rajaraman A, Tannier E, Chauve C (2012) ANGES: reconstructing ANcestral GEnomeS maps. *Bioinformatics*, **28**, 2388–2390.

Katoh K, Misawa K, Kuma K, Miyata T (2002) MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. *Nucleic Acids Research*, **30**, 3059–3066.

Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*, **30**, 772–780.

Li L, Stoeckert Jr. CJ, Roos DS (2003) OrthoMCL: identification of ortholog groups for eukaryotic genomes. *Genome Research*, **13**, 2178–2189.

Stamatakis A (2014) RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. *Bioinformatics*, **30**, 1312–1313.

