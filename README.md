# rearrvisr_dataprep

The scripts and code snippets in this repository were used to prepare example input files for the R package [`rearrvisr`](https://github.com/dorolin/rearrvisr). Some of them are more sophisticated Perl scripts that have multiple command line arguments and can be adjusted to work on a variety of input files. Others are rather unsophisticated ("quick-and-dirty") BASH code snippets that might easily fail on non-Linux operation systems, and where a lot of improvement in efficiency could be done. There is absolutely no warranty that the workflow presented below will produce accurate results when repeated on fasta or gff3 files that are formatted differently than the ones for which the scripts were written.


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

    &rarr; *Returns file MEL_genome.txt. (In case that genes with overlapping midpoints still exist, a second file \*_overlapping.txt is returned; this is the case for ANA, ERE, GRI, PER, PSE, and YAK, containing one to two overlapping genes per species).*
    
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

Orthologs were identified for retained sequences of all 12 species with OMA standalone v2.2.0 (Altenhoff et al. 2015), using as guidance tree the phylogeny published in Drosophila 12 Genomes Consortium (2007), and default settings otherwise. More information on running OMA and usage examples (e.g., for parallelization) are in the [OMA standalone documentation](https://omabrowser.org/standalone/#downloads).

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

* Place the parameters/parameters.drw file in the directory where OMA will be run (path/to/OMA/) and run OMA in three steps.
    
    ```bash
    OMA -c 
    OMA -s
    OMA
    ```

    &rarr; *Creates the directory Output/ that contains files that will be used in the steps below.*


## Phylogeny

A species tree with branch lengths is required as input file for the ancestral genome reconstruction software ANGES v1.01 (Jones et al. 2012). To generate such a tree, sequences of OMA orthologous groups that only included one-to-one orthologous genes present in all 12 *Drosophila* species were extracted and aligned separately for each orthologous group with MAFFT v7.407 (Katoh et al. 2002; Katoh and Standley 2013). More information on running MAFFT is in the [MAFFT manual](https://mafft.cbrc.jp/alignment/software/manual/manual.html). Alignments with no more than 20% missing data were concatenated, and a phylogenetic tree was computed with RAxML v8.2.12 (Stamatakis 2014). More information on running RAxML and usage examples are in the [RAxML manual](https://github.com/stamatak/standard-RAxML/blob/master/manual/NewManual.pdf) and the [RAxML hands-on session](https://cme.h-its.org/exelixis/web/software/raxml/hands_on.html).

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
    
* Align sequences for each orthologous group with MAFFT, using the iterative refinement method incorporating local pairwise alignment information, as shown for one example below. 

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
    
* Copy path/to/MAFFT/concatOGAlignment.phy to the directory where RAxML will be run (path/to/RAxML/), and change to this directory. Then remove columns from the alignment that only contain undetermined values by performing a test run with RAxML.

    ```bash
    raxmlHPC-AVX -f c -m PROTGAMMAAUTO -s concatOGAlignment.phy -n Test
    ```
    
    &rarr; *Returns file concatOGAlignment.phy.reduced.*

* Compute phylogenetic tree with RAxML, using the rapid bootstrapping algorithm and automatic determination of the best protein substitution model.

    ```bash
    raxmlHPC-PTHREADS-AVX -f a -x 6438990 -p 74378439 -# autoMRE -m PROTGAMMAAUTO \
                          --auto-prot=ml -s concatOGAlignment.phy.reduced -n One -T 24
    ```
    
    &rarr; *Returns unrooted phylogenetic tree RAxML_bestTree.One.*
    
* Root phylogenetic tree at the branch that best balances the subtree lengths.

    ```bash
    raxmlHPC-AVX -f I -m PROTGAMMAJTTF -t RAxML_bestTree.One -n OneRooted
    ```
    
    &rarr; *Returns rooted phylogenetic tree RAxML_rootedTree.OneRooted.*


## Preparation of genome maps

Genome maps for all retained genes (i.e., excluding low quality sequences and non-representative splicing variants) for all 12 species were created as shown below for *D. melanogaster* and *D. ananassae*.

* Identify alternative splice variants that were not used by OMA for indentifying orthologs and exclude them from the initial gene position files for each species. The used splice variants for all species are listed in the path/to/OMA/Output/used_splicing_variants.txt file generated by OMA. *Note:* the bash commands below are not very efficient and take a long time to complete with many alternative transcripts.

    ```bash
    ## make list with unused splice variants using the MEL.splice file created above
    cat path/to/OMA/Output/used_splicing_variants.txt | grep ^MEL | cut -f2 | sed "s/MEL_//g" >tmp
    cat tmp | while read line; do
        cat MEL.splice | grep $line | sed 's/ //g' | tr ';' '\n' | grep -v $line >>MEL_unused_splicing_variants.txt
    done
    rm -f tmp
    ## remove unused splice variants from the initial MEL_info.txt file created above
    cat MEL_info.txt | grep -v -w -f MEL_unused_splicing_variants.txt >MEL_info_1splice.txt
    ```
    
    &rarr; *Creates file MEL_info_1splice.txt.*

* Prepare genome maps with gene start and end positions calculated as CDS midpoints -/+ 1 base pair. This largely avoids the occurrence of overlapping gene positions in the genome maps, as overlaps are not supported by the genome reconstruction software ANGES v1.01 (Jones et al. 2012). The script below outputs a list of genes with overlapping midpoint positions (not taking strand into account) in case that overlaps still exist. 

    ```bash
    Rscript scripts/makeGeneListsOC.R MEL_info_1splice.txt MEL
    ```

    &rarr; *Returns file MEL_genome.txt. (In case that genes with overlapping midpoints still exist, a second file \*_overlapping.txt is returned; this is the case for ANA, ERE, GRI, PER, and YAK, containing one to two overlapping genes per species).*

* Figure out more about the overlapping genes, if any (not applicable to MEL).

    ```bash
    grep -f <(cat ANA_overlapping.txt | tr ' ' '\n') ANA_genome.txt
    ```
    
    &rarr; *Overlapping genes are all on opposite strands for ANA, ERE, GRI, PER, and YAK.*
    
* In case that overlapping gene positions still exist, genome maps need to be adjusted either manually or by running the script below that automatically shifts gene positions by a few base pairs to solve overlaps.

    ```bash
    Rscript scripts/makeGeneListsOC_beta.R ANA_info_1splice.txt ANA
    ```

    &rarr; *Returns file ANA_genome.txt.*

* Test whether species or chromosome names contain characters `.`, `:`, or `-`, which are not supported by ANGES. *Note:* some other special characters, like `@` or `|`, can also result in ANGES to crash without specific error messages.

    ```bash
    cut -d ' ' -f2 MEL_genome.txt | grep "\.\|:\|-" | head -n 1
    ```
    
    &rarr; *Nothing should be returned. Requires replacement of unsupported characters in the \*_genome.txt files otherwise.*


## Ancestral genome reconstruction

Input files for the genome reconstruction software ANGES v1.01 (Jones et al. 2012) were created based on the genome maps for all 12 *Drosophila* species and the phylogenetic tree prepared above. More information on running ANGES is in the [ANGES manual](https://github.com/cchauve/ANGeS/blob/master/anges_1_01_v2.pdf).

* Copy the genome maps for all 12 species (\*_genome.txt) to a new directory (path/to/ANGES/input/) that does not contain other files that match the pattern \*_genome.txt. Generate the 'Markers' input file for ANGES, as well as modified versions of \*_genome.txt (i.e., \*_markers.txt), which are required as input files for `rearrvisr`. This step is based on data in the path/to/OMA/Output/OrthologousGroups.txt file generated by OMA and the genome maps prepared above.

    ```bash
    Rscript scripts/oma2anges.R path/to/ANGES/input path/to/OMA/Output
    ```
    
    &rarr; *Returns the file Markers and files \*_markers.txt for each species.*

* Generate the 'Tree' input file for ANGES. Simply copy the rooted phylogenetic tree generated by RAxML (path/to/RAxML/RAxML_rootedTree.OneRooted) to path/to/ANGES/input/ and change the name of the file to Tree. Then use a text editor to insert an `@` character at the ancestral node for which a genome map will be reconstructed (i.e., for the ancestor of the *melanogaster* subgroup *'MSSYE'*, separating *D. melanogaster*, *D. simulans*, *D. sechellia*, *D. yakuba*, and *D. erecta* from the remainder of the *Drosophila* species). *Note:* The original tree had one branch of length 0.000000, which was edited to 0.000001.

    `(((GRI:0.107008,(MOJ:0.097938,VIR:0.058702):0.026270):0.121882,WIL:0.173794):0.000001,((ANA:0.102866,(((SEC:0.007967,SIM:0.006847):0.007368,MEL:0.012766):0.014760,(ERE:0.021357,YAK:0.019588):0.006552)@:0.077823):0.053434,(PSE:0.002489,PER:0.007924):0.114863):0.049467);`

* Create a new directory where ANGES will be run (path/to/ANGES/run1/) and change to this directory. Create the directory INPUT/ and copy or move the ../input/Markers and the ../input/Tree file to INPUT/, and also place the parameters/Parameters file in the INPUT/ directory. Then run ANGES by calling the master script anges_CAR.py of the software.

    ```bash
    python path/to/anges_1.01/src/MASTER/anges_CAR.py INPUT/Parameters &
    ```
    
    &rarr; *Creates the directory CARS/ that contains the file MSSYE_PQTREE_HEUR.*


## Identification of rearrangements

Genome rearrangements in *D. melanogaster* that occurred after divergence from the *melanogaster* subgroup ancestor *'MSSYE'* were identified and visualized with the R package [`rearrvisr`](https://github.com/dorolin/rearrvisr). More information on running `rearrvisr` and a description of the output are in the [package vignette](https://github.com/dorolin/rearrvisr/blob/master/vignettes/rearrvisr.pdf) and the function documentation.

* Read and verify input files in R.

    ```R
    library(rearrvisr)
    
    ## read markers file for D. melanogaster
    MEL_markers <- read.table(file = "path/to/ANGES/input/MEL_markers.txt", header = TRUE, as.is = TRUE)
    MEL_markers$scaff <- as.character(MEL_markers$scaff)
    ## verify file format
    checkInfile(MEL_markers, "focalgenome", checkorder = TRUE)
    
    ## read raw PQ-tree of reconstucted genome for 'MSSYE'
    MSSYE_PQTREE_HEUR <- read.table("path/to/ANGES/run1/CARS/MSSYE_PQTREE_HEUR", sep = ",", 
                                    comment.char = "", as.is = TRUE)  
    ## convert and verify file format
    MSSYE_compgenome <- convertPQtree(MSSYE_PQTREE_HEUR)
    checkInfile(MSSYE_compgenome, "compgenome", checkorder = TRUE)
    ```
    
* Identify rearrangements and summarize information for blocks of conserved marker order.

    ```R
    ## identify rearrangements (runs for some seconds)
    SYNT_MEL_MSSYE <- computeRearrs(MEL_markers, MSSYE_compgenome, doubled = TRUE)
    ## summarize blocks
    BLOCKS_MEL_MSSYE <- summarizeBlocks(SYNT_MEL_MSSYE, MEL_markers, MSSYE_compgenome,
                                        c("2L", "2R", "3L", "3R", "X"))
    ```
    
    &rarr; *Returns lists SYNT_MEL_MSSYE and BLOCKS_MEL_MSSYE that store rearrangements and additional information.*
    
* Make plots in pdf format to visualize rearrangements.

    ```R
    genomeImagePlot(SYNT_MEL_MSSYE, MEL_markers, c("2L", "2R", "3L", "3R", "X"),
                main = "D. melanogaster - MSSYE", makepdf = TRUE, filename = "MEL_genome.pdf")
    genomeRearrPlot(BLOCKS_MEL_MSSYE, MSSYE_compgenome, c("2L", "2R", "3L", "3R", "X"),
                main = "D. melanogaster - MSSYE", blockwidth = 1.15, y0pad = 3,
                makepdf = TRUE, filename = "MEL_rearr.pdf")
    ```
    
    &rarr; *Creates pdfs MEL_genome.pdf and MEL_rearr.pdf that show rearrangements along five D. melanogaster chromosomes.*



# References

Altenhoff AM, Skunca N, Glover N, et al. (2015) The OMA orthology database in 2015: function predictions, better plant support, synteny view and other improvements. *Nucleic Acids Research*, **43**, D240–D249.

Drosophila 12 Genomes Consortium (2007) Evolution of genes and genomes on the *Drosophila* phylogeny. *Nature*, **450**, 203–218.

Jones BR, Rajaraman A, Tannier E, Chauve C (2012) ANGES: reconstructing ANcestral GEnomeS maps. *Bioinformatics*, **28**, 2388–2390.

Katoh K, Misawa K, Kuma K, Miyata T (2002) MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. *Nucleic Acids Research*, **30**, 3059–3066.

Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*, **30**, 772–780.

Li L, Stoeckert Jr. CJ, Roos DS (2003) OrthoMCL: identification of ortholog groups for eukaryotic genomes. *Genome Research*, **13**, 2178–2189.

Stamatakis A (2014) RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. *Bioinformatics*, **30**, 1312–1313.

