BioJava 5.0.0
=============

unreleased

- For short structure selections (e.g. 1abc.A:1-100), ligands within 5A will be included


BioJava 4.2.5 
=============

release date: December 7th 2016
This release contains [30](https://github.com/biojava/biojava/compare/biojava-4.2.4...biojava-4.2.5) commits from 7 contributors.

### Bug fixes

* Fix for new phosphositeplus.org format, #610
* org.biojava.nbio.genome.parsers.gff.Location union() and intersect() now work correctly, #355
* Minor addition of crystallographic metadata fields to handle legacy PDB entries
* Jmol interchange format is now mmCIF, allowing for multiletter chain ids
* Update to latest jmol 14.6.2_2016.08.28
* A few minor bug fixes

BioJava 4.2.4
=============

release date: July 29th 2016
This release contains over [17](https://github.com/biojava/biojava/compare/biojava-4.2.3...biojava-4.2.4) commits from 4 contributors.

### Bug fixes

* NCBI links now using https (see [NCBI's announcement](http://www.ncbi.nlm.nih.gov/news/06-10-2016-ncbi-https/) )
* CATH links redirected to new server http://release.cathdb.info/
* SCOP default location now points to the Berkeley server after demise of Scop at MRC LMB
* Fixed important bug in mmCIF writing where structures with multiple models were written with identical coordinates
* Fixed bug in Group cloning where chemical components weren't cloned
* Added utility class for Chromosome mapping

BioJava 4.2.3
=============

release date: July 28th 2016
This release contains over [13](https://github.com/biojava/biojava/compare/biojava-4.2.2...biojava-4.2.3) commits from 2 contributors.

### Bug fixes

* mmCIF file writing: special fields (e.g. containing hyphens) are now correctly written
* General improvements in mmCIF file read and write

BioJava 4.2.2
=============

release date: June 14th 2016
This release contains over [31](https://github.com/biojava/biojava/compare/biojava-4.2.1...biojava-4.2.2) commits from 5 contributors.

### Bug fixes

This is a bug-fix release

* CE-Symm features and bug fixes
 - Better data structures for symmetry axes (particularly for hierarchical symmetry)
 - Fix bug with symmetry axis positioning
 - Optimization includes all symmetry repeats for hierarchical symmetry
* Update of protein modifications to latest version,
  - including new glycans and chromophores
  - Updating naming definitions to latest conventions


BioJava 4.2.1
=============

release date: May 3rd 2016
This release contains over [31](https://github.com/biojava/biojava/compare/biojava-4.2.0...biojava-4.2.1) commits from 7 contributors.

### Bug fixes

Biojava-structure

- Nucleotide bonds are now generated
- BIO: identifiers are now correctly handled
- Several fixes for CE-Symm
- Substructures now contain seqres groups (isse #449)
- Structures containing insertion codes are now written correctly to mmCIF
- AtomCache now uses the correct default parsing parameters (issue #455)
- Fixed problem with some atom charges that weren't being added 
- CATH updated to 4.0.0
- Better ECOD javadocs (issue #452)

Biojava-structure-gui
- Removed javaws dependency (issue #459)

BioJava 4.2.0
=============

release date: March 10th 2016

This release contains over [750](https://github.com/biojava/biojava/compare/6f8d796fee92edbbcd001c33cdae4f15c5480741...biojava-4.2.0) commits from 16 contributors.

### New Features

BioJava 4.2.0 offers many new features, as well several bug-fixes.

General  

-   Requires Java 7
-   Better logging with SLF4J

Biojava-Core  

-   New SearchIO framework including blast xml parser

Biojava-structure  

-   Secondary structure assignment (DSSP compatible)
-   Multiple Structure Alignments
    -   New MultipleStructureAlignment datastructure supporting flexible and order-independent alignments
    -   MultipleMC algorithm
        -   Can use any pairwise StructureAlignment implementation
    -   serialize and parse multiple structure alignments as XML files, output as Text, FatCat, FASTA, Rotation Matrices, etc.
-   More complete mmCIF and cif parsing
    -   Parse bonds, sites, charges
    -   Better support for non-deposited pdb and mmcif files
-   Include CE-Symm algorithm for finding internal symmetry (Myers-Turnbull, 2014)
-   Replaced internal graph datastructures with Jgraph
-   Unified StructureIdentifier framework
-   Improved chemical component framework, now by default providing full chemical description by using DownloadChemCompProvider
-   Optimised memory usage of Residue/Atoms

Biojava-structure-gui  

-   MultipleAlignmentGUI for visualizing Multiple Structure Alignments with Jmol
-   SymmetryDisplay for visualizing internal symmetry

Biojava-Phylo  

-   Use Forester 1.038
-   Significant bug fixes
-   use SubstitutionMatrices in the core module (instead of imported Jalview matrices),
-   use Sequence and Compound classes from the alignment module
-   provide some Wrapper methods to communicate with forester,
-   decouple distance matrix calculation from tree constructor,
-   provide methods for common distance matrix calculations and framework for user-defined distances,
-   update the forester version to have the correct NJ tree constructor
    AND
-   correct some of the tree evaluator statistics.









