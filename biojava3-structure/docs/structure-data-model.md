# The BioJava-structure data model

A biologically and chemically meaningful data representation of PDB/mmCIF.

## The basics   

BioJava at its core is a collection of file parsers and (in some cases) data models to represent frequently used biological data.  The protein-structure modules represent macromolecular data in a way that should make it easy to work with. The representation is essentially independ of the underlying file format and the user can chose to work with either PDB or mmCIF files and still get an almost identical data representation.

## The main hierarchy

BioJava provides a flexible data structure for managing protein structural data. The 
[http://www.biojava.org/docs/api/org/biojava/bio/structure/Structure.html Structure] class is the main container. 

A Structure has a hierarchy of sub-objects:

<pre>
Structure
   |
   Model(s)
        |
        Chain(s)
            |
             Group(s)
                 |
                 Atom(s)
</pre>

All structure objects contain one or more "models". That means also X-ray structures contain an "virtual" model which serves as a container for the chains. The most common way to access chains will be via

<pre>
        List<Chain>chains = structure.getChains();
</pre>

This works for both NMR and X-ray based structures and by default the first model is getting accessed.


## Working with atoms

Different ways are provided how to access the data contained in a [Structure](http://www.biojava.org/docs/api/org/biojava/bio/structure/Structure.html).
If you want to directly access an array of [Atoms](http://www.biojava.org/docs/api/org/biojava/bio/structure/Atom.html) you can use the utility class called [StructureTools](http://www.biojava.org/docs/api/org/biojava/bio/structure/StructureTools.html)

<pre>

    // get all C-alpha atoms in the structure
    Atom[] caAtoms = StructureTools.getAtomCAArray(structure);
</pre>


## Working with groups

The [Group](http://www.biojava.org/docs/api/org/biojava/bio/structure/Group.html) interface defines all methods common to a group of atoms. There are 3 types of Groups:

* [AminoAcid](http://www.biojava.org/docs/api/org/biojava/bio/structure/AminoAcid.html)
* [Nucleotide](http://www.biojava.org/docs/api/org/biojava/bio/structure/NucleotideImpl.html) 
* [Hetatom](http://www.biojava.org/docs/api/org/biojava/bio/structure/HetatomImpl.html) 

In order to get all amino acids that have been observed in a PDB chain, you can use the following utility method:

<pre>
            Chain chain = s.getChainByPDB("A");
            List<Group> groups = chain.getAtomGroups("amino");
            for (Group group : groups) {
                AminoAcid aa = (AminoAcid) group;

                // do something amino acid specific, e.g. print the secondary structure assignment
                System.out.println(aa + " " + aa.getSecStruc());
            }
</pre>


In a similar way you can access all nucleotide groups by
<pre>
            chain.getAtomGroups("nucleotide");
</pre>

The Hetatom groups are access in a similar fashion:
<pre>
            chain.getAtomGroups("hetatm");
</pre>


Since all 3 types of groups are implementing the Group interface, you can also iterate over all groups and check for the instance type:

<pre>
            List<Group> allgroups = chain.getAtomGroups();
            for (Group group : groups) {
                if ( group instanceof AminoAcid) {
                    AminoAcid aa = (AminoAcid) group;
                    System.out.println(aa.getSecStruc());
                }
            }
</pre>


### How does BioJava decide what groups are amino acids?

BioJava supports the [Chemical Components Dictionary](http://www.wwpdb.org/ccd.html) which specifies the correct representation of each group. Let's take a look at how [Selenomethionine](http://en.wikipedia.org/wiki/Selenomethionine) and water is dealt with:

<pre>
            Structure structure = StructureIO.getStructure("1A62");
                    
            for (Chain chain : structure.getChains()){
                for (Group group : chain.getAtomGroups()){
                    if ( group.getPDBName().equals("MSE") || group.getPDBName().equals("HOH")){
                        System.out.println(group.getPDBName() + " is a group of type " + group.getType());
                    }
                }
            }
</pre>

This should give this output:

<pre>
MSE is a group of type amino
MSE is a group of type amino
MSE is a group of type amino
HOH is a group of type hetatm
HOH is a group of type hetatm
HOH is a group of type hetatm
...
</pre>

As you can see, although MSE is flaged as HETATM in the PDB file, BioJava still represents it correctly as an amino acid. They key is that the [definition file for MSE](http://www.rcsb.org/pdb/files/ligand/MSE.cif) flags it as "L-PEPTIDE LINKING", which is being used by BioJava.


### How to access Chemical Component definitions
Bye default BioJava ships with a minimal representation of standard amino acids, however if you want to parse the whole PDB archive, it is good to tell the library to either

1. fetch missing Chemical Component definitions on the fly (small download and parsing delays every time a new chemical compound is found), or
2. Load all definitions at startup (slow startup, but then no further delays later on)

You can enable the first behaviour by doing using the [FileParsingParameters](http://www.biojava.org/docs/api/org/biojava/bio/structure/io/FileParsingParameters.html) class:

<pre>
            AtomCache cache = new AtomCache();
            
             // by default all files are stored at a temporary location.
            // you can set this either via at startup with -DPDB_DIR=/path/to/files/
            // or hard code it this way:
            cache.setPath("/tmp/");
            
            FileParsingParameters params = new FileParsingParameters();
            
            params.setLoadChemCompInfo(true);
            cache.setFileParsingParams(params);
            
            StructureIO.setAtomCache(cache);
            
            Structure structure = StructureIO.getStructure(...);
</pre>

If you want to enable the second behaviour (slow loading of all chem comps at startup, but no further small delays later on) you can use the same code but change the behaviour by switching the [ChemCompProvider](http://www.biojava.org/docs/api/org/biojava/bio/structure/io/mmcif/ChemCompProvider.html) implementation in the [ChemCompGroupFactory](http://www.biojava.org/docs/api/org/biojava/bio/structure/io/mmcif/ChemCompGroupFactory.html)

<pre>
     ChemCompGroupFactory.setChemCompProvider(new AllChemCompProvider());
</pre>






