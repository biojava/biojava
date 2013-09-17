# How to parse mmCIF files using BioJava

A quick tutorial how to work with mmCIF files.

## What is mmCIF?

The Protein Data Bank (PDB) has been distributing its archival files as PDB files for a long time. The PDB file format is based on "punchcard"-style rules how to store data in a flat file. With the increasing complexity of macromolecules that have are being resolved experimentally, this file format can not be used any more to represent some or the more complex structures. As such, the wwPDB recently announced the transition from PDB to mmCIF/PDBx as  the principal deposition and dissemination file format (see 
[here](http://www.wwpdb.org/news/news_2013.html#22-May-2013) and 
[here](http://wwpdb.org/workshop/wgroup.html)). 

The mmCIF file format has been around for some time (see [Westbrook 2000][] and [Westbrook 2003][] ) [BioJava](http://www.biojava.org) has been supporting mmCIF already for several years. This tutorial is meant to provide a quick introduction into how to parse mmCIF files using [BioJava](http://www.biojava.org)

## The basics

BioJava provides you with both a mmCIF parser and a data model that reads PDB and mmCIF files into a biological and chemically  meaningful data model (BioJava supports the [Chemical Components Dictionary](http://www.wwpdb.org/ccd.html)). If you don't want to use that data model, you can still use BioJava's file parsers, and more on that later, let's start first with the most basic way of loading a protein structure.

## Quick Installation

Before we start, just one quick paragraph of how to get access to BioJava.

BioJava is open source and you can get the code from [Github](https://github.com/biojava/biojava), however it might be easier this way:

BioJava uses [Maven](http://maven.apache.org/) as a build and distribution system. If you are new to Maven, take a look at the [Getting Started with Maven](http://maven.apache.org/guides/getting-started/index.html)  guide.

We are providing a BioJava specific Maven repository at (http://biojava.org/download/maven/) .

You can add the BioJava repository by adding the following XML to your project pom.xml file:
```xml
        <repositories>
            ...
            <repository>
                <id>biojava-maven-repo</id>
                <name>BioJava repository</name>
                <url>http://www.biojava.org/download/maven/</url>           
            </repository>
        </repositories>
        <dependencies>
                ...
                <dependency>
                        <!-- This imports the latest SNAPSHOT builds from the protein structure modules of BioJava
                        -->                        
                        <groupId>org.biojava</groupId>
                        <artifactId>biojava3-structure</artifactId>
                        <version>3.0.7-SNAPSHOT</version>
                </dependency>
                <!-- other biojava jars as needed -->
        </dependencies> 
```

If you run 'mvn package' on your project, the BioJava dependencies will be automatically downloaded and installed for you.

## First steps

The simplest way to load a PDB file is by using the [StructureIO](http://www.biojava.org/docs/api/org/biojava3/structure/StructureIO.html) class.

<pre>
    Structure structure = StructureIO.getStructure("4HHB");
    // and let's print out how many atoms are in this structure
    System.out.println(StructureTools.getNrAtoms(structure));
</pre>



BioJava  automatically downloaded the PDB file for hemoglobin [4HHB](http://www.rcsb.org/pdb/explore.do?structureId=4HHB) and copied it into a temporary location. This demonstrates two things:

+ BioJava can automatically download and install files locally
+ BioJava by default writes those files into a temporary location (The system temp directory "java.io.tempdir"). 

If you already have a local PDB installation, you can configure where BioJava should read the files from by setting the PDB_DIR system property

<pre>
    -DPDB_DIR=/wherever/you/want/
</pre>

## From PDB to mmCIF

By default BioJava is using the PDB file format for parsing data. In order to switch it to use mmCIF, we can take control over the underlying [AtomCache](http://www.biojava.org/docs/api/org/biojava/bio/structure/align/util/AtomCache.html) which manages your PDB (and btw. also SCOP, CATH) installations.

<pre>
        AtomCache cache = new AtomCache();
            
        cache.setUseMmCif(true);
            
        // if you struggled to set the PDB_DIR property correctly in the previous step, 
        // you could set it manually like this:
        cache.setPath("/tmp/");
            
        StructureIO.setAtomCache(cache);
            
        Structure structure = StructureIO.getStructure("4HHB");
                    
        // and let's count how many chains are in this structure.
        System.out.println(structure.getChains().size());
</pre>

As you can see, the AtomCache will again download the missing mmCIF file for 4HHB in the background. 

## Low level access

If you want to learn how to use the BioJava mmCIF parser to populate your own data structure, let's first take a look this lower-level code:

<pre>
        InputStream inStream =  new FileInputStream(fileName);
 
        MMcifParser parser = new SimpleMMcifParser();
 
        SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
 
        // The Consumer builds up the BioJava - structure object.
        // you could also hook in your own and build up you own data model.          
        parser.addMMcifConsumer(consumer);
 
        try {
            parser.parse(new BufferedReader(new InputStreamReader(inStream)));
        } catch (IOException e){
            e.printStackTrace();
        }
 
        // now get the protein structure.
        Structure cifStructure = consumer.getStructure();
</pre>

The parser operates similar to a XML parser by triggering "events". The [SimpleMMcifConsumer](http://www.biojava.org/docs/api/org/biojava/bio/structure/io/mmcif/SimpleMMcifConsumer.html) listens to new categories being read from the file and then builds up the BioJava data model.

To re-use the parser for your own datamodel, just implement the [MMcifConsumer](http://www.biojava.org/docs/api/org/biojava/bio/structure/io/mmcif/MMcifConsumer.html) interface and add it to the [SimpleMMcifParser](http://www.biojava.org/docs/api/org/biojava/bio/structure/io/mmcif/SimpleMMcifParser.html).
<pre>
        parser.addMMcifConsumer(myOwnConsumerImplementation);
</pre>

## I loaded a Structure object, what now?

BioJava provides a number of algorithms and visualisation tools that you can use to further analyse the structure, or look at it. Here a couple of suggestions for further reads:

+ [The BioJava Cookbook for protein structures](http://biojava.org/wiki/BioJava:CookBook#Protein_Structure)
+ How does BioJava [represent the content](structure-data-model.md) of a PDB/mmCIF file?
+ [How to calculate a protein structure alignment using BioJava](http://biojava.org/wiki/BioJava:CookBook:PDB:align)
+ [How to work with Groups (AminoAcid, Nucleotide, Hetatom)](http://biojava.org/wiki/BioJava:CookBook:PDB:groups)



<!-- References -->


[Westbrook 2000]: http://www.ncbi.nlm.nih.gov/pubmed/10842738 "Westbrook JD and Bourne PE. STAR/mmCIF: an ontology for macromolecular structure. Bioinformatics 2000 Feb; 16(2) 159-68. pmid:10842738." 

[Westbrook 2003]: http://www.ncbi.nlm.nih.gov/pubmed/12647386 "Westbrook JD and Fitzgerald PM. The PDB format, mmCIF, and other data formats. Methods Biochem Anal 2003; 44 161-79. pmid:12647386."

