# How to parse mmCIF files using BioJava

## What is mmCIF?

The Protein Data Bank (PDB) has been distributing its archival files as PDB files for a long time. The PDB file format is based on "punchcard"-style rules how to store data in a flat file. With the increasing complexity of macromolecules that have are being resolved experimentally, this file format can not be used any more to represent some or the more complex structures. As such, the wwPDB recently announced the transition from PDB to mmCIF/PDBx as  the principal deposition and dissemination file format (see 
[here](http://www.wwpdb.org/news/news_2013.html#22-May-2013) and 
[here](http://wwpdb.org/workshop/wgroup.html)). 

The mmCIF file format has been around for some time (see [Westbrook 2000][] and [Westbrook 2003][] ) [BioJava](http://www.biojava.org) has been supporting mmCIF already for several years. This tutorial is meant to provide a quick introduction into how to parse mmCIF files using [BioJava](http://www.biojava.org)

## The basics

BioJava provides you with both a mmCIF parser and a data model that reads PDB and mmCIF files into a biological and chemically meaningful data model. If you don't want to use that data model, you can still use BioJava's file parsers, and more on that later, let's start first with the most basic way of loading a protein structure.

## Quick Installation

Before we start, just one quick paragraph of how to get access to BioJava.

BioJava is open source and you can get the code from [Github](https://github.com/biojava/biojava), however it might be easier this way:

BioJava uses [Maven](http://maven.apache.org/) as a build and distribution system. If you are new to Maven, take a look at the [Getting Started with Maven](http://maven.apache.org/guides/getting-started/index.html)  guide.

We are providing a BioJava specific Maven repository at (http://biojava.org/download/maven/) .

You can add the BioJava repository by adding the following XML to your project pom.xml file:
<myxml>
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
    
</myxml>

If you run 'mvn package' on your project, the BioJava dependencies will be automatically downloaded and installed for you.

## First steps

The simplest way to load a PDB file is by using the [StructureIO](http://www.biojava.org/docs/api/org/biojava3/structure/StructureIO.html) class.

<pre>
    Structure structure = StructureIO.getStructure("4HHB");
</pre>

BioJava  automatically downloaded the PDB file for [4HHB](http://www.rcsb.org/pdb/explore.do?structureId=4HHB) and copied it into a temporary location. This demonstrates two things:

+ BioJava can automatically download and install files locally
+ BioJava by default writes those files into a temporary location (The system temp directory "java.io.tempdir"). 

You can configure where BioJava should read the files from by setting the PDB_DIR system property

<pre>
    -DPDB_DIR=/wherever/you/want/
</pre>


<!-- References -->


[Westbrook 2000]: http://www.ncbi.nlm.nih.gov/pubmed/10842738 "Westbrook JD and Bourne PE. STAR/mmCIF: an ontology for macromolecular structure. Bioinformatics 2000 Feb; 16(2) 159-68. pmid:10842738." 

[Westbrook 2003]: http://www.ncbi.nlm.nih.gov/pubmed/12647386 "Westbrook JD and Fitzgerald PM. The PDB format, mmCIF, and other data formats. Methods Biochem Anal 2003; 44 161-79. pmid:12647386."

