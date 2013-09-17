# How to parse mmCIF files using BioJava

## mmCIF to be come the new principle deposition and dissemination file format

The wwPDB recently announced the transition from PDB to mmCIF/PDBx as  the principal deposition and dissemination file format (see 
[http://www.wwpdb.org/news/news_2013.html#22-May-2013](http://www.wwpdb.org/news/news_2013.html#22-May-2013) and 
[http://wwpdb.org/workshop/wgroup.html](http://wwpdb.org/workshop/wgroup.html)). [BioJava](http://www.biojava.org) has been containing support for mmCIF already for several years. This tutorial is meant to provide a quick introduction into how to parse mmCIF files using [BioJava](http://www.biojava.org)

## The basics

BioJava provides you with both a mmCIF parser and a data model that reads PDB and mmCIF files into a biological and chemically meaningful data model. If you don't want to use that data model, you can still use BioJava's file parsers, and more on that later, let's start first with the most basic way of loading a protein structure.

## Quick Installation

Just one paragraph of how to get access to BioJava. 
BioJava uses [Maven](http://maven.apache.org/) as a build and distribution system. If you are new to Maven, take a look at the [Getting Started with Maven](http://maven.apache.org/guides/getting-started/index.html)  guide.

We are providing a BioJava specific Maven repository at (http://biojava.org/download/maven/) .

You can add the BioJava repository by adding the following XML to your project pom.xml file:
<pre>
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
                        <groupId>org.biojava</groupId>
                        <artifactId>biojava3-core</artifactId>
                        <version>3.0.6</version>
                </dependency>
                <!-- other biojava jars as needed -->
        </dependencies>
    
</pre>

