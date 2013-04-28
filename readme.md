# Welcome to BioJava v.3

BioJava is an open-source project dedicated to providing a Java framework for **processing biological data**. It provides analytical and statistical routines, parsers for common file formats and allows the manipulation of sequences and 3D structures. The goal of the biojava project is to facilitate rapid application development for bioinformatics.

BioJava is licensed under LGPL 2.1.

Please visit our [homepage](http://www.biojava.org/).

### Documentation

The [BioJava Cookbook](http://biojava.org/wiki/BioJava:CookBook) is a collection of simple examples that teach the basics for how to work with BioJava.

### Maven Repository

The [Maven Repository](http://biojava.org/download/maven/) contains the jars of all releases, as well as SNAPSHOT builds of the latest code base.

### Quick Installation

If you are using Maven you can add the BioJava repository by adding the following XML to your project pom.xml file:

``` 
    <repositories>
      <repository>
        <id>biojava-maven-repo</id>
        <name>BioJava repository</name>
        <url>http://www.biojava.org/download/maven/</url>			
      </repository>
    </repositories>

    <dependencies>
      <dependency>
        <groupId>org.biojava</groupId>
        <artifactId>biojava3-core</artifactId>
        <version>3.0.5</version>
      </dependency>
      <!-- other biojava jars as needed -->
    </dependencies>
```
### Build Status
[![Build Status](https://travis-ci.org/biojava/biojava.png)](https://travis-ci.org/biojava/biojava)

