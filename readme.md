# Welcome to BioJava v.3

[![Build Status](https://travis-ci.org/biojava/biojava.svg?branch=master)](https://travis-ci.org/biojava/biojava) [![Dependency Status](http://www.versioneye.com/user/projects/53f8b9f0e09da3dcbb000394/badge.svg?style=flat)](http://www.versioneye.com/user/projects/53f8b9f0e09da3dcbb000394) [![Version](http://img.shields.io/badge/version-3.1.0-blue.svg?style=flat)](http://biojava.org/wiki/BioJava:Download) [![License](http://img.shields.io/badge/license-LGPL_2.1-blue.svg?style=flat)](https://github.com/biojava/biojava/blob/master/LICENSE)


BioJava is an open-source project dedicated to providing a Java framework for **processing biological data**. It provides analytical and statistical routines, parsers for common file formats and allows the manipulation of sequences and 3D structures. The goal of the biojava project is to facilitate rapid application development for bioinformatics.

Please visit our [homepage](http://www.biojava.org/).

### Documentation

The [BioJava Cookbook](http://biojava.org/wiki/BioJava:CookBook) is a collection of simple examples that teach the basics for how to work with BioJava.
We are also working on a [tutorial](https://github.com/biojava/biojava3-tutorial) for biojava (currently only for protein structure modules). 

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
        <version>3.1.0</version>
      </dependency>
      <!-- other biojava jars as needed -->
    </dependencies>
```

Note: for the upcoming BioJava 4.0.0 release we are planning to migrate the hosting to Maven Central and this step should not be required.


### Mailing Lists

BioJava has two main mailing lists. In order to avoid SPAM both lists only accept postings from list members. Anybody can become a list member, so please subscribe before you post. If you send without being subscribed your mail might get stuck in the moderation loop, which can cause several weeks of delay (no fun to read through all that spam).

#### biojava-l general discussion list

* [biojava-l@biojava.org](http://lists.open-bio.org/mailman/listinfo/biojava-l)

This list is intended for general discussion, advice, questions, offers of help, announcements, expressions of appreciation, bugs found in release code and requests for features.

#### biojava-dev developers list
 
* [biojava-dev@biojava.org](http://lists.open-bio.org/mailman/listinfo/biojava-dev)

This list is intended for more technical discussions about API design, bugs in git development code, performance issues and things that might not be of interest to the more casual user.
