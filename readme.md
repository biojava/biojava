# Welcome to <img src="logo-full.png" height="35"/>

[![Build Status](https://travis-ci.org/biojava/biojava.svg?branch=master)](https://travis-ci.org/biojava/biojava) [![Dependency Status](http://www.versioneye.com/user/projects/5776f02568ee07004137f521/badge.svg?style=flat)](http://www.versioneye.com/user/projects/5776f02568ee07004137f521) [![Version](http://img.shields.io/badge/version-4.2.4-blue.svg?style=flat)](https://github.com/biojava/biojava/releases/tag/biojava-4.2.4) [![License](http://img.shields.io/badge/license-LGPL_2.1-blue.svg?style=flat)](https://github.com/biojava/biojava/blob/master/LICENSE) [![Join the chat at https://gitter.im/biojava/biojava](https://badges.gitter.im/biojava/biojava.svg)](https://gitter.im/biojava/biojava?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


BioJava is an open-source project dedicated to providing a Java framework for **processing biological data**. It provides analytical and statistical routines, parsers for common file formats, reference implementations of popular algorithms, and allows the manipulation of sequences and 3D structures. The goal of the biojava project is to facilitate rapid application development for bioinformatics.

Please visit our [homepage](http://www.biojava.org/).

### Documentation

The [BioJava Cookbook](http://biojava.org/wiki/BioJava:CookBook4.0/) is a collection of simple examples that teach the basics for how to work with BioJava.
We are also working on a [tutorial](https://github.com/biojava/biojava3-tutorial) for biojava (currently only for protein structure modules). 

### Maven Repository

BioJava release are available from Maven Central. Snapshot builds are distributed using [OSS Sonatype](https://oss.sonatype.org/content/repositories/snapshots/org/biojava)

### Quick Installation

If you are using Maven you can add the BioJava repository by adding the following XML to your project pom.xml file:

```xml
    <dependencies>
      <dependency>
        <groupId>org.biojava</groupId>
        <artifactId>biojava-core</artifactId>
        <version>4.2.4</version>
      </dependency>
      <!-- other biojava jars as needed -->
    </dependencies>
```

### Snapshot builds

To use the latest builds from BioJava, you can add the following config your project's pom.xml:

```xml
<repositories>
    <repository>
      <id>oss.sonatype.org-snapshot</id>
      <url>http://oss.sonatype.org/content/repositories/snapshots</url>
      <releases>
        <enabled>false</enabled>
      </releases>
      <snapshots>
        <enabled>true</enabled>
      </snapshots>
    </repository>
  </repositories>
  ```

### Mailing Lists

BioJava has two main mailing lists. In order to avoid SPAM both lists only accept postings from list members. Anybody can become a list member, so please subscribe before you post. If you send without being subscribed your mail might get stuck in the moderation loop, which can cause several weeks of delay (no fun to read through all that spam).

#### biojava-l general discussion list

* [biojava-l@biojava.org](http://lists.open-bio.org/mailman/listinfo/biojava-l)

This list is intended for general discussion, advice, questions, offers of help, announcements, expressions of appreciation, bugs found in release code and requests for features.

#### biojava-dev developers list
 
* [biojava-dev@biojava.org](http://lists.open-bio.org/mailman/listinfo/biojava-dev)

This list is intended for more technical discussions about API design, bugs in git development code, performance issues and things that might not be of interest to the more casual user.

### Please cite


**BioJava: an open-source framework for bioinformatics in 2012**<br/>
*Andreas Prlic; Andrew Yates; Spencer E. Bliven; Peter W. Rose; Julius Jacobsen; Peter V. Troshin; Mark Chapman; Jianjiong Gao; Chuan Hock Koh; Sylvain Foisy; Richard Holland; Gediminas Rimsa; Michael L. Heuer; H. Brandstatter-Muller; Philip E. Bourne; Scooter Willis* <br/>
[Bioinformatics (2012) 28 (20): 2693-2695.](http://bioinformatics.oxfordjournals.org/content/28/20/2693.abstract) <br/>
[![doi](http://img.shields.io/badge/doi-10.1093%2Fbioinformatics%2Fbts494-blue.svg?style=flat)](http://bioinformatics.oxfordjournals.org/content/28/20/2693.abstract) [![pubmed](http://img.shields.io/badge/pubmed-22877863-blue.svg?style=flat)](https://www.ncbi.nlm.nih.gov/pubmed/22877863)
