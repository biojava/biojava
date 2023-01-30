# Welcome to <img src="logo-full.png" height="35"/>

![Build](https://github.com/biojava/biojava/actions/workflows/master.yml/badge.svg)
[![Version](http://img.shields.io/badge/version-6.1.0-blue.svg?style=flat)](https://github.com/biojava/biojava/releases/tag/biojava-6.1.0) [![License](http://img.shields.io/badge/license-LGPL_2.1-blue.svg?style=flat)](https://github.com/biojava/biojava/blob/master/LICENSE) [![Join the chat at https://gitter.im/biojava/biojava](https://badges.gitter.im/biojava/biojava.svg)](https://gitter.im/biojava/biojava?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


BioJava is an open-source project dedicated to providing a Java framework for **processing biological data**. It provides analytical and statistical routines, parsers for common file formats, reference implementations of popular algorithms, and allows the manipulation of sequences and 3D structures. The goal of the biojava project is to facilitate rapid application development for bioinformatics.

Please visit our [homepage](http://biojava.org/).

### Documentation

The [BioJava tutorial](https://github.com/biojava/biojava-tutorial) is a great place to get started. It is most complete for the biojava-structure module. 

The [BioJava Cookbook](http://biojava.org/wiki/BioJava:CookBook4.0/) contains an older and slightly outdated collection of simple examples that teach the basics for how to work with BioJava.

Full javadocs are available at the [BioJava website](http://biojava.org/docs/api).

### Maven Repository

BioJava release are available from Maven Central.

### Quick Installation

If you are using Maven you can add the BioJava repository by adding the following XML to your project pom.xml file:

```xml
    <dependencies>
      <dependency>
        <groupId>org.biojava</groupId>
        <artifactId>biojava-core</artifactId>
        <version>6.1.0</version>
      </dependency>
      <!-- other biojava modules as needed -->
    </dependencies>
```

### Mailing Lists

BioJava has one main mailing list. In order to avoid SPAM the list only accepts postings from list members. Anybody can become a list member, so please subscribe before you post. If you send without being subscribed your mail might get stuck in the moderation loop, which can cause several weeks of delay (no fun to read through all that spam).

#### biojava-l general discussion list

* [biojava-l@biojava.org](http://lists.open-bio.org/mailman/listinfo/biojava-l)

This list is intended for general discussion, advice, questions, offers of help, announcements, expressions of appreciation, bugs found in release code and requests for features.

#### biojava-dev developers list
 
A [dev mailing list](http://lists.open-bio.org/mailman/listinfo/biojava-dev) used to exist, but it has now been shut down. For dev discussions we now use github issues. Please search existing issues and if you don't find the answer to your question submit a new issue.

### Please cite

**BioJava 5: A community driven open-source bioinformatics library**<br/>
*Aleix Lafita, Spencer Bliven, Andreas Prlić, Dmytro Guzenko, Peter W. Rose, Anthony Bradley, Paolo Pavan, Douglas Myers-Turnbull, Yana Valasatava, Michael Heuer, Matt Larson, Stephen K. Burley, Jose M. Duarte* <br/>
[PLOS Computational Biology 15(2): e1006791](http://dx.plos.org/10.1371/journal.pcbi.1006791) <br/>
[![doi](http://img.shields.io/badge/doi-10.1371%2Fjournal.pcbi.1006791-blue.svg?style=flat)](https://doi.org/10.1371/journal.pcbi.1006791)
