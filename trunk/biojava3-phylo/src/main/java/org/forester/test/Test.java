// $Id: Test.java,v 1.253 2009/12/30 04:33:45 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Set;

import org.forester.application.support_transfer;
import org.forester.development.SupportCount;
import org.forester.go.TestGo;
import org.forester.io.parsers.HmmscanPerDomainTableParser;
import org.forester.io.parsers.HmmscanPerDomainTableParser.INDIVIDUAL_SCORE_CUTOFF;
import org.forester.io.parsers.nexus.NexusBinaryStatesMatrixParser;
import org.forester.io.parsers.nexus.NexusCharactersParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.parsers.tol.TolParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.pccx.TestPccx;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyBranch;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.BinaryCharacters;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogenyinference.CharacterStateMatrix;
import org.forester.phylogenyinference.TestPhylogenyReconstruction;
import org.forester.phylogenyinference.CharacterStateMatrix.BinaryStates;
import org.forester.sdi.SDI;
import org.forester.sdi.SDIR;
import org.forester.sdi.SDIse;
import org.forester.sdi.TaxonomyAssigner;
import org.forester.sdi.TestGSDI;
import org.forester.surfacing.Protein;
import org.forester.surfacing.TestSurfacing;
import org.forester.tools.ConfidenceAssessor;
import org.forester.tools.TreeSplitMatrix;
import org.forester.util.AsciiHistogram;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.forester.util.GeneralTable;

@SuppressWarnings( "unused")
public final class Test {

    private final static double  ZERO_DIFF                 = 1.0E-9;
    private final static String  PATH_TO_TEST_DATA         = System.getProperty( "user.dir" )
                                                                   + ForesterUtil.getFileSeparator() + "test_data"
                                                                   + ForesterUtil.getFileSeparator();
    private final static String  PATH_TO_RESOURCES         = System.getProperty( "user.dir" )
                                                                   + ForesterUtil.getFileSeparator() + "resources"
                                                                   + ForesterUtil.getFileSeparator();
    private final static boolean USE_LOCAL_PHYLOXML_SCHEMA = true;
    private static final String  PHYLOXML_REMOTE_XSD       = ForesterConstants.PHYLO_XML_LOCATION + "/"
                                                                   + ForesterConstants.PHYLO_XML_VERSION + "/"
                                                                   + ForesterConstants.PHYLO_XML_XSD;
    private static final String  PHYLOXML_LOCAL_XSD        = PATH_TO_RESOURCES + "phyloxml_schema/"
                                                                   + ForesterConstants.PHYLO_XML_VERSION + "/"
                                                                   + ForesterConstants.PHYLO_XML_XSD;

    private final static Phylogeny createPhylogeny( final String nhx ) throws IOException {
        final Phylogeny p = ParserBasedPhylogenyFactory.getInstance().create( nhx, new NHXParser() )[ 0 ];
        return p;
    }

    private final static Event getEvent( final Phylogeny p, final String n1, final String n2 ) {
        final PhylogenyMethods pm = PhylogenyMethods.getInstance();
        return pm.getLCA( p.getNode( n1 ), p.getNode( n2 ) ).getNodeData().getEvent();
    }

    public static boolean isEqual( final double a, final double b ) {
        return ( ( Math.abs( a - b ) ) < Test.ZERO_DIFF );
    }

    public static void main( final String[] args ) {
        System.out.println( "[Java version: " + ForesterUtil.JAVA_VERSION + " " + ForesterUtil.JAVA_VENDOR + "]" );
        System.out.println( "[OS: " + ForesterUtil.OS_NAME + " " + ForesterUtil.OS_ARCH + " " + ForesterUtil.OS_VERSION
                + "]" );
        Locale.setDefault( Locale.US );
        System.out.println( "[Locale: " + Locale.getDefault() + "]" );
        int failed = 0;
        int succeeded = 0;
        System.out.print( "[Test if directory with files for testing exists/is readable: " );
        if ( Test.testDir( PATH_TO_TEST_DATA ) ) {
            System.out.println( "OK.]" );
        }
        else {
            System.out.println( "could not find/read from directory \"" + PATH_TO_TEST_DATA + "\".]" );
            System.out.println( "Testing aborted." );
            System.exit( -1 );
        }
        System.out.print( "[Test if resources directory exists/is readable: " );
        if ( testDir( PATH_TO_RESOURCES ) ) {
            System.out.println( "OK.]" );
        }
        else {
            System.out.println( "could not find/read from directory \"" + Test.PATH_TO_RESOURCES + "\".]" );
            System.out.println( "Testing aborted." );
            System.exit( -1 );
        }
        final long start_time = new Date().getTime();
        System.out.print( "Hmmscan output parser: " );
        if ( testHmmscanOutputParser() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Basic node methods: " );
        if ( Test.testBasicNodeMethods() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Basic node construction and parsing of NHX (node level): " );
        if ( Test.testNHXNodeParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "NH parsing: " );
        if ( Test.testNHParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Conversion to NHX (node level): " );
        if ( Test.testNHXconversion() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "NHX parsing: " );
        if ( Test.testNHXParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "NHX parsing with quotes: " );
        if ( Test.testNHXParsingQuotes() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Nexus characters parsing: " );
        if ( Test.testNexusCharactersParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Nexus tree parsing: " );
        if ( Test.testNexusTreeParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Nexus tree parsing (translating): " );
        if ( Test.testNexusTreeParsingTranslating() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Nexus matrix parsing: " );
        if ( Test.testNexusMatrixParsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Basic phyloXML parsing: " );
        if ( Test.testBasicPhyloXMLparsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Basic phyloXML parsing (validating against schema): " );
        if ( testBasicPhyloXMLparsingValidating() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Roundtrip phyloXML parsing (validating against schema): " );
        if ( Test.testBasicPhyloXMLparsingRoundtrip() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Tol XML parsing: " );
        if ( Test.testBasicTolXMLparsing() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Copying of node data: " );
        if ( Test.testCopyOfNodeData() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Basic tree methods: " );
        if ( Test.testBasicTreeMethods() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Postorder Iterator: " );
        if ( Test.testPostOrderIterator() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Preorder Iterator: " );
        if ( Test.testPreOrderIterator() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Levelorder Iterator: " );
        if ( Test.testLevelOrderIterator() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Re-id methods: " );
        if ( Test.testReIdMethods() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Methods on last external nodes: " );
        if ( Test.testLastExternalNodeMethods() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Methods on external nodes: " );
        if ( Test.testExternalNodeRelatedMethods() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Deletion of external nodes: " );
        if ( Test.testDeletionOfExternalNodes() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Subtree deletion: " );
        if ( Test.testSubtreeDeletion() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Phylogeny branch: " );
        if ( Test.testPhylogenyBranch() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Rerooting: " );
        if ( Test.testRerooting() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Mipoint rooting: " );
        if ( Test.testMidpointrooting() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Support count: " );
        if ( Test.testSupportCount() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Support transfer: " );
        if ( Test.testSupportTransfer() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Finding of LCA: " );
        if ( Test.testGetLCA() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Calculation of distance between nodes: " );
        if ( Test.testGetDistance() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "SDIse: " );
        if ( Test.testSDIse() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Taxonomy assigner: " );
        if ( Test.testTaxonomyAssigner() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "SDIunrooted: " );
        if ( Test.testSDIunrooted() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "GSDI: " );
        if ( TestGSDI.test() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Descriptive statistics: " );
        if ( Test.testDescriptiveStatistics() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Data objects and methods: " );
        if ( Test.testDataObjects() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Phylogeny reconstruction:" );
        System.out.println();
        if ( TestPhylogenyReconstruction.test( new File( PATH_TO_TEST_DATA ) ) ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Analysis of domain architectures: " );
        System.out.println();
        if ( TestSurfacing.test( new File( PATH_TO_TEST_DATA ) ) ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "GO: " );
        System.out.println();
        if ( TestGo.test( new File( PATH_TO_TEST_DATA ) ) ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Modeling tools: " );
        if ( TestPccx.test() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Split Matrix: " );
        if ( Test.testSplit() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Split Matrix: " );
        if ( Test.testConfidenceAssessor() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "Basic table: " );
        if ( Test.testBasicTable() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.print( "General table: " );
        if ( Test.testGeneralTable() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.println();
        final Runtime rt = java.lang.Runtime.getRuntime();
        final long free_memory = rt.freeMemory() / 1000000;
        final long total_memory = rt.totalMemory() / 1000000;
        System.out.println( "Running time    : " + ( new Date().getTime() - start_time ) + "ms " + "(free memory: "
                + free_memory + "MB, total memory: " + total_memory + "MB)" );
        System.out.println();
        System.out.println( "Successful tests: " + succeeded );
        System.out.println( "Failed     tests: " + failed );
        System.out.println();
        if ( failed < 1 ) {
            System.out.println( "OK." );
        }
        else {
            System.out.println( "Not OK." );
        }
        // System.out.println();
        // Development.setTime( true );
        //try {
        //  final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        //  final String clc = System.getProperty( "user.dir" ) + ForesterUtil.getFileSeparator()
        //          + "examples" + ForesterUtil.getFileSeparator() + "CLC.nhx";
        // final String multi = Test.PATH_TO_EXAMPLE_FILES +
        // "multifurcations_ex_1.nhx";
        // final String domains = Test.PATH_TO_EXAMPLE_FILES + "domains1.nhx";
        // final Phylogeny t1 = factory.create( new File( domains ), new
        // NHXParser() )[ 0 ];
        //  final Phylogeny t2 = factory.create( new File( clc ), new NHXParser() )[ 0 ];
        // }
        // catch ( final Exception e ) {
        //     e.printStackTrace();
        // }
        // t1.getRoot().preorderPrint();
        // final PhylogenyFactory factory = ParserBasedPhylogenyFactory
        // .getInstance();
        // try {
        //            
        // Helper.readNHtree( new File( PATH_TO_EXAMPLE_FILES
        // + "\\AtNBSpos.nhx" ) );
        // factory.create(
        // new File( PATH_TO_EXAMPLE_FILES + "\\AtNBSpos.nhx" ),
        // new NHXParser() );
        // Helper.readNHtree( new File( PATH_TO_EXAMPLE_FILES
        // + "\\AtNBSpos.nhx" ) );
        // factory.create(
        // new File( PATH_TO_EXAMPLE_FILES + "\\AtNBSpos.nhx" ),
        // new NHXParser() );
        //            
        //
        // Helper.readNHtree( new File( PATH_TO_EXAMPLE_FILES
        // + "\\big_tree.nhx" ) );
        // Helper.readNHtree( new File( PATH_TO_EXAMPLE_FILES
        // + "\\big_tree.nhx" ) );
        // factory.create(
        // new File( PATH_TO_EXAMPLE_FILES + "\\big_tree.nhx" ),
        // new NHXParser() );
        // factory.create(
        // new File( PATH_TO_EXAMPLE_FILES + "\\big_tree.nhx" ),
        // new NHXParser() );
        //
        // Helper.readNHtree( new File( PATH_TO_EXAMPLE_FILES
        // + "\\big_tree.nhx" ) );
        // Helper.readNHtree( new File( PATH_TO_EXAMPLE_FILES
        // + "\\big_tree.nhx" ) );
        //
        // factory.create(
        // new File( PATH_TO_EXAMPLE_FILES + "\\big_tree.nhx" ),
        // new NHXParser() );
        // factory.create(
        // new File( PATH_TO_EXAMPLE_FILES + "\\big_tree.nhx" ),
        // new NHXParser() );
        //
        // Helper.readNHtree( new File( PATH_TO_EXAMPLE_FILES
        // + "\\AtNBSpos.nhx" ) );
        // factory.create(
        // new File( PATH_TO_EXAMPLE_FILES + "\\AtNBSpos.nhx" ),
        // new NHXParser() );
        //
        // }
        // catch ( IOException e ) {
        // // TODO Auto-generated catch block
        // e.printStackTrace();
        // }
    }

    private static boolean testConfidenceAssessor() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0 = factory.create( "((((A,B)ab,C)abc,D)abcd,E)abcde", new NHXParser() )[ 0 ];
            final Phylogeny[] ev0 = factory
                    .create( "((((A,B),C),D),E);((((A,B),C),D),E);((((A,B),C),D),E);((((A,B),C),D),E);",
                             new NHXParser() );
            ConfidenceAssessor.evaluate( "bootstrap", ev0, t0, false, 1, 0, 2 );
            if ( !isEqual( t0.getNode( "ab" ).getBranchData().getConfidence( 0 ).getValue(), 3 ) ) {
                return false;
            }
            if ( !isEqual( t0.getNode( "abc" ).getBranchData().getConfidence( 0 ).getValue(), 3 ) ) {
                return false;
            }
            final Phylogeny t1 = factory.create( "((((A,B)ab,C)abc,D)abcd,E)abcde", new NHXParser() )[ 0 ];
            final Phylogeny[] ev1 = factory
                    .create( "((((A,B),C),D),E);((A,B),((E,D),C));(((A,B),C),(E,D));(A,(((E,D),C),B));(B,(A,((E,D),C)));(C,((E,D),(A,B)));(D,(E,((A,B),C)));",
                             new NHXParser() );
            ConfidenceAssessor.evaluate( "bootstrap", ev1, t1, false, 1 );
            if ( !isEqual( t1.getNode( "ab" ).getBranchData().getConfidence( 0 ).getValue(), 7 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "abc" ).getBranchData().getConfidence( 0 ).getValue(), 7 ) ) {
                return false;
            }
            final Phylogeny[] t2 = factory
                    .create( "((((a,b),c),d),e);(((a,b),c),(d,e));(((((a,b),c),d),e),f);((((a,b),c),(d,e)),f);(((a,b),c),d,e);((a,b,c),d,e);",
                             new NHXParser() );
            final Phylogeny[] ev2 = factory
                    .create( "((((a,b),c),d),e);((((a,b),c),d),e);((((a,b),e),d),c);((((a,b),e),d),c);(((a,b),(c,d)),e);((a,b),x);((a,b),(x,y));(a,b);(a,e);(a,b,c);",
                             new NHXParser() );
            for( final Phylogeny target : t2 ) {
                ConfidenceAssessor.evaluate( "bootstrap", ev2, target, false, 1 );
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testSplit() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p0 = factory.create( "(((A,B,C),D),(E,(F,G)))R", new NHXParser() )[ 0 ];
            //Archaeopteryx.createApplication( p0 );
            final Set<PhylogenyNode> ex = new HashSet<PhylogenyNode>();
            ex.add( new PhylogenyNode( "A" ) );
            ex.add( new PhylogenyNode( "B" ) );
            ex.add( new PhylogenyNode( "C" ) );
            ex.add( new PhylogenyNode( "D" ) );
            ex.add( new PhylogenyNode( "E" ) );
            ex.add( new PhylogenyNode( "F" ) );
            ex.add( new PhylogenyNode( "G" ) );
            final TreeSplitMatrix s0 = new TreeSplitMatrix( p0, false, ex );
            // System.out.println( s0.toString() );
            //
            Set<PhylogenyNode> query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "B" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            // query_nodes = new HashSet<PhylogenyNode>();
            // query_nodes.add( new PhylogenyNode( "A" ) );
            // if ( !s0.match( query_nodes ) ) {
            //      return false;
            //  }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "B" ) );
            query_nodes.add( new PhylogenyNode( "C" ) );
            query_nodes.add( new PhylogenyNode( "D" ) );
            query_nodes.add( new PhylogenyNode( "E" ) );
            query_nodes.add( new PhylogenyNode( "F" ) );
            query_nodes.add( new PhylogenyNode( "G" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "B" ) );
            query_nodes.add( new PhylogenyNode( "C" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "D" ) );
            query_nodes.add( new PhylogenyNode( "E" ) );
            query_nodes.add( new PhylogenyNode( "F" ) );
            query_nodes.add( new PhylogenyNode( "G" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "B" ) );
            query_nodes.add( new PhylogenyNode( "C" ) );
            query_nodes.add( new PhylogenyNode( "D" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "E" ) );
            query_nodes.add( new PhylogenyNode( "F" ) );
            query_nodes.add( new PhylogenyNode( "G" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "F" ) );
            query_nodes.add( new PhylogenyNode( "G" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "E" ) );
            query_nodes.add( new PhylogenyNode( "D" ) );
            query_nodes.add( new PhylogenyNode( "C" ) );
            query_nodes.add( new PhylogenyNode( "B" ) );
            query_nodes.add( new PhylogenyNode( "A" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "F" ) );
            query_nodes.add( new PhylogenyNode( "G" ) );
            query_nodes.add( new PhylogenyNode( "E" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "F" ) );
            query_nodes.add( new PhylogenyNode( "G" ) );
            query_nodes.add( new PhylogenyNode( "E" ) );
            query_nodes.add( new PhylogenyNode( "D" ) );
            if ( !s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "F" ) );
            query_nodes.add( new PhylogenyNode( "A" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "E" ) );
            query_nodes.add( new PhylogenyNode( "B" ) );
            query_nodes.add( new PhylogenyNode( "C" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "F" ) );
            query_nodes.add( new PhylogenyNode( "G" ) );
            query_nodes.add( new PhylogenyNode( "E" ) );
            query_nodes.add( new PhylogenyNode( "D" ) );
            query_nodes.add( new PhylogenyNode( "C" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "B" ) );
            query_nodes.add( new PhylogenyNode( "D" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "D" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "B" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "C" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "E" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "F" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "G" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "F" ) );
            query_nodes.add( new PhylogenyNode( "G" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "B" ) );
            query_nodes.add( new PhylogenyNode( "D" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "E" ) );
            query_nodes.add( new PhylogenyNode( "D" ) );
            query_nodes.add( new PhylogenyNode( "A" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            //
            query_nodes = new HashSet<PhylogenyNode>();
            query_nodes.add( new PhylogenyNode( "E" ) );
            query_nodes.add( new PhylogenyNode( "D" ) );
            query_nodes.add( new PhylogenyNode( "A" ) );
            query_nodes.add( new PhylogenyNode( "G" ) );
            if ( s0.match( query_nodes ) ) {
                return false;
            }
            /////////
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "A" ) );
            //            query_nodes.add( new PhylogenyNode( "B" ) );
            //            query_nodes.add( new PhylogenyNode( "C" ) );
            //            query_nodes.add( new PhylogenyNode( "D" ) );
            //            query_nodes.add( new PhylogenyNode( "E" ) );
            //            query_nodes.add( new PhylogenyNode( "F" ) );
            //            query_nodes.add( new PhylogenyNode( "G" ) );
            //            if ( !s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "A" ) );
            //            query_nodes.add( new PhylogenyNode( "B" ) );
            //            query_nodes.add( new PhylogenyNode( "C" ) );
            //            if ( !s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //            //
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "D" ) );
            //            query_nodes.add( new PhylogenyNode( "E" ) );
            //            query_nodes.add( new PhylogenyNode( "F" ) );
            //            query_nodes.add( new PhylogenyNode( "G" ) );
            //            if ( !s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //            //
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "A" ) );
            //            query_nodes.add( new PhylogenyNode( "B" ) );
            //            query_nodes.add( new PhylogenyNode( "C" ) );
            //            query_nodes.add( new PhylogenyNode( "D" ) );
            //            if ( !s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //            //
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "E" ) );
            //            query_nodes.add( new PhylogenyNode( "F" ) );
            //            query_nodes.add( new PhylogenyNode( "G" ) );
            //            if ( !s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //            //
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "F" ) );
            //            query_nodes.add( new PhylogenyNode( "G" ) );
            //            if ( !s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //            //
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "E" ) );
            //            query_nodes.add( new PhylogenyNode( "G" ) );
            //            if ( s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //            //
            //            query_nodes = new HashSet<PhylogenyNode>();
            //            query_nodes.add( new PhylogenyNode( "X" ) );
            //            query_nodes.add( new PhylogenyNode( "Y" ) );
            //            query_nodes.add( new PhylogenyNode( "A" ) );
            //            query_nodes.add( new PhylogenyNode( "B" ) );
            //            if ( s0.match( query_nodes ) ) {
            //                return false;
            //            }
            //
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testBasicNodeMethods() {
        try {
            if ( PhylogenyNode.getNodeCount() != 0 ) {
                return false;
            }
            final PhylogenyNode n1 = new PhylogenyNode();
            final PhylogenyNode n2 = new PhylogenyNode( "", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            final PhylogenyNode n3 = new PhylogenyNode( "n3", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            final PhylogenyNode n4 = new PhylogenyNode( "n4:0.01", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( n1.isHasAssignedEvent() ) {
                return false;
            }
            if ( PhylogenyNode.getNodeCount() != 4 ) {
                return false;
            }
            if ( n3.getIndicator() != 0 ) {
                return false;
            }
            if ( n3.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            if ( !n3.isExternal() ) {
                return false;
            }
            if ( !n3.isRoot() ) {
                return false;
            }
            if ( !n4.getNodeName().equals( "n4" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicPhyloXMLparsing() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhyloXmlParser xml_parser = new PhyloXmlParser();
            final Phylogeny[] phylogenies_0 = factory.create( Test.PATH_TO_TEST_DATA + "phyloxml_test_t1.xml",
                                                              xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_0.length != 4 ) {
                return false;
            }
            final Phylogeny t1 = phylogenies_0[ 0 ];
            final Phylogeny t2 = phylogenies_0[ 1 ];
            final Phylogeny t3 = phylogenies_0[ 2 ];
            final Phylogeny t4 = phylogenies_0[ 3 ];
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            if ( !t1.isRooted() ) {
                return false;
            }
            if ( t1.isRerootable() ) {
                return false;
            }
            if ( !t1.getType().equals( "gene_tree" ) ) {
                return false;
            }
            if ( t2.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "node a" ).getDistanceToParent(), 1.0 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "node b" ).getDistanceToParent(), 2.0 ) ) {
                return false;
            }
            if ( t3.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            if ( !t1.getName().equals( "t1" ) ) {
                return false;
            }
            if ( !t2.getName().equals( "t2" ) ) {
                return false;
            }
            if ( !t3.getName().equals( "t3" ) ) {
                return false;
            }
            if ( !t4.getName().equals( "t4" ) ) {
                return false;
            }
            if ( !t3.getIdentifier().getValue().equals( "1-1" ) ) {
                return false;
            }
            if ( !t3.getIdentifier().getProvider().equals( "treebank" ) ) {
                return false;
            }
            if ( !t3.getNode( "root node" ).getNodeData().getSequence().getType().equals( "protein" ) ) {
                return false;
            }
            if ( !t3.getNode( "root node" ).getNodeData().getSequence().getName()
                    .equals( "Apoptosis facilitator Bcl-2-like 14 protein" ) ) {
                return false;
            }
            if ( !t3.getNode( "root node" ).getNodeData().getSequence().getSymbol().equals( "BCL2L14" ) ) {
                return false;
            }
            if ( !t3.getNode( "root node" ).getNodeData().getSequence().getAccession().getValue().equals( "Q9BZR8" ) ) {
                return false;
            }
            if ( !t3.getNode( "root node" ).getNodeData().getSequence().getAccession().getSource().equals( "UniProtKB" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) ).getDesc()
                    .equals( "apoptosis" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) ).getRef()
                    .equals( "GO:0006915" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getSource().equals( "UniProtKB" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getEvidence().equals( "experimental" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) ).getType()
                    .equals( "function" ) ) {
                return false;
            }
            if ( ( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getConfidence().getValue() != 1 ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getConfidence().getType().equals( "ml" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) ).getDesc()
                    .equals( "apoptosis" ) ) {
                return false;
            }
            if ( ( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getProperties().getProperty( "AFFY:expression" ).getAppliesTo() != AppliesTo.ANNOTATION ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getProperties().getProperty( "AFFY:expression" ).getDataType().equals( "xsd:double" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getProperties().getProperty( "AFFY:expression" ).getRef().equals( "AFFY:expression" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getProperties().getProperty( "AFFY:expression" ).getUnit().equals( "AFFY:x" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getProperties().getProperty( "AFFY:expression" ).getValue().equals( "0.2" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getProperties().getProperty( "MED:disease" ).getValue().equals( "lymphoma" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 1 ) ).getRef()
                    .equals( "GO:0005829" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) ).getDesc()
                    .equals( "intracellular organelle" ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getUri().getType().equals( "source" ) ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getUri().getDescription()
                    .equals( "UniProt link" ) ) ) {
                return false;
            }
            if ( !( t3.getNode( "root node" ).getNodeData().getSequence().getLocation().equals( "12p13-p12" ) ) ) {
                return false;
            }
            //if ( !( t3.getNode( "root node" ).getNodeData().getDistribution().getDesc().equals( "irgendwo" ) ) ) {
            //     return false;
            //}
            //            if ( !( t3.getNode( "root node" ).getNodeData().getReference().getDoi().equals( "10.1074/jbc.M005889200" ) ) ) {
            //                return false;
            //            }
            //            if ( !t3.getNode( "root node" ).getNodeData().getTaxonomy().getType().equals( "host" ) ) {
            //                return false;
            //            }
            //            if ( !t3.getNode( "root node" ).getNodeData().getTaxonomy().getTaxonomyCode().equals( "ECDYS" ) ) {
            //                return false;
            //            }
            //            if ( !t3.getNode( "root node" ).getNodeData().getTaxonomy().getScientificName().equals( "ecdysozoa" ) ) {
            //                return false;
            //            }
            //            if ( !t3.getNode( "root node" ).getNodeData().getTaxonomy().getCommonName().equals( "molting animals" ) ) {
            //                return false;
            //            }
            //            if ( !t3.getNode( "root node" ).getNodeData().getTaxonomy().getIdentifier().getValue().equals( "1" ) ) {
            //                return false;
            //            }
            //            if ( !t3.getNode( "root node" ).getNodeData().getTaxonomy().getIdentifier().getType().equals( "ncbi" ) ) {
            //                return false;
            //            }
            //            if ( t3.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getTotalLength() != 124 ) {
            //                return false;
            //            }
            //            if ( !t3.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 ).getName()
            //                    .equals( "B" ) ) {
            //                return false;
            //            }
            //            if ( t3.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 ).getFrom() != 21 ) {
            //                return false;
            //            }
            //            if ( t3.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 ).getTo() != 44 ) {
            //                return false;
            //            }
            //            if ( t3.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 ).getLength() != 24 ) {
            //                return false;
            //            }
            //            if ( t3.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 )
            //                    .getConfidence() != 2144 ) {
            //                return false;
            //            }
            //            if ( !t3.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 ).getId()
            //                    .equals( "pfam" ) ) {
            //                return false;
            //            }
            //            if ( t3.getNode( "node bb" ).getNodeData().getBinaryCharacters().getGainedCharacters().size() != 3 ) {
            //                return false;
            //            }
            //            if ( t3.getNode( "node bb" ).getNodeData().getBinaryCharacters().getPresentCharacters().size() != 2 ) {
            //                return false;
            //            }
            //            if ( t3.getNode( "node bb" ).getNodeData().getBinaryCharacters().getLostCharacters().size() != 1 ) {
            //                return false;
            //            }
            //            if ( !t3.getNode( "node bb" ).getNodeData().getBinaryCharacters().getType().equals( "domains" ) ) {
            //                return false;
            //            }
            //            if ( ( ( BinaryCharacters ) t3.getNode( "node bb" ).getNodeData().getBinaryCharacters().copy() )
            //                    .getLostCount() != BinaryCharacters.COUNT_DEFAULT ) {
            //                ;
            //                return false;
            //            }
            //            if ( t3.getNode( "node b" ).getNodeData().getBinaryCharacters().getGainedCount() != 1 ) {
            //                return false;
            //            }
            //            if ( t3.getNode( "node b" ).getNodeData().getBinaryCharacters().getGainedCharacters().size() != 1 ) {
            //                return false;
            //            }
            //            if ( t3.getNode( "node b" ).getNodeData().getBinaryCharacters().getLostCount() != 3 ) {
            //                return false;
            //            }
            //            if ( t3.getNode( "node b" ).getNodeData().getBinaryCharacters().getLostCharacters().size() != 3 ) {
            //                return false;
            //            }
            //            if ( t3.getNode( "node b" ).getNodeData().getBinaryCharacters().getPresentCount() != 2 ) {
            //                return false;
            //            }
            //            if ( t3.getNode( "node b" ).getNodeData().getBinaryCharacters().getPresentCharacters().size() != 2 ) {
            //                return false;
            //            }
            //            if ( !t3.getNode( "node b" ).getNodeData().getBinaryCharacters().getType().equals( "characters" ) ) {
            //                return false;
            //            }
            //            final Phylogeny[] phylogenies_1 = factory.create( Test.PATH_TO_TEST_DATA + "phyloxml_test_t4.xml",
            //                                                              xml_parser );
            //            if ( xml_parser.getErrorCount() > 0 ) {
            //                System.out.println( xml_parser.getErrorMessages().toString() );
            //                return false;
            //            }
            //            if ( phylogenies_1.length != 2 ) {
            //                return false;
            //            }
            //            final Phylogeny a = phylogenies_1[ 0 ];
            //            if ( !a.getName().equals( "tree 4" ) ) {
            //                return false;
            //            }
            //            if ( a.getNumberOfExternalNodes() != 3 ) {
            //                return false;
            //            }
            //            if ( !a.getNode( "node b1" ).getNodeData().getSequence().getName().equals( "b1 gene" ) ) {
            //                return false;
            //            }
            //            if ( !a.getNode( "node b1" ).getNodeData().getTaxonomy().getCommonName().equals( "b1 species" ) ) {
            //                return false;
            //            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicPhyloXMLparsingRoundtrip() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhyloXmlParser xml_parser = new PhyloXmlParser();
            if ( USE_LOCAL_PHYLOXML_SCHEMA ) {
                xml_parser.setValidateAgainstSchema( PHYLOXML_LOCAL_XSD );
            }
            else {
                xml_parser.setValidateAgainstSchema( PHYLOXML_REMOTE_XSD );
            }
            final Phylogeny[] phylogenies_0 = factory.create( Test.PATH_TO_TEST_DATA + "phyloxml_test_t1.xml",
                                                              xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_0.length != 4 ) {
                return false;
            }
            final StringBuffer t1_sb = new StringBuffer( phylogenies_0[ 0 ].toPhyloXML( 0 ) );
            final Phylogeny[] phylogenies_t1 = factory.create( t1_sb, xml_parser );
            if ( phylogenies_t1.length != 1 ) {
                return false;
            }
            final Phylogeny t1_rt = phylogenies_t1[ 0 ];
            if ( !t1_rt.getDistanceUnit().equals( "cc" ) ) {
                return false;
            }
            if ( !t1_rt.isRooted() ) {
                return false;
            }
            if ( t1_rt.isRerootable() ) {
                return false;
            }
            if ( !t1_rt.getType().equals( "gene_tree" ) ) {
                return false;
            }
            final StringBuffer t3_sb_0 = new StringBuffer( phylogenies_0[ 2 ].toPhyloXML( 0 ) );
            final Phylogeny[] phylogenies_1_0 = factory.create( t3_sb_0, xml_parser );
            final StringBuffer t3_sb = new StringBuffer( phylogenies_1_0[ 0 ].toPhyloXML( 0 ) );
            final Phylogeny[] phylogenies_1 = factory.create( t3_sb, xml_parser );
            if ( phylogenies_1.length != 1 ) {
                return false;
            }
            final Phylogeny t3_rt = phylogenies_1[ 0 ];
            if ( !t3_rt.getName().equals( "t3" ) ) {
                return false;
            }
            if ( t3_rt.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            if ( !t3_rt.getIdentifier().getValue().equals( "1-1" ) ) {
                return false;
            }
            if ( !t3_rt.getIdentifier().getProvider().equals( "treebank" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getSequence().getType().equals( "protein" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getSequence().getName()
                    .equals( "Apoptosis facilitator Bcl-2-like 14 protein" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getSequence().getSymbol().equals( "BCL2L14" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getSequence().getAccession().getValue().equals( "Q9BZR8" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getSequence().getAccession().getSource()
                    .equals( "UniProtKB" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getDesc().equals( "apoptosis" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getRef().equals( "GO:0006915" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getSource().equals( "UniProtKB" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getEvidence().equals( "experimental" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getType().equals( "function" ) ) {
                return false;
            }
            if ( ( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getConfidence().getValue() != 1 ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getConfidence().getType().equals( "ml" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getDesc().equals( "apoptosis" ) ) {
                return false;
            }
            if ( ( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getProperties().getProperty( "AFFY:expression" ).getAppliesTo() != AppliesTo.ANNOTATION ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getProperties().getProperty( "AFFY:expression" ).getDataType().equals( "xsd:double" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getProperties().getProperty( "AFFY:expression" ).getRef().equals( "AFFY:expression" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getProperties().getProperty( "AFFY:expression" ).getUnit().equals( "AFFY:x" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getProperties().getProperty( "AFFY:expression" ).getValue().equals( "0.2" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 0 ) )
                    .getProperties().getProperty( "MED:disease" ).getValue().equals( "lymphoma" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 1 ) )
                    .getRef().equals( "GO:0005829" ) ) {
                return false;
            }
            if ( !( ( Annotation ) t3_rt.getNode( "root node" ).getNodeData().getSequence().getAnnotation( 2 ) )
                    .getDesc().equals( "intracellular organelle" ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getUri().getType().equals( "source" ) ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getUri().getDescription()
                    .equals( "UniProt link" ) ) ) {
                return false;
            }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getSequence().getLocation().equals( "12p13-p12" ) ) ) {
                return false;
            }
            //  if ( !( t3_rt.getNode( "root node" ).getNodeData().getDistribution().getDesc().equals( "irgendwo" ) ) ) {
            //      return false;
            //  }
            if ( !( t3_rt.getNode( "root node" ).getNodeData().getReference().getDoi()
                    .equals( "10.1074/jbc.M005889200" ) ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getTaxonomy().getTaxonomyCode().equals( "ECDYS" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getTaxonomy().getScientificName().equals( "ecdysozoa" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getTaxonomy().getCommonName().equals( "molting animals" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getTaxonomy().getIdentifier().getValue().equals( "1" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "root node" ).getNodeData().getTaxonomy().getIdentifier().getProvider()
                    .equals( "ncbi" ) ) {
                return false;
            }
            if ( t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getTotalLength() != 124 ) {
                return false;
            }
            if ( !t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 )
                    .getName().equals( "B" ) ) {
                return false;
            }
            if ( t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 )
                    .getFrom() != 21 ) {
                return false;
            }
            if ( t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 ).getTo() != 44 ) {
                return false;
            }
            if ( t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 )
                    .getLength() != 24 ) {
                return false;
            }
            if ( t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 )
                    .getConfidence() != 2144 ) {
                return false;
            }
            if ( !t3_rt.getNode( "node bc" ).getNodeData().getSequence().getDomainArchitecture().getDomain( 0 ).getId()
                    .equals( "pfam" ) ) {
                return false;
            }
            if ( t3_rt.getNode( "node bb" ).getNodeData().getBinaryCharacters().getGainedCharacters().size() != 3 ) {
                return false;
            }
            if ( t3_rt.getNode( "node bb" ).getNodeData().getBinaryCharacters().getPresentCharacters().size() != 2 ) {
                return false;
            }
            if ( t3_rt.getNode( "node bb" ).getNodeData().getBinaryCharacters().getLostCharacters().size() != 1 ) {
                return false;
            }
            if ( !t3_rt.getNode( "node bb" ).getNodeData().getBinaryCharacters().getType().equals( "domains" ) ) {
                return false;
            }
            final Taxonomy taxbb = t3_rt.getNode( "node bb" ).getNodeData().getTaxonomy();
            if ( !taxbb.getAuthority().equals( "Stephenson, 1935" ) ) {
                return false;
            }
            if ( !taxbb.getCommonName().equals( "starlet sea anemone" ) ) {
                return false;
            }
            if ( !taxbb.getIdentifier().getProvider().equals( "EOL" ) ) {
                return false;
            }
            if ( !taxbb.getIdentifier().getValue().equals( "704294" ) ) {
                return false;
            }
            if ( !taxbb.getTaxonomyCode().equals( "NEMVE" ) ) {
                return false;
            }
            if ( !taxbb.getScientificName().equals( "Nematostella vectensis" ) ) {
                return false;
            }
            if ( taxbb.getSynonyms().size() != 2 ) {
                return false;
            }
            if ( !taxbb.getSynonyms().contains( "Nematostella vectensis Stephenson1935" ) ) {
                return false;
            }
            if ( !taxbb.getSynonyms().contains( "See Anemone" ) ) {
                return false;
            }
            if ( !taxbb.getUri().getDescription().equals( "EOL" ) ) {
                return false;
            }
            if ( !taxbb.getUri().getType().equals( "linkout" ) ) {
                return false;
            }
            if ( !taxbb.getUri().getValue().toString().equals( "http://www.eol.org/pages/704294" ) ) {
                return false;
            }
            if ( ( ( BinaryCharacters ) t3_rt.getNode( "node bb" ).getNodeData().getBinaryCharacters().copy() )
                    .getLostCount() != BinaryCharacters.COUNT_DEFAULT ) {
                ;
                return false;
            }
            if ( t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getGainedCount() != 1 ) {
                return false;
            }
            if ( t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getGainedCharacters().size() != 1 ) {
                return false;
            }
            if ( t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getLostCount() != 3 ) {
                return false;
            }
            if ( t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getLostCharacters().size() != 3 ) {
                return false;
            }
            if ( t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getPresentCount() != 2 ) {
                return false;
            }
            if ( t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getPresentCharacters().size() != 2 ) {
                return false;
            }
            if ( !t3_rt.getNode( "node b" ).getNodeData().getBinaryCharacters().getType().equals( "characters" ) ) {
                return false;
            }
            //
            if ( !t3_rt.getNode( "node ba" ).getNodeData().getDate().getDesc().equals( "Silurian" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "node ba" ).getNodeData().getDate().getValue().toPlainString()
                    .equalsIgnoreCase( "435" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "node ba" ).getNodeData().getDate().getMin().toPlainString().equalsIgnoreCase( "416" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "node ba" ).getNodeData().getDate().getMax().toPlainString()
                    .equalsIgnoreCase( "443.7" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "node ba" ).getNodeData().getDate().getUnit().equals( "mya" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "node bb" ).getNodeData().getDate().getDesc().equals( "Triassic" ) ) {
                return false;
            }
            if ( !t3_rt.getNode( "node bc" ).getNodeData().getDate().getValue().toPlainString()
                    .equalsIgnoreCase( "433" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicPhyloXMLparsingValidating() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            PhyloXmlParser xml_parser = null;
            try {
                xml_parser = PhyloXmlParser.createPhyloXmlParserXsdValidating();
            }
            catch ( final Exception e ) {
                // Do nothing -- means were not running from jar.
            }
            if ( xml_parser == null ) {
                xml_parser = new PhyloXmlParser();
                if ( USE_LOCAL_PHYLOXML_SCHEMA ) {
                    xml_parser.setValidateAgainstSchema( PHYLOXML_LOCAL_XSD );
                }
                else {
                    xml_parser.setValidateAgainstSchema( PHYLOXML_REMOTE_XSD );
                }
            }
            final Phylogeny[] phylogenies_0 = factory.create( Test.PATH_TO_TEST_DATA + "phyloxml_test_t1.xml",
                                                              xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_0.length != 4 ) {
                return false;
            }
            final Phylogeny t1 = phylogenies_0[ 0 ];
            final Phylogeny t2 = phylogenies_0[ 1 ];
            final Phylogeny t3 = phylogenies_0[ 2 ];
            final Phylogeny t4 = phylogenies_0[ 3 ];
            if ( !t1.getName().equals( "t1" ) ) {
                return false;
            }
            if ( !t2.getName().equals( "t2" ) ) {
                return false;
            }
            if ( !t3.getName().equals( "t3" ) ) {
                return false;
            }
            if ( !t4.getName().equals( "t4" ) ) {
                return false;
            }
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            if ( t2.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            if ( t3.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            final String x2 = Test.PATH_TO_TEST_DATA + "phyloxml_test_t1.xml";
            final Phylogeny[] phylogenies_1 = factory.create( x2, xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( "errors:" );
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_1.length != 4 ) {
                return false;
            }
            final Phylogeny[] phylogenies_2 = factory.create( Test.PATH_TO_TEST_DATA + "phyloxml_test_t3.xml",
                                                              xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( "errors:" );
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_2.length != 1 ) {
                return false;
            }
            if ( phylogenies_2[ 0 ].getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            final Phylogeny[] phylogenies_3 = factory.create( Test.PATH_TO_TEST_DATA + "phyloxml_test_t4.xml",
                                                              xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_3.length != 2 ) {
                return false;
            }
            final Phylogeny a = phylogenies_3[ 0 ];
            if ( !a.getName().equals( "tree 4" ) ) {
                return false;
            }
            if ( a.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !a.getNode( "node b1" ).getNodeData().getSequence().getName().equals( "b1 gene" ) ) {
                return false;
            }
            if ( !a.getNode( "node b1" ).getNodeData().getTaxonomy().getCommonName().equals( "b1 species" ) ) {
                return false;
            }
            final Phylogeny[] phylogenies_4 = factory.create( Test.PATH_TO_TEST_DATA + "special_characters.xml",
                                                              xml_parser );
            if ( xml_parser.getErrorCount() > 0 ) {
                System.out.println( xml_parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_4.length != 1 ) {
                return false;
            }
            final Phylogeny s = phylogenies_4[ 0 ];
            if ( s.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            s.getNode( "first" );
            s.getNode( "<>" );
            s.getNode( "\"<a'b&c'd\">\"" );
            s.getNode( "'''\"" );
            s.getNode( "\"\"\"" );
            s.getNode( "dick & doof" );
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicTable() {
        try {
            final BasicTable<String> t0 = new BasicTable<String>();
            if ( t0.getNumberOfColumns() != 0 ) {
                return false;
            }
            if ( t0.getNumberOfRows() != 0 ) {
                return false;
            }
            t0.setValue( 3, 2, "23" );
            t0.setValue( 10, 1, "error" );
            t0.setValue( 10, 1, "110" );
            t0.setValue( 9, 1, "19" );
            t0.setValue( 1, 10, "101" );
            t0.setValue( 10, 10, "1010" );
            t0.setValue( 100, 10, "10100" );
            t0.setValue( 0, 0, "00" );
            if ( !t0.getValue( 3, 2 ).equals( "23" ) ) {
                return false;
            }
            if ( !t0.getValue( 10, 1 ).equals( "110" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 1, 10 ).equals( "101" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 10, 10 ).equals( "1010" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 100, 10 ).equals( "10100" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 9, 1 ).equals( "19" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 0, 0 ).equals( "00" ) ) {
                return false;
            }
            if ( t0.getNumberOfColumns() != 101 ) {
                return false;
            }
            if ( t0.getNumberOfRows() != 11 ) {
                return false;
            }
            if ( t0.getValueAsString( 49, 4 ) != null ) {
                return false;
            }
            final String l = ForesterUtil.getLineSeparator();
            final StringBuffer source = new StringBuffer();
            source.append( "" + l );
            source.append( "# 1 1 1 1 1 1 1 1" + l );
            source.append( " 00 01 02 03" + l );
            source.append( "   10 11 12 13  " + l );
            source.append( "20 21 22 23 " + l );
            source.append( "    30  31    32 33" + l );
            source.append( "40 41 42 43" + l );
            source.append( "  # 1 1 1 1 1 " + l );
            source.append( "50 51 52 53 54" + l );
            final BasicTable<String> t1 = BasicTableParser.parse( source.toString(), " " );
            if ( t1.getNumberOfColumns() != 5 ) {
                return false;
            }
            if ( t1.getNumberOfRows() != 6 ) {
                return false;
            }
            if ( !t1.getValueAsString( 0, 0 ).equals( "00" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( 1, 0 ).equals( "01" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( 3, 0 ).equals( "03" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( 4, 5 ).equals( "54" ) ) {
                return false;
            }
            final StringBuffer source1 = new StringBuffer();
            source1.append( "" + l );
            source1.append( "# 1; 1; 1; 1 ;1 ;1; 1 ;1;" + l );
            source1.append( " 00; 01 ;02;03" + l );
            source1.append( "   10; 11; 12; 13  " + l );
            source1.append( "20; 21; 22; 23 " + l );
            source1.append( "    30;  31;    32; 33" + l );
            source1.append( "40;41;42;43" + l );
            source1.append( "  # 1 1 1 1 1 " + l );
            source1.append( ";;;50  ;  ;52; 53;;54   " + l );
            final BasicTable<String> t2 = BasicTableParser.parse( source1.toString(), ";" );
            if ( t2.getNumberOfColumns() != 5 ) {
                return false;
            }
            if ( t2.getNumberOfRows() != 6 ) {
                return false;
            }
            if ( !t2.getValueAsString( 0, 0 ).equals( "00" ) ) {
                return false;
            }
            if ( !t2.getValueAsString( 1, 0 ).equals( "01" ) ) {
                return false;
            }
            if ( !t2.getValueAsString( 3, 0 ).equals( "03" ) ) {
                return false;
            }
            if ( !t2.getValueAsString( 3, 3 ).equals( "33" ) ) {
                return false;
            }
            if ( !t2.getValueAsString( 3, 5 ).equals( "53" ) ) {
                return false;
            }
            if ( !t2.getValueAsString( 1, 5 ).equals( "" ) ) {
                return false;
            }
            final StringBuffer source2 = new StringBuffer();
            source2.append( "" + l );
            source2.append( "comment: 1; 1; 1; 1 ;1 ;1; 1 ;1;" + l );
            source2.append( " 00; 01 ;02;03" + l );
            source2.append( "   10; 11; 12; 13  " + l );
            source2.append( "20; 21; 22; 23 " + l );
            source2.append( "                     " + l );
            source2.append( "    30;  31;    32; 33" + l );
            source2.append( "40;41;42;43" + l );
            source2.append( "  comment: 1 1 1 1 1 " + l );
            source2.append( ";;;50  ;   52; 53;;54   " + l );
            final List<BasicTable<String>> tl = BasicTableParser.parse( source2.toString(),
                                                                        ";",
                                                                        false,
                                                                        "comment:",
                                                                        false );
            if ( tl.size() != 2 ) {
                return false;
            }
            final BasicTable<String> t3 = tl.get( 0 );
            final BasicTable<String> t4 = tl.get( 1 );
            if ( t3.getNumberOfColumns() != 4 ) {
                return false;
            }
            if ( t3.getNumberOfRows() != 3 ) {
                return false;
            }
            if ( t4.getNumberOfColumns() != 4 ) {
                return false;
            }
            if ( t4.getNumberOfRows() != 3 ) {
                return false;
            }
            if ( !t3.getValueAsString( 0, 0 ).equals( "00" ) ) {
                return false;
            }
            if ( !t4.getValueAsString( 0, 0 ).equals( "30" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicTolXMLparsing() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final TolParser parser = new TolParser();
            final Phylogeny[] phylogenies_0 = factory.create( Test.PATH_TO_TEST_DATA + "tol_2484.tol", parser );
            if ( parser.getErrorCount() > 0 ) {
                System.out.println( parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_0.length != 1 ) {
                return false;
            }
            final Phylogeny t1 = phylogenies_0[ 0 ];
            if ( t1.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            if ( !t1.isRooted() ) {
                return false;
            }
            if ( !t1.getRoot().getNodeData().getTaxonomy().getCommonName().equals( "Mesozoa" ) ) {
                return false;
            }
            if ( !t1.getRoot().getNodeData().getTaxonomy().getIdentifier().getValue().equals( "2484" ) ) {
                return false;
            }
            if ( !t1.getRoot().getChildNode( 0 ).getNodeData().getTaxonomy().getCommonName().equals( "Rhombozoa" ) ) {
                return false;
            }
            if ( t1.getRoot().getChildNode( 0 ).getNumberOfDescendants() != 3 ) {
                return false;
            }
            final Phylogeny[] phylogenies_1 = factory.create( Test.PATH_TO_TEST_DATA + "tol_2.tol", parser );
            if ( parser.getErrorCount() > 0 ) {
                System.out.println( parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_1.length != 1 ) {
                return false;
            }
            final Phylogeny t2 = phylogenies_1[ 0 ];
            if ( t2.getNumberOfExternalNodes() != 664 ) {
                return false;
            }
            if ( !t2.isRooted() ) {
                return false;
            }
            if ( !t2.getRoot().getNodeData().getTaxonomy().getCommonName().equals( "Eubacteria" ) ) {
                return false;
            }
            if ( !t2.getRoot().getNodeData().getTaxonomy().getIdentifier().getValue().equals( "2" ) ) {
                return false;
            }
            if ( t2.getRoot().getNumberOfDescendants() != 24 ) {
                return false;
            }
            if ( t2.getRoot().getNumberOfDescendants() != 24 ) {
                return false;
            }
            if ( !t2.getRoot().getChildNode( 0 ).getNodeData().getTaxonomy().getCommonName().equals( "Aquificae" ) ) {
                return false;
            }
            if ( !t2.getRoot().getChildNode( 0 ).getChildNode( 0 ).getNodeData().getTaxonomy().getScientificName()
                    .equals( "Aquifex" ) ) {
                return false;
            }
            final Phylogeny[] phylogenies_2 = factory.create( Test.PATH_TO_TEST_DATA + "tol_5.tol", parser );
            if ( parser.getErrorCount() > 0 ) {
                System.out.println( parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_2.length != 1 ) {
                return false;
            }
            final Phylogeny t3 = phylogenies_2[ 0 ];
            if ( t3.getNumberOfExternalNodes() != 184 ) {
                return false;
            }
            if ( !t3.getRoot().getNodeData().getTaxonomy().getCommonName().equals( "Viruses" ) ) {
                return false;
            }
            if ( !t3.getRoot().getNodeData().getTaxonomy().getIdentifier().getValue().equals( "5" ) ) {
                return false;
            }
            if ( t3.getRoot().getNumberOfDescendants() != 6 ) {
                return false;
            }
            final Phylogeny[] phylogenies_3 = factory.create( Test.PATH_TO_TEST_DATA + "tol_4567.tol", parser );
            if ( parser.getErrorCount() > 0 ) {
                System.out.println( parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_3.length != 1 ) {
                return false;
            }
            final Phylogeny t4 = phylogenies_3[ 0 ];
            if ( t4.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            if ( !t4.getRoot().getNodeData().getTaxonomy().getCommonName().equals( "Marpissa decorata" ) ) {
                return false;
            }
            if ( !t4.getRoot().getNodeData().getTaxonomy().getIdentifier().getValue().equals( "4567" ) ) {
                return false;
            }
            if ( t4.getRoot().getNumberOfDescendants() != 0 ) {
                return false;
            }
            final Phylogeny[] phylogenies_4 = factory.create( Test.PATH_TO_TEST_DATA + "tol_16299.tol", parser );
            if ( parser.getErrorCount() > 0 ) {
                System.out.println( parser.getErrorMessages().toString() );
                return false;
            }
            if ( phylogenies_4.length != 1 ) {
                return false;
            }
            final Phylogeny t5 = phylogenies_4[ 0 ];
            if ( t5.getNumberOfExternalNodes() != 13 ) {
                return false;
            }
            if ( !t5.getRoot().getNodeData().getTaxonomy().getCommonName().equals( "Hominidae" ) ) {
                return false;
            }
            if ( !t5.getRoot().getNodeData().getTaxonomy().getIdentifier().getValue().equals( "16299" ) ) {
                return false;
            }
            if ( t5.getRoot().getNumberOfDescendants() != 2 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicTreeMethods() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t1 = factory.create();
            if ( !t1.isEmpty() ) {
                return false;
            }
            final Phylogeny t2 = factory.create( "((A:1,B:2)AB:1,(C:3,D:5)CD:3)ABCD:0.5", new NHXParser() )[ 0 ];
            if ( t2.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            if ( t2.getHeight() != 8.5 ) {
                return false;
            }
            if ( !t2.isCompletelyBinary() ) {
                return false;
            }
            if ( t2.isEmpty() ) {
                return false;
            }
            final Phylogeny t3 = factory.create( "((A:1,B:2,C:10)ABC:1,(D:3,E:5)DE:3)", new NHXParser() )[ 0 ];
            if ( t3.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            if ( t3.getHeight() != 11 ) {
                return false;
            }
            if ( t3.isCompletelyBinary() ) {
                return false;
            }
            final PhylogenyNode n = t3.getNode( "ABC" );
            PhylogenyNodeIterator it;
            for( it = n.iterateChildNodesForward(); it.hasNext(); ) {
                it.next();
            }
            for( it.reset(); it.hasNext(); ) {
                it.next();
            }
            final PhylogenyNodeIterator it2 = n.iterateChildNodesForward();
            if ( !it2.next().getNodeName().equals( "A" ) ) {
                return false;
            }
            if ( !it2.next().getNodeName().equals( "B" ) ) {
                return false;
            }
            if ( !it2.next().getNodeName().equals( "C" ) ) {
                return false;
            }
            if ( it2.hasNext() ) {
                return false;
            }
            final Phylogeny t4 = factory.create( "((A:1,B:2,C:10)ABC:1,(D:3,E:5)DE:3,(F,G,H,I))", new NHXParser() )[ 0 ];
            if ( t4.getNumberOfExternalNodes() != 9 ) {
                return false;
            }
            if ( t4.getHeight() != 11 ) {
                return false;
            }
            if ( t4.isCompletelyBinary() ) {
                return false;
            }
            final StringBuffer sb5 = new StringBuffer( "(((A11:2)A1:2,(A21:1,A22:2,A23)A2:11,A3:2)A:2,B:10,C:3,D:8)" );
            final Phylogeny t5 = factory.create( sb5, new NHXParser() )[ 0 ];
            if ( t5.getNumberOfExternalNodes() != 8 ) {
                return false;
            }
            if ( t5.getHeight() != 15 ) {
                return false;
            }
            final StringBuffer sb6 = new StringBuffer( "(X,Y,Z,(((A111)A11:2)A1:2,(X,Y,Z,A21:1,A22:2,A23)A2:11,A3:2)A:2,B:10,C:3,D:8)" );
            final Phylogeny t6 = factory.create( sb6, new NHXParser() )[ 0 ];
            if ( t6.getHeight() != 15 ) {
                return false;
            }
            final StringBuffer sb7 = new StringBuffer( "(((A11:2)A1:2,(A21:1,A22:2,A23)A2:11,A3:2)A:2,B:10,C:15,D:8)" );
            final Phylogeny t7 = factory.create( sb7, new NHXParser() )[ 0 ];
            if ( t7.getHeight() != 15 ) {
                return false;
            }
            final StringBuffer sb8 = new StringBuffer( "(((A11:11)A1:2,(A21:2,A22:2,A23,A24,AA:)A2:11,A3:2)A:2,B:15,C:15,D:15)" );
            final Phylogeny t8 = factory.create( sb8, new NHXParser() )[ 0 ];
            if ( t8.getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( t8.getHeight() != 15 ) {
                return false;
            }
            final char[] a9 = new char[] {};
            final Phylogeny t9 = factory.create( a9, new NHXParser() )[ 0 ];
            if ( t9.getHeight() != 0 ) {
                return false;
            }
            final char[] a10 = new char[] { 'a', ':', '6' };
            final Phylogeny t10 = factory.create( a10, new NHXParser() )[ 0 ];
            if ( t10.getHeight() != 6 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testCopyOfNodeData() {
        try {
            final PhylogenyNode n1 = new PhylogenyNode( "n5:0.1[&&NHX:S=Ecoli:E=1.1.1.1:D=Y:Co=Y:B=56:T=1:O=22:SO=33:SN=44:W=2:C=10.20.30:XN=S=tag1=value1=unit1]" );
            final PhylogenyNode n2 = n1.copyNodeData();
            if ( !n1.toNewHampshireX().equals( n2.toNewHampshireX() ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDataObjects() {
        try {
            final Confidence s0 = new Confidence();
            final Confidence s1 = new Confidence();
            if ( !s0.isEqual( s1 ) ) {
                return false;
            }
            final Confidence s2 = new Confidence( 0.23, "bootstrap" );
            final Confidence s3 = new Confidence( 0.23, "bootstrap" );
            if ( s2.isEqual( s1 ) ) {
                return false;
            }
            if ( !s2.isEqual( s3 ) ) {
                return false;
            }
            final Confidence s4 = ( Confidence ) s3.copy();
            if ( !s4.isEqual( s3 ) ) {
                return false;
            }
            s3.asSimpleText();
            s3.asText();
            // Taxonomy
            // ----------
            final Taxonomy t1 = new Taxonomy();
            final Taxonomy t2 = new Taxonomy();
            final Taxonomy t3 = new Taxonomy();
            final Taxonomy t4 = new Taxonomy();
            final Taxonomy t5 = new Taxonomy();
            t1.setIdentifier( new Identifier( "ecoli" ) );
            t1.setTaxonomyCode( "ECOLI" );
            t1.setScientificName( "E. coli" );
            t1.setCommonName( "coli" );
            final Taxonomy t0 = ( Taxonomy ) t1.copy();
            if ( !t1.isEqual( t0 ) ) {
                return false;
            }
            t2.setIdentifier( new Identifier( "ecoli" ) );
            t2.setTaxonomyCode( "other" );
            t2.setScientificName( "what" );
            t2.setCommonName( "something" );
            if ( !t1.isEqual( t2 ) ) {
                return false;
            }
            t2.setIdentifier( new Identifier( "nemve" ) );
            if ( t1.isEqual( t2 ) ) {
                return false;
            }
            t1.setIdentifier( null );
            t3.setTaxonomyCode( "ECOLI" );
            t3.setScientificName( "what" );
            t3.setCommonName( "something" );
            if ( !t1.isEqual( t3 ) ) {
                return false;
            }
            t1.setIdentifier( null );
            t1.setTaxonomyCode( "" );
            t4.setScientificName( "E. ColI" );
            t4.setCommonName( "something" );
            if ( !t1.isEqual( t4 ) ) {
                return false;
            }
            t4.setScientificName( "B. subtilis" );
            t4.setCommonName( "something" );
            if ( t1.isEqual( t4 ) ) {
                return false;
            }
            t1.setIdentifier( null );
            t1.setTaxonomyCode( "" );
            t1.setScientificName( "" );
            t5.setCommonName( "COLI" );
            if ( !t1.isEqual( t5 ) ) {
                return false;
            }
            t5.setCommonName( "vibrio" );
            if ( t1.isEqual( t5 ) ) {
                return false;
            }
            // Identifier
            // ----------
            final Identifier id0 = new Identifier( "123", "pfam" );
            final Identifier id1 = ( Identifier ) id0.copy();
            if ( !id1.isEqual( id1 ) ) {
                return false;
            }
            if ( !id1.isEqual( id0 ) ) {
                return false;
            }
            if ( !id0.isEqual( id1 ) ) {
                return false;
            }
            id1.asSimpleText();
            id1.asText();
            // ProteinDomain
            // ---------------
            final ProteinDomain pd0 = new ProteinDomain( "abc", 100, 200 );
            final ProteinDomain pd1 = ( ProteinDomain ) pd0.copy();
            if ( !pd1.isEqual( pd1 ) ) {
                return false;
            }
            if ( !pd1.isEqual( pd0 ) ) {
                return false;
            }
            pd1.asSimpleText();
            pd1.asText();
            final ProteinDomain pd2 = new ProteinDomain( pd0.getName(), pd0.getFrom(), pd0.getTo(), "id" );
            final ProteinDomain pd3 = ( ProteinDomain ) pd2.copy();
            if ( !pd3.isEqual( pd3 ) ) {
                return false;
            }
            if ( !pd2.isEqual( pd3 ) ) {
                return false;
            }
            if ( !pd0.isEqual( pd3 ) ) {
                return false;
            }
            pd3.asSimpleText();
            pd3.asText();
            // DomainArchitecture
            // ------------------
            final ProteinDomain d0 = new ProteinDomain( "domain0", 10, 20 );
            final ProteinDomain d1 = new ProteinDomain( "domain1", 30, 40 );
            final ProteinDomain d2 = new ProteinDomain( "domain2", 50, 60 );
            final ProteinDomain d3 = new ProteinDomain( "domain3", 70, 80 );
            final ProteinDomain d4 = new ProteinDomain( "domain4", 90, 100 );
            final ArrayList<PhylogenyData> domains0 = new ArrayList<PhylogenyData>();
            domains0.add( d2 );
            domains0.add( d0 );
            domains0.add( d3 );
            domains0.add( d1 );
            final DomainArchitecture ds0 = new DomainArchitecture( domains0, 110 );
            if ( ds0.getNumberOfDomains() != 4 ) {
                return false;
            }
            final DomainArchitecture ds1 = ( DomainArchitecture ) ds0.copy();
            if ( !ds0.isEqual( ds0 ) ) {
                return false;
            }
            if ( !ds0.isEqual( ds1 ) ) {
                return false;
            }
            if ( ds1.getNumberOfDomains() != 4 ) {
                return false;
            }
            final ArrayList<PhylogenyData> domains1 = new ArrayList<PhylogenyData>();
            domains1.add( d1 );
            domains1.add( d2 );
            domains1.add( d4 );
            domains1.add( d0 );
            final DomainArchitecture ds2 = new DomainArchitecture( domains1, 200 );
            if ( ds0.isEqual( ds2 ) ) {
                return false;
            }
            ds1.asSimpleText();
            ds1.asText();
            ds1.toNHX();
            final DomainArchitecture ds3 = new DomainArchitecture( "120>30>40>0.9>b>50>60>0.4>c>10>20>0.1>a" );
            if ( !ds3.toNHX().toString().equals( ":DS=120>10>20>0.1>a>30>40>0.9>b>50>60>0.4>c" ) ) {
                System.out.println( ds3.toNHX() );
                return false;
            }
            if ( ds3.getNumberOfDomains() != 3 ) {
                return false;
            }
            // Event
            // -----
            final Event e1 = new Event( Event.EventType.fusion );
            if ( e1.isDuplication() ) {
                return false;
            }
            if ( !e1.isFusion() ) {
                return false;
            }
            if ( !e1.asText().toString().equals( "fusion" ) ) {
                return false;
            }
            if ( !e1.asSimpleText().toString().equals( "fusion" ) ) {
                return false;
            }
            final Event e11 = new Event( Event.EventType.fusion );
            if ( !e11.isEqual( e1 ) ) {
                return false;
            }
            if ( !e11.toNHX().toString().equals( "" ) ) {
                return false;
            }
            final Event e2 = new Event( Event.EventType.speciation_or_duplication );
            if ( e2.isDuplication() ) {
                return false;
            }
            if ( !e2.isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !e2.asText().toString().equals( "speciation_or_duplication" ) ) {
                return false;
            }
            if ( !e2.asSimpleText().toString().equals( "?" ) ) {
                return false;
            }
            if ( !e2.toNHX().toString().equals( ":D=?" ) ) {
                return false;
            }
            if ( e11.isEqual( e2 ) ) {
                return false;
            }
            final Event e2c = ( Event ) e2.copy();
            if ( !e2c.isEqual( e2 ) ) {
                return false;
            }
            Event e3 = new Event( 1, 2, 3 );
            if ( e3.isDuplication() ) {
                return false;
            }
            if ( e3.isSpeciation() ) {
                return false;
            }
            if ( e3.isGeneLoss() ) {
                return false;
            }
            if ( !e3.asText().toString().equals( "duplications [1] speciations [2] gene-losses [3]" ) ) {
                return false;
            }
            final Event e3c = ( Event ) e3.copy();
            final Event e3cc = ( Event ) e3c.copy();
            if ( !e3c.asSimpleText().toString().equals( "D2S3L" ) ) {
                return false;
            }
            e3 = null;
            if ( !e3c.isEqual( e3cc ) ) {
                return false;
            }
            Event e4 = new Event( 1, 2, 3 );
            if ( !e4.asText().toString().equals( "duplications [1] speciations [2] gene-losses [3]" ) ) {
                return false;
            }
            if ( !e4.asSimpleText().toString().equals( "D2S3L" ) ) {
                return false;
            }
            final Event e4c = ( Event ) e4.copy();
            e4 = null;
            final Event e4cc = ( Event ) e4c.copy();
            if ( !e4cc.asText().toString().equals( "duplications [1] speciations [2] gene-losses [3]" ) ) {
                return false;
            }
            if ( !e4c.isEqual( e4cc ) ) {
                return false;
            }
            final Event e5 = new Event();
            if ( !e5.isUnassigned() ) {
                return false;
            }
            if ( !e5.asText().toString().equals( "unassigned" ) ) {
                return false;
            }
            if ( !e5.asSimpleText().toString().equals( "" ) ) {
                return false;
            }
            final Event e6 = new Event( 1, 0, 0 );
            if ( !e6.asText().toString().equals( "duplication" ) ) {
                return false;
            }
            if ( !e6.asSimpleText().toString().equals( "D" ) ) {
                return false;
            }
            final Event e7 = new Event( 0, 1, 0 );
            if ( !e7.asText().toString().equals( "speciation" ) ) {
                return false;
            }
            if ( !e7.asSimpleText().toString().equals( "S" ) ) {
                return false;
            }
            final Event e8 = new Event( 0, 0, 1 );
            if ( !e8.asText().toString().equals( "gene-loss" ) ) {
                return false;
            }
            if ( !e8.asSimpleText().toString().equals( "L" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDeletionOfExternalNodes() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0 = factory.create( "A", new NHXParser() )[ 0 ];
            final PhylogenyWriter w = new PhylogenyWriter();
            if ( t0.isEmpty() ) {
                return false;
            }
            if ( t0.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            t0.deleteSubtree( t0.getNode( "A" ), false );
            if ( t0.getNumberOfExternalNodes() != 0 ) {
                return false;
            }
            if ( !t0.isEmpty() ) {
                return false;
            }
            final Phylogeny t1 = factory.create( "(A,B)r", new NHXParser() )[ 0 ];
            if ( t1.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "A" ), false );
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            if ( !t1.getNode( "B" ).getNodeName().equals( "B" ) ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "B" ), false );
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "r" ), false );
            if ( !t1.isEmpty() ) {
                return false;
            }
            final Phylogeny t2 = factory.create( "((A,B),C)", new NHXParser() )[ 0 ];
            if ( t2.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            t2.deleteSubtree( t2.getNode( "B" ), false );
            if ( t2.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            t2.toNewHampshireX();
            PhylogenyNode n = t2.getNode( "A" );
            if ( !n.getNextExternalNode().getNodeName().equals( "C" ) ) {
                return false;
            }
            t2.deleteSubtree( t2.getNode( "A" ), false );
            if ( t2.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            t2.deleteSubtree( t2.getNode( "C" ), true );
            if ( t2.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            final Phylogeny t3 = factory.create( "((A,B),(C,D))", new NHXParser() )[ 0 ];
            if ( t3.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            t3.deleteSubtree( t3.getNode( "B" ), true );
            if ( t3.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            n = t3.getNode( "A" );
            if ( !n.getNextExternalNode().getNodeName().equals( "C" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getNextExternalNode().getNodeName().equals( "D" ) ) {
                return false;
            }
            t3.deleteSubtree( t3.getNode( "A" ), true );
            if ( t3.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            n = t3.getNode( "C" );
            if ( !n.getNextExternalNode().getNodeName().equals( "D" ) ) {
                return false;
            }
            t3.deleteSubtree( t3.getNode( "C" ), true );
            if ( t3.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            t3.deleteSubtree( t3.getNode( "D" ), true );
            if ( t3.getNumberOfExternalNodes() != 0 ) {
                return false;
            }
            final Phylogeny t4 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            if ( t4.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            t4.deleteSubtree( t4.getNode( "B2" ), true );
            if ( t4.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            String s = w.toNewHampshire( t4, false, true ).toString();
            if ( !s.equals( "((A,(B11,B12)),(C,D));" ) ) {
                return false;
            }
            t4.deleteSubtree( t4.getNode( "B11" ), true );
            if ( t4.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            t4.deleteSubtree( t4.getNode( "C" ), true );
            if ( t4.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            n = t4.getNode( "A" );
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "B12" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "D" ) ) {
                return false;
            }
            s = w.toNewHampshire( t4, false, true ).toString();
            if ( !s.equals( "((A,B12),D);" ) ) {
                return false;
            }
            final Phylogeny t5 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            t5.deleteSubtree( t5.getNode( "A" ), true );
            if ( t5.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t5, false, true ).toString();
            if ( !s.equals( "(((B11,B12),B2),(C,D));" ) ) {
                return false;
            }
            final Phylogeny t6 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            t6.deleteSubtree( t6.getNode( "B11" ), true );
            if ( t6.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t6, false, false ).toString();
            if ( !s.equals( "((A,(B12,B2)),(C,D));" ) ) {
                return false;
            }
            final Phylogeny t7 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            t7.deleteSubtree( t7.getNode( "B12" ), true );
            if ( t7.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t7, false, true ).toString();
            if ( !s.equals( "((A,(B11,B2)),(C,D));" ) ) {
                return false;
            }
            final Phylogeny t8 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            t8.deleteSubtree( t8.getNode( "B2" ), true );
            if ( t8.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t8, false, false ).toString();
            if ( !s.equals( "((A,(B11,B12)),(C,D));" ) ) {
                return false;
            }
            final Phylogeny t9 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            t9.deleteSubtree( t9.getNode( "C" ), true );
            if ( t9.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t9, false, true ).toString();
            if ( !s.equals( "((A,((B11,B12),B2)),D);" ) ) {
                return false;
            }
            final Phylogeny t10 = factory.create( "((A,((B11,B12),B2)),(C,D))", new NHXParser() )[ 0 ];
            t10.deleteSubtree( t10.getNode( "D" ), true );
            if ( t10.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t10, false, true ).toString();
            if ( !s.equals( "((A,((B11,B12),B2)),C);" ) ) {
                return false;
            }
            final Phylogeny t11 = factory.create( "(A,B,C)", new NHXParser() )[ 0 ];
            t11.deleteSubtree( t11.getNode( "A" ), true );
            if ( t11.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            s = w.toNewHampshire( t11, false, true ).toString();
            if ( !s.equals( "(B,C);" ) ) {
                return false;
            }
            t11.deleteSubtree( t11.getNode( "C" ), true );
            if ( t11.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            s = w.toNewHampshire( t11, false, false ).toString();
            if ( !s.equals( "B;" ) ) {
                return false;
            }
            final Phylogeny t12 = factory.create( "((A1,A2,A3),(B1,B2,B3),(C1,C2,C3))", new NHXParser() )[ 0 ];
            t12.deleteSubtree( t12.getNode( "B2" ), true );
            if ( t12.getNumberOfExternalNodes() != 8 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "((A1,A2,A3),(B1,B3),(C1,C2,C3));" ) ) {
                return false;
            }
            t12.deleteSubtree( t12.getNode( "B3" ), true );
            if ( t12.getNumberOfExternalNodes() != 7 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "((A1,A2,A3),B1,(C1,C2,C3));" ) ) {
                return false;
            }
            t12.deleteSubtree( t12.getNode( "C3" ), true );
            if ( t12.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "((A1,A2,A3),B1,(C1,C2));" ) ) {
                return false;
            }
            t12.deleteSubtree( t12.getNode( "A1" ), true );
            if ( t12.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "((A2,A3),B1,(C1,C2));" ) ) {
                return false;
            }
            t12.deleteSubtree( t12.getNode( "B1" ), true );
            if ( t12.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "((A2,A3),(C1,C2));" ) ) {
                return false;
            }
            t12.deleteSubtree( t12.getNode( "A3" ), true );
            if ( t12.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "(A2,(C1,C2));" ) ) {
                return false;
            }
            t12.deleteSubtree( t12.getNode( "A2" ), true );
            if ( t12.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            s = w.toNewHampshire( t12, false, true ).toString();
            if ( !s.equals( "(C1,C2);" ) ) {
                return false;
            }
            final Phylogeny t13 = factory.create( "(A,B,C,(D:1.0,E:2.0):3.0)", new NHXParser() )[ 0 ];
            t13.deleteSubtree( t13.getNode( "D" ), true );
            if ( t13.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            s = w.toNewHampshire( t13, false, true ).toString();
            if ( !s.equals( "(A,B,C,E:5.0);" ) ) {
                return false;
            }
            final Phylogeny t14 = factory.create( "((A,B,C,(D:0.1,E:0.4):1.0),F)", new NHXParser() )[ 0 ];
            t14.deleteSubtree( t14.getNode( "E" ), true );
            if ( t14.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            s = w.toNewHampshire( t14, false, true ).toString();
            if ( !s.equals( "((A,B,C,D:1.1),F);" ) ) {
                return false;
            }
            final Phylogeny t15 = factory.create( "((A1,A2,A3,A4),(B1,B2,B3,B4),(C1,C2,C3,C4))", new NHXParser() )[ 0 ];
            t15.deleteSubtree( t15.getNode( "B2" ), true );
            if ( t15.getNumberOfExternalNodes() != 11 ) {
                return false;
            }
            t15.deleteSubtree( t15.getNode( "B1" ), true );
            if ( t15.getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            t15.deleteSubtree( t15.getNode( "B3" ), true );
            if ( t15.getNumberOfExternalNodes() != 9 ) {
                return false;
            }
            t15.deleteSubtree( t15.getNode( "B4" ), true );
            if ( t15.getNumberOfExternalNodes() != 8 ) {
                return false;
            }
            t15.deleteSubtree( t15.getNode( "A1" ), true );
            if ( t15.getNumberOfExternalNodes() != 7 ) {
                return false;
            }
            t15.deleteSubtree( t15.getNode( "C4" ), true );
            if ( t15.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDescriptiveStatistics() {
        try {
            final DescriptiveStatistics dss1 = new BasicDescriptiveStatistics();
            dss1.addValue( 82 );
            dss1.addValue( 78 );
            dss1.addValue( 70 );
            dss1.addValue( 58 );
            dss1.addValue( 42 );
            if ( dss1.getN() != 5 ) {
                return false;
            }
            if ( !Test.isEqual( dss1.getMin(), 42 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.getMax(), 82 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.arithmeticMean(), 66 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.sampleStandardDeviation(), 16.24807680927192 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.median(), 70 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.midrange(), 62 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.sampleVariance(), 264 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.pearsonianSkewness(), -0.7385489458759964 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.coefficientOfVariation(), 0.24618298195866547 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.sampleStandardUnit( 66 - 16.24807680927192 ), -1.0 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.getValue( 1 ), 78 ) ) {
                return false;
            }
            dss1.addValue( 123 );
            if ( !Test.isEqual( dss1.arithmeticMean(), 75.5 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.getMax(), 123 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss1.standardErrorOfMean(), 11.200446419674531 ) ) {
                return false;
            }
            final DescriptiveStatistics dss2 = new BasicDescriptiveStatistics();
            dss2.addValue( -1.85 );
            dss2.addValue( 57.5 );
            dss2.addValue( 92.78 );
            dss2.addValue( 57.78 );
            if ( !Test.isEqual( dss2.median(), 57.64 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss2.sampleStandardDeviation(), 39.266984753946495 ) ) {
                return false;
            }
            final double[] a = dss2.getDataAsDoubleArray();
            if ( !Test.isEqual( a[ 3 ], 57.78 ) ) {
                return false;
            }
            dss2.addValue( -100 );
            if ( !Test.isEqual( dss2.sampleStandardDeviation(), 75.829111296388 ) ) {
                return false;
            }
            if ( !Test.isEqual( dss2.sampleVariance(), 5750.05412 ) ) {
                return false;
            }
            final double[] ds = new double[ 14 ];
            ds[ 0 ] = 34;
            ds[ 1 ] = 23;
            ds[ 2 ] = 1;
            ds[ 3 ] = 32;
            ds[ 4 ] = 11;
            ds[ 5 ] = 2;
            ds[ 6 ] = 12;
            ds[ 7 ] = 33;
            ds[ 8 ] = 13;
            ds[ 9 ] = 22;
            ds[ 10 ] = 21;
            ds[ 11 ] = 35;
            ds[ 12 ] = 24;
            ds[ 13 ] = 31;
            final int[] bins = BasicDescriptiveStatistics.performBinning( ds, 0, 40, 4 );
            if ( bins.length != 4 ) {
                return false;
            }
            if ( bins[ 0 ] != 2 ) {
                return false;
            }
            if ( bins[ 1 ] != 3 ) {
                return false;
            }
            if ( bins[ 2 ] != 4 ) {
                return false;
            }
            if ( bins[ 3 ] != 5 ) {
                return false;
            }
            final double[] ds1 = new double[ 9 ];
            ds1[ 0 ] = 10.0;
            ds1[ 1 ] = 19.0;
            ds1[ 2 ] = 9.999;
            ds1[ 3 ] = 0.0;
            ds1[ 4 ] = 39.9;
            ds1[ 5 ] = 39.999;
            ds1[ 6 ] = 30.0;
            ds1[ 7 ] = 19.999;
            ds1[ 8 ] = 30.1;
            final int[] bins1 = BasicDescriptiveStatistics.performBinning( ds1, 0, 40, 4 );
            if ( bins1.length != 4 ) {
                return false;
            }
            if ( bins1[ 0 ] != 2 ) {
                return false;
            }
            if ( bins1[ 1 ] != 3 ) {
                return false;
            }
            if ( bins1[ 2 ] != 0 ) {
                return false;
            }
            if ( bins1[ 3 ] != 4 ) {
                return false;
            }
            final int[] bins1_1 = BasicDescriptiveStatistics.performBinning( ds1, 0, 40, 3 );
            if ( bins1_1.length != 3 ) {
                return false;
            }
            if ( bins1_1[ 0 ] != 3 ) {
                return false;
            }
            if ( bins1_1[ 1 ] != 2 ) {
                return false;
            }
            if ( bins1_1[ 2 ] != 4 ) {
                return false;
            }
            final int[] bins1_2 = BasicDescriptiveStatistics.performBinning( ds1, 1, 39, 3 );
            if ( bins1_2.length != 3 ) {
                return false;
            }
            if ( bins1_2[ 0 ] != 2 ) {
                return false;
            }
            if ( bins1_2[ 1 ] != 2 ) {
                return false;
            }
            if ( bins1_2[ 2 ] != 2 ) {
                return false;
            }
            final DescriptiveStatistics dss3 = new BasicDescriptiveStatistics();
            dss3.addValue( 1 );
            dss3.addValue( 1 );
            dss3.addValue( 1 );
            dss3.addValue( 2 );
            dss3.addValue( 3 );
            dss3.addValue( 4 );
            dss3.addValue( 5 );
            dss3.addValue( 5 );
            dss3.addValue( 5 );
            dss3.addValue( 6 );
            dss3.addValue( 7 );
            dss3.addValue( 8 );
            dss3.addValue( 9 );
            dss3.addValue( 10 );
            dss3.addValue( 10 );
            dss3.addValue( 10 );
            final AsciiHistogram histo = new AsciiHistogram( dss3 );
            histo.toStringBuffer( 10, '=', 40, 5 );
            histo.toStringBuffer( 3, 8, 10, '=', 40, 5 );
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDir( final String file ) {
        try {
            final File f = new File( file );
            if ( !f.exists() ) {
                return false;
            }
            if ( !f.isDirectory() ) {
                return false;
            }
            if ( !f.canRead() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            return false;
        }
        return true;
    }

    private static boolean testExternalNodeRelatedMethods() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t1 = factory.create( "((A,B),(C,D))", new NHXParser() )[ 0 ];
            PhylogenyNode n = t1.getNode( "A" );
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "B" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "C" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "D" ) ) {
                return false;
            }
            n = t1.getNode( "B" );
            while ( !n.isLastExternalNode() ) {
                n = n.getNextExternalNode();
            }
            final Phylogeny t2 = factory.create( "(((A,B),C),D)", new NHXParser() )[ 0 ];
            n = t2.getNode( "A" );
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "B" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "C" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "D" ) ) {
                return false;
            }
            n = t2.getNode( "B" );
            while ( !n.isLastExternalNode() ) {
                n = n.getNextExternalNode();
            }
            final Phylogeny t3 = factory.create( "(((A,B),(C,D)),((E,F),(G,H)))", new NHXParser() )[ 0 ];
            n = t3.getNode( "A" );
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "B" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "C" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "D" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "E" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "F" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "G" ) ) {
                return false;
            }
            n = n.getNextExternalNode();
            if ( !n.getNodeName().equals( "H" ) ) {
                return false;
            }
            n = t3.getNode( "B" );
            while ( !n.isLastExternalNode() ) {
                n = n.getNextExternalNode();
            }
            final Phylogeny t4 = factory.create( "((A,B),(C,D))", new NHXParser() )[ 0 ];
            for( final PhylogenyNodeIterator iter = t4.iteratorExternalForward(); iter.hasNext(); ) {
                final PhylogenyNode node = iter.next();
            }
            final Phylogeny t5 = factory.create( "(((A,B),(C,D)),((E,F),(G,H)))", new NHXParser() )[ 0 ];
            for( final PhylogenyNodeIterator iter = t5.iteratorExternalForward(); iter.hasNext(); ) {
                final PhylogenyNode node = iter.next();
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testGeneralTable() {
        try {
            final GeneralTable<Integer, String> t0 = new GeneralTable<Integer, String>();
            t0.setValue( 3, 2, "23" );
            t0.setValue( 10, 1, "error" );
            t0.setValue( 10, 1, "110" );
            t0.setValue( 9, 1, "19" );
            t0.setValue( 1, 10, "101" );
            t0.setValue( 10, 10, "1010" );
            t0.setValue( 100, 10, "10100" );
            t0.setValue( 0, 0, "00" );
            if ( !t0.getValue( 3, 2 ).equals( "23" ) ) {
                return false;
            }
            if ( !t0.getValue( 10, 1 ).equals( "110" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 1, 10 ).equals( "101" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 10, 10 ).equals( "1010" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 100, 10 ).equals( "10100" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 9, 1 ).equals( "19" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 0, 0 ).equals( "00" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 49, 4 ).equals( "" ) ) {
                return false;
            }
            if ( !t0.getValueAsString( 22349, 3434344 ).equals( "" ) ) {
                return false;
            }
            final GeneralTable<String, String> t1 = new GeneralTable<String, String>();
            t1.setValue( "3", "2", "23" );
            t1.setValue( "10", "1", "error" );
            t1.setValue( "10", "1", "110" );
            t1.setValue( "9", "1", "19" );
            t1.setValue( "1", "10", "101" );
            t1.setValue( "10", "10", "1010" );
            t1.setValue( "100", "10", "10100" );
            t1.setValue( "0", "0", "00" );
            t1.setValue( "qwerty", "zxcvbnm", "asdef" );
            if ( !t1.getValue( "3", "2" ).equals( "23" ) ) {
                return false;
            }
            if ( !t1.getValue( "10", "1" ).equals( "110" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "1", "10" ).equals( "101" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "10", "10" ).equals( "1010" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "100", "10" ).equals( "10100" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "9", "1" ).equals( "19" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "0", "0" ).equals( "00" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "qwerty", "zxcvbnm" ).equals( "asdef" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "49", "4" ).equals( "" ) ) {
                return false;
            }
            if ( !t1.getValueAsString( "22349", "3434344" ).equals( "" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testGetDistance() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p1 = factory.create( "(((A:1,B:2,X:100)ab:3,C:4)abc:5,(D:7,(E:9,F:10)ef:8)def:6)r",
                                                 new NHXParser() )[ 0 ];
            final PhylogenyMethods pm = PhylogenyMethods.getInstance();
            if ( pm.calculateDistance( p1.getNode( "C" ), p1.getNode( "C" ) ) != 0 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "def" ), p1.getNode( "def" ) ) != 0 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "ef" ), p1.getNode( "ef" ) ) != 0 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "r" ), p1.getNode( "r" ) ) != 0 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "A" ), p1.getNode( "A" ) ) != 0 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "A" ), p1.getNode( "B" ) ) != 3 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "B" ), p1.getNode( "A" ) ) != 3 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "A" ), p1.getNode( "C" ) ) != 8 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "C" ), p1.getNode( "A" ) ) != 8 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "A" ), p1.getNode( "D" ) ) != 22 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "A" ), p1.getNode( "E" ) ) != 32 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "E" ), p1.getNode( "A" ) ) != 32 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "A" ), p1.getNode( "F" ) ) != 33 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "F" ), p1.getNode( "A" ) ) != 33 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "A" ), p1.getNode( "ab" ) ) != 1 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "ab" ), p1.getNode( "A" ) ) != 1 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "A" ), p1.getNode( "abc" ) ) != 4 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "abc" ), p1.getNode( "A" ) ) != 4 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "A" ), p1.getNode( "r" ) ) != 9 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "r" ), p1.getNode( "A" ) ) != 9 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "A" ), p1.getNode( "def" ) ) != 15 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "def" ), p1.getNode( "A" ) ) != 15 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "A" ), p1.getNode( "ef" ) ) != 23 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "ef" ), p1.getNode( "A" ) ) != 23 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "ef" ), p1.getNode( "def" ) ) != 8 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "def" ), p1.getNode( "ef" ) ) != 8 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "ef" ), p1.getNode( "r" ) ) != 14 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "ef" ), p1.getNode( "abc" ) ) != 19 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "ef" ), p1.getNode( "ab" ) ) != 22 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "ab" ), p1.getNode( "ef" ) ) != 22 ) {
                return false;
            }
            if ( pm.calculateDistance( p1.getNode( "def" ), p1.getNode( "abc" ) ) != 11 ) {
                return false;
            }
            final Phylogeny p2 = factory.create( "((A:4,B:5,C:6)abc:1,(D:7,E:8,F:9)def:2,(G:10,H:11,I:12)ghi:3)r",
                                                 new NHXParser() )[ 0 ];
            if ( pm.calculateDistance( p2.getNode( "A" ), p2.getNode( "B" ) ) != 9 ) {
                return false;
            }
            if ( pm.calculateDistance( p2.getNode( "A" ), p2.getNode( "C" ) ) != 10 ) {
                return false;
            }
            if ( pm.calculateDistance( p2.getNode( "A" ), p2.getNode( "D" ) ) != 14 ) {
                return false;
            }
            if ( pm.calculateDistance( p2.getNode( "A" ), p2.getNode( "ghi" ) ) != 8 ) {
                return false;
            }
            if ( pm.calculateDistance( p2.getNode( "A" ), p2.getNode( "I" ) ) != 20 ) {
                return false;
            }
            if ( pm.calculateDistance( p2.getNode( "G" ), p2.getNode( "ghi" ) ) != 10 ) {
                return false;
            }
            if ( pm.calculateDistance( p2.getNode( "r" ), p2.getNode( "r" ) ) != 0 ) {
                return false;
            }
            if ( pm.calculateDistance( p2.getNode( "r" ), p2.getNode( "G" ) ) != 13 ) {
                return false;
            }
            if ( pm.calculateDistance( p2.getNode( "G" ), p2.getNode( "r" ) ) != 13 ) {
                return false;
            }
            if ( pm.calculateDistance( p2.getNode( "G" ), p2.getNode( "H" ) ) != 21 ) {
                return false;
            }
            if ( pm.calculateDistance( p2.getNode( "G" ), p2.getNode( "I" ) ) != 22 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testGetLCA() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p1 = factory.create( "((((((A,B)ab,C)abc,D)abcd,E)abcde,F)abcdef,(G,H)gh)abcdefgh",
                                                 new NHXParser() )[ 0 ];
            final PhylogenyMethods pm = PhylogenyMethods.getInstance();
            final PhylogenyNode A = pm.getLCA( p1.getNode( "A" ), p1.getNode( "A" ) );
            if ( !A.getNodeName().equals( "A" ) ) {
                return false;
            }
            final PhylogenyNode gh = pm.getLCA( p1.getNode( "gh" ), p1.getNode( "gh" ) );
            if ( !gh.getNodeName().equals( "gh" ) ) {
                return false;
            }
            final PhylogenyNode ab = pm.getLCA( p1.getNode( "A" ), p1.getNode( "B" ) );
            if ( !ab.getNodeName().equals( "ab" ) ) {
                return false;
            }
            final PhylogenyNode ab2 = pm.getLCA( p1.getNode( "B" ), p1.getNode( "A" ) );
            if ( !ab2.getNodeName().equals( "ab" ) ) {
                return false;
            }
            final PhylogenyNode gh2 = pm.getLCA( p1.getNode( "H" ), p1.getNode( "G" ) );
            if ( !gh2.getNodeName().equals( "gh" ) ) {
                return false;
            }
            final PhylogenyNode gh3 = pm.getLCA( p1.getNode( "G" ), p1.getNode( "H" ) );
            if ( !gh3.getNodeName().equals( "gh" ) ) {
                return false;
            }
            final PhylogenyNode abc = pm.getLCA( p1.getNode( "C" ), p1.getNode( "A" ) );
            if ( !abc.getNodeName().equals( "abc" ) ) {
                return false;
            }
            final PhylogenyNode abc2 = pm.getLCA( p1.getNode( "A" ), p1.getNode( "C" ) );
            if ( !abc2.getNodeName().equals( "abc" ) ) {
                return false;
            }
            final PhylogenyNode abcd = pm.getLCA( p1.getNode( "A" ), p1.getNode( "D" ) );
            if ( !abcd.getNodeName().equals( "abcd" ) ) {
                return false;
            }
            final PhylogenyNode abcd2 = pm.getLCA( p1.getNode( "D" ), p1.getNode( "A" ) );
            if ( !abcd2.getNodeName().equals( "abcd" ) ) {
                return false;
            }
            final PhylogenyNode abcdef = pm.getLCA( p1.getNode( "A" ), p1.getNode( "F" ) );
            if ( !abcdef.getNodeName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode abcdef2 = pm.getLCA( p1.getNode( "F" ), p1.getNode( "A" ) );
            if ( !abcdef2.getNodeName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode abcdef3 = pm.getLCA( p1.getNode( "ab" ), p1.getNode( "F" ) );
            if ( !abcdef3.getNodeName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode abcdef4 = pm.getLCA( p1.getNode( "F" ), p1.getNode( "ab" ) );
            if ( !abcdef4.getNodeName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode abcde = pm.getLCA( p1.getNode( "A" ), p1.getNode( "E" ) );
            if ( !abcde.getNodeName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode abcde2 = pm.getLCA( p1.getNode( "E" ), p1.getNode( "A" ) );
            if ( !abcde2.getNodeName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode r = pm.getLCA( p1.getNode( "abcdefgh" ), p1.getNode( "abcdefgh" ) );
            if ( !r.getNodeName().equals( "abcdefgh" ) ) {
                return false;
            }
            final PhylogenyNode r2 = pm.getLCA( p1.getNode( "A" ), p1.getNode( "H" ) );
            if ( !r2.getNodeName().equals( "abcdefgh" ) ) {
                return false;
            }
            final PhylogenyNode r3 = pm.getLCA( p1.getNode( "H" ), p1.getNode( "A" ) );
            if ( !r3.getNodeName().equals( "abcdefgh" ) ) {
                return false;
            }
            final PhylogenyNode abcde3 = pm.getLCA( p1.getNode( "E" ), p1.getNode( "abcde" ) );
            if ( !abcde3.getNodeName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode abcde4 = pm.getLCA( p1.getNode( "abcde" ), p1.getNode( "E" ) );
            if ( !abcde4.getNodeName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode ab3 = pm.getLCA( p1.getNode( "ab" ), p1.getNode( "B" ) );
            if ( !ab3.getNodeName().equals( "ab" ) ) {
                return false;
            }
            final PhylogenyNode ab4 = pm.getLCA( p1.getNode( "B" ), p1.getNode( "ab" ) );
            if ( !ab4.getNodeName().equals( "ab" ) ) {
                return false;
            }
            final Phylogeny p2 = factory.create( "(a,b,(((c,d)cd,e)cde,f)cdef)r", new NHXParser() )[ 0 ];
            final PhylogenyNode cd = pm.getLCA( p2.getNode( "c" ), p2.getNode( "d" ) );
            if ( !cd.getNodeName().equals( "cd" ) ) {
                return false;
            }
            final PhylogenyNode cd2 = pm.getLCA( p2.getNode( "d" ), p2.getNode( "c" ) );
            if ( !cd2.getNodeName().equals( "cd" ) ) {
                return false;
            }
            final PhylogenyNode cde = pm.getLCA( p2.getNode( "c" ), p2.getNode( "e" ) );
            if ( !cde.getNodeName().equals( "cde" ) ) {
                return false;
            }
            final PhylogenyNode cde2 = pm.getLCA( p2.getNode( "e" ), p2.getNode( "c" ) );
            if ( !cde2.getNodeName().equals( "cde" ) ) {
                return false;
            }
            final PhylogenyNode cdef = pm.getLCA( p2.getNode( "c" ), p2.getNode( "f" ) );
            if ( !cdef.getNodeName().equals( "cdef" ) ) {
                return false;
            }
            final PhylogenyNode cdef2 = pm.getLCA( p2.getNode( "d" ), p2.getNode( "f" ) );
            if ( !cdef2.getNodeName().equals( "cdef" ) ) {
                return false;
            }
            final PhylogenyNode cdef3 = pm.getLCA( p2.getNode( "f" ), p2.getNode( "d" ) );
            if ( !cdef3.getNodeName().equals( "cdef" ) ) {
                return false;
            }
            final PhylogenyNode rt = pm.getLCA( p2.getNode( "c" ), p2.getNode( "a" ) );
            if ( !rt.getNodeName().equals( "r" ) ) {
                return false;
            }
            final Phylogeny p3 = factory
                    .create( "((((a,(b,c)bc)abc,(d,e)de)abcde,f)abcdef,(((g,h)gh,(i,j)ij)ghij,k)ghijk,l)",
                             new NHXParser() )[ 0 ];
            final PhylogenyNode bc_3 = pm.getLCA( p3.getNode( "b" ), p3.getNode( "c" ) );
            if ( !bc_3.getNodeName().equals( "bc" ) ) {
                return false;
            }
            final PhylogenyNode ac_3 = pm.getLCA( p3.getNode( "a" ), p3.getNode( "c" ) );
            if ( !ac_3.getNodeName().equals( "abc" ) ) {
                return false;
            }
            final PhylogenyNode ad_3 = pm.getLCA( p3.getNode( "a" ), p3.getNode( "d" ) );
            if ( !ad_3.getNodeName().equals( "abcde" ) ) {
                return false;
            }
            final PhylogenyNode af_3 = pm.getLCA( p3.getNode( "a" ), p3.getNode( "f" ) );
            if ( !af_3.getNodeName().equals( "abcdef" ) ) {
                return false;
            }
            final PhylogenyNode ag_3 = pm.getLCA( p3.getNode( "a" ), p3.getNode( "g" ) );
            if ( !ag_3.getNodeName().equals( "" ) ) {
                return false;
            }
            if ( !ag_3.isRoot() ) {
                return false;
            }
            final PhylogenyNode al_3 = pm.getLCA( p3.getNode( "a" ), p3.getNode( "l" ) );
            if ( !al_3.getNodeName().equals( "" ) ) {
                return false;
            }
            if ( !al_3.isRoot() ) {
                return false;
            }
            final PhylogenyNode kl_3 = pm.getLCA( p3.getNode( "k" ), p3.getNode( "l" ) );
            if ( !kl_3.getNodeName().equals( "" ) ) {
                return false;
            }
            if ( !kl_3.isRoot() ) {
                return false;
            }
            final PhylogenyNode fl_3 = pm.getLCA( p3.getNode( "f" ), p3.getNode( "l" ) );
            if ( !fl_3.getNodeName().equals( "" ) ) {
                return false;
            }
            if ( !fl_3.isRoot() ) {
                return false;
            }
            final PhylogenyNode gk_3 = pm.getLCA( p3.getNode( "g" ), p3.getNode( "k" ) );
            if ( !gk_3.getNodeName().equals( "ghijk" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testHmmscanOutputParser() {
        final String test_dir = Test.PATH_TO_TEST_DATA;
        try {
            final HmmscanPerDomainTableParser parser1 = new HmmscanPerDomainTableParser( new File( test_dir
                    + ForesterUtil.getFileSeparator() + "hmmscan30b3_output_1" ), "MONBR", INDIVIDUAL_SCORE_CUTOFF.NONE );
            parser1.parse();
            final HmmscanPerDomainTableParser parser2 = new HmmscanPerDomainTableParser( new File( test_dir
                    + ForesterUtil.getFileSeparator() + "hmmscan30b3_output_2" ), "MONBR", INDIVIDUAL_SCORE_CUTOFF.NONE );
            final List<Protein> domain_collections = parser2.parse();
            if ( parser2.getProteinsEncountered() != 4 ) {
                return false;
            }
            if ( domain_collections.size() != 4 ) {
                return false;
            }
            if ( parser2.getDomainsEncountered() != 69 ) {
                return false;
            }
            if ( parser2.getDomainsIgnoredDueToDuf() != 0 ) {
                return false;
            }
            if ( parser2.getDomainsIgnoredDueToEval() != 0 ) {
                return false;
            }
            final Protein p1 = domain_collections.get( 0 );
            if ( p1.getNumberOfProteinDomains() != 15 ) {
                return false;
            }
            final Protein p4 = domain_collections.get( 3 );
            if ( p4.getNumberOfProteinDomains() != 1 ) {
                return false;
            }
            if ( !p4.getProteinDomain( 0 ).getDomainId().toString().equals( "DNA_pol_B_new" ) ) {
                return false;
            }
            if ( p4.getProteinDomain( 0 ).getFrom() != 51 ) {
                return false;
            }
            if ( p4.getProteinDomain( 0 ).getTo() != 395 ) {
                return false;
            }
            if ( !Test.isEqual( p4.getProteinDomain( 0 ).getPerDomainEvalue(), 1.2e-39 ) ) {
                return false;
            }
            if ( !Test.isEqual( p4.getProteinDomain( 0 ).getPerDomainScore(), 135.7 ) ) {
                return false;
            }
            if ( !Test.isEqual( p4.getProteinDomain( 0 ).getPerSequenceEvalue(), 8.3e-40 ) ) {
                return false;
            }
            if ( !Test.isEqual( p4.getProteinDomain( 0 ).getPerSequenceScore(), 136.3 ) ) {
                return false;
            }
            if ( !Test.isEqual( p4.getProteinDomain( 0 ).getNumber(), 1 ) ) {
                return false;
            }
            if ( !Test.isEqual( p4.getProteinDomain( 0 ).getTotalCount(), 1 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testLastExternalNodeMethods() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final char[] a0 = { '(', '(', 'A', ',', 'B', ')', ',', '(', 'C', ',', 'D', ')', ')', };
            final Phylogeny t0 = factory.create( a0, new NHXParser() )[ 0 ];
            final PhylogenyNode n1 = t0.getNode( "A" );
            if ( n1.isLastExternalNode() ) {
                return false;
            }
            final PhylogenyNode n2 = t0.getNode( "B" );
            if ( n2.isLastExternalNode() ) {
                return false;
            }
            final PhylogenyNode n3 = t0.getNode( "C" );
            if ( n3.isLastExternalNode() ) {
                return false;
            }
            final PhylogenyNode n4 = t0.getNode( "D" );
            if ( !n4.isLastExternalNode() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testLevelOrderIterator() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0 = factory.create( "((A,B)ab,(C,D)cd)r", new NHXParser() )[ 0 ];
            PhylogenyNodeIterator it0;
            for( it0 = t0.iteratorLevelOrder(); it0.hasNext(); ) {
                it0.next();
            }
            for( it0.reset(); it0.hasNext(); ) {
                it0.next();
            }
            final PhylogenyNodeIterator it = t0.iteratorLevelOrder();
            if ( !it.next().getNodeName().equals( "r" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "ab" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "cd" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "A" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "B" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "C" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "D" ) ) {
                return false;
            }
            if ( it.hasNext() ) {
                return false;
            }
            final Phylogeny t2 = factory.create( "(((1,2,(a,(X,Y,Z)b)3,4,5,6)A,B,C)abc,(D,E,(f1,(f21)f2,f3)F,G)defg)r",
                                                 new NHXParser() )[ 0 ];
            PhylogenyNodeIterator it2;
            for( it2 = t2.iteratorLevelOrder(); it2.hasNext(); ) {
                it2.next();
            }
            for( it2.reset(); it2.hasNext(); ) {
                it2.next();
            }
            final PhylogenyNodeIterator it3 = t2.iteratorLevelOrder();
            if ( !it3.next().getNodeName().equals( "r" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "abc" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "defg" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "A" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "B" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "C" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "D" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "E" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "F" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "G" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "1" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "2" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "3" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "4" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "5" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "6" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "f1" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "f2" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "f3" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "a" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "b" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "f21" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "X" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "Y" ) ) {
                return false;
            }
            if ( !it3.next().getNodeName().equals( "Z" ) ) {
                return false;
            }
            if ( it3.hasNext() ) {
                return false;
            }
            final Phylogeny t4 = factory.create( "((((D)C)B)A)r", new NHXParser() )[ 0 ];
            PhylogenyNodeIterator it4;
            for( it4 = t4.iteratorLevelOrder(); it4.hasNext(); ) {
                it4.next();
            }
            for( it4.reset(); it4.hasNext(); ) {
                it4.next();
            }
            final PhylogenyNodeIterator it5 = t4.iteratorLevelOrder();
            if ( !it5.next().getNodeName().equals( "r" ) ) {
                return false;
            }
            if ( !it5.next().getNodeName().equals( "A" ) ) {
                return false;
            }
            if ( !it5.next().getNodeName().equals( "B" ) ) {
                return false;
            }
            if ( !it5.next().getNodeName().equals( "C" ) ) {
                return false;
            }
            if ( !it5.next().getNodeName().equals( "D" ) ) {
                return false;
            }
            final Phylogeny t5 = factory.create( "A", new NHXParser() )[ 0 ];
            PhylogenyNodeIterator it6;
            for( it6 = t5.iteratorLevelOrder(); it6.hasNext(); ) {
                it6.next();
            }
            for( it6.reset(); it6.hasNext(); ) {
                it6.next();
            }
            final PhylogenyNodeIterator it7 = t5.iteratorLevelOrder();
            if ( !it7.next().getNodeName().equals( "A" ) ) {
                return false;
            }
            if ( it.hasNext() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testMidpointrooting() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t1 = factory.create( "((A:1,B:2)AB:1[&&NHX:B=55],(C:3,D:4)CD:3[&&NHX:B=10])ABCD:0.5",
                                                 new NHXParser() )[ 0 ];
            if ( !t1.isRooted() ) {
                return false;
            }
            PhylogenyMethods.midpointRoot( t1 );
            if ( !isEqual( t1.getNode( "A" ).getDistanceToParent(), 1 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "B" ).getDistanceToParent(), 2 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "C" ).getDistanceToParent(), 3 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "D" ).getDistanceToParent(), 4 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "CD" ).getDistanceToParent(), 1 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "AB" ).getDistanceToParent(), 3 ) ) {
                return false;
            }
            t1.reRoot( t1.getNode( "A" ) );
            PhylogenyMethods.midpointRoot( t1 );
            if ( !isEqual( t1.getNode( "A" ).getDistanceToParent(), 1 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "B" ).getDistanceToParent(), 2 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "C" ).getDistanceToParent(), 3 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "D" ).getDistanceToParent(), 4 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "CD" ).getDistanceToParent(), 1 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "AB" ).getDistanceToParent(), 3 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNexusCharactersParsing() {
        try {
            final NexusCharactersParser parser = new NexusCharactersParser();
            parser.setSource( new File( Test.PATH_TO_TEST_DATA + "nexus_test_7.nex" ) );
            parser.parse();
            String[] labels = parser.getCharStateLabels();
            if ( labels.length != 7 ) {
                return false;
            }
            if ( !labels[ 0 ].equals( "14-3-3" ) ) {
                return false;
            }
            if ( !labels[ 1 ].equals( "2-Hacid_dh" ) ) {
                return false;
            }
            if ( !labels[ 2 ].equals( "2-Hacid_dh_C" ) ) {
                return false;
            }
            if ( !labels[ 3 ].equals( "2-oxoacid_dh" ) ) {
                return false;
            }
            if ( !labels[ 4 ].equals( "2OG-FeII_Oxy" ) ) {
                return false;
            }
            if ( !labels[ 5 ].equals( "3-HAO" ) ) {
                return false;
            }
            if ( !labels[ 6 ].equals( "3_5_exonuc" ) ) {
                return false;
            }
            parser.setSource( new File( Test.PATH_TO_TEST_DATA + "nexus_test_8.nex" ) );
            parser.parse();
            labels = parser.getCharStateLabels();
            if ( labels.length != 7 ) {
                return false;
            }
            if ( !labels[ 0 ].equals( "14-3-3" ) ) {
                return false;
            }
            if ( !labels[ 1 ].equals( "2-Hacid_dh" ) ) {
                return false;
            }
            if ( !labels[ 2 ].equals( "2-Hacid_dh_C" ) ) {
                return false;
            }
            if ( !labels[ 3 ].equals( "2-oxoacid_dh" ) ) {
                return false;
            }
            if ( !labels[ 4 ].equals( "2OG-FeII_Oxy" ) ) {
                return false;
            }
            if ( !labels[ 5 ].equals( "3-HAO" ) ) {
                return false;
            }
            if ( !labels[ 6 ].equals( "3_5_exonuc" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNexusMatrixParsing() {
        try {
            final NexusBinaryStatesMatrixParser parser = new NexusBinaryStatesMatrixParser();
            parser.setSource( new File( Test.PATH_TO_TEST_DATA + "nexus_test_9.nex" ) );
            parser.parse();
            final CharacterStateMatrix<BinaryStates> m = parser.getMatrix();
            if ( m.getNumberOfCharacters() != 9 ) {
                return false;
            }
            if ( m.getNumberOfIdentifiers() != 5 ) {
                return false;
            }
            if ( m.getState( 0, 0 ) != BinaryStates.PRESENT ) {
                return false;
            }
            if ( m.getState( 0, 1 ) != BinaryStates.ABSENT ) {
                return false;
            }
            if ( m.getState( 1, 0 ) != BinaryStates.PRESENT ) {
                return false;
            }
            if ( m.getState( 2, 0 ) != BinaryStates.ABSENT ) {
                return false;
            }
            if ( m.getState( 4, 8 ) != BinaryStates.PRESENT ) {
                return false;
            }
            if ( !m.getIdentifier( 0 ).equals( "MOUSE" ) ) {
                return false;
            }
            if ( !m.getIdentifier( 4 ).equals( "ARATH" ) ) {
                return false;
            }
            //            if ( labels.length != 7 ) {
            //                return false;
            //            }
            //            if ( !labels[ 0 ].equals( "14-3-3" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 1 ].equals( "2-Hacid_dh" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 2 ].equals( "2-Hacid_dh_C" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 3 ].equals( "2-oxoacid_dh" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 4 ].equals( "2OG-FeII_Oxy" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 5 ].equals( "3-HAO" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 6 ].equals( "3_5_exonuc" ) ) {
            //                return false;
            //            }
            //            parser.setSource( new File( Test.PATH_TO_TEST_DATA + "nexus_test_8.nex" ) );
            //            parser.parse();
            //            labels = parser.getCharStateLabels();
            //            if ( labels.length != 7 ) {
            //                return false;
            //            }
            //            if ( !labels[ 0 ].equals( "14-3-3" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 1 ].equals( "2-Hacid_dh" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 2 ].equals( "2-Hacid_dh_C" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 3 ].equals( "2-oxoacid_dh" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 4 ].equals( "2OG-FeII_Oxy" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 5 ].equals( "3-HAO" ) ) {
            //                return false;
            //            }
            //            if ( !labels[ 6 ].equals( "3_5_exonuc" ) ) {
            //                return false;
            //            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNexusTreeParsing() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final NexusPhylogeniesParser parser = new NexusPhylogeniesParser();
            Phylogeny[] phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_1.nex", parser );
            if ( phylogenies.length != 1 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 25 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "" ) ) {
                return false;
            }
            phylogenies = null;
            phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_2.nex", parser );
            if ( phylogenies.length != 1 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "name" ) ) {
                return false;
            }
            phylogenies = null;
            phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_3.nex", parser );
            if ( phylogenies.length != 1 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "" ) ) {
                return false;
            }
            if ( phylogenies[ 0 ].isRooted() ) {
                return false;
            }
            phylogenies = null;
            phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_4.nex", parser );
            if ( phylogenies.length != 18 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "tree 0" ) ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getName().equals( "tree 1" ) ) {
                return false;
            }
            if ( phylogenies[ 1 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( phylogenies[ 2 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( phylogenies[ 3 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( phylogenies[ 4 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( phylogenies[ 5 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( phylogenies[ 6 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( phylogenies[ 7 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 8 ].getName().equals( "tree 8" ) ) {
                return false;
            }
            if ( phylogenies[ 8 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 8 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 9 ].getName().equals( "tree 9" ) ) {
                return false;
            }
            if ( !phylogenies[ 9 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 9 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 10 ].getName().equals( "tree 10" ) ) {
                return false;
            }
            if ( !phylogenies[ 10 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 10 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 11 ].getName().equals( "tree 11" ) ) {
                return false;
            }
            if ( phylogenies[ 11 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 11 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 12 ].getName().equals( "tree 12" ) ) {
                return false;
            }
            if ( !phylogenies[ 12 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 12 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 13 ].getName().equals( "tree 13" ) ) {
                return false;
            }
            if ( !phylogenies[ 13 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 13 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 14 ].getName().equals( "tree 14" ) ) {
                return false;
            }
            if ( !phylogenies[ 14 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 14 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phylogenies[ 15 ].getName().equals( "tree 15" ) ) {
                return false;
            }
            if ( phylogenies[ 15 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 15 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phylogenies[ 16 ].getName().equals( "tree 16" ) ) {
                return false;
            }
            if ( !phylogenies[ 16 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 16 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
            if ( !phylogenies[ 17 ].getName().equals( "tree 17" ) ) {
                return false;
            }
            if ( phylogenies[ 17 ].isRooted() ) {
                return false;
            }
            if ( phylogenies[ 17 ].getNumberOfExternalNodes() != 10 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNexusTreeParsingTranslating() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final NexusPhylogeniesParser parser = new NexusPhylogeniesParser();
            Phylogeny[] phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_5.nex", parser );
            if ( phylogenies.length != 1 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "Tree0" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNodeName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNextExternalNode().getNodeName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getNodeName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
            phylogenies = null;
            phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_6.nex", parser );
            if ( phylogenies.length != 3 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "Tree0" ) ) {
                return false;
            }
            if ( phylogenies[ 0 ].isRooted() ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNodeName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNextExternalNode().getNodeName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getNodeName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
            if ( phylogenies[ 1 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getName().equals( "Tree1" ) ) {
                return false;
            }
            if ( phylogenies[ 1 ].isRooted() ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getFirstExternalNode().getNodeName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getFirstExternalNode().getNextExternalNode().getNodeName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getNodeName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
            if ( phylogenies[ 2 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getName().equals( "Tree2" ) ) {
                return false;
            }
            if ( !phylogenies[ 2 ].isRooted() ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getFirstExternalNode().getNodeName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getFirstExternalNode().getNextExternalNode().getNodeName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getNodeName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
            phylogenies = null;
            phylogenies = factory.create( Test.PATH_TO_TEST_DATA + "nexus_test_7.nex", parser );
            if ( phylogenies.length != 3 ) {
                return false;
            }
            if ( phylogenies[ 0 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getName().equals( "Tree0" ) ) {
                return false;
            }
            if ( phylogenies[ 0 ].isRooted() ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNodeName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNextExternalNode().getNodeName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 0 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getNodeName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
            if ( phylogenies[ 1 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getName().equals( "Tree1" ) ) {
                return false;
            }
            if ( phylogenies[ 1 ].isRooted() ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getFirstExternalNode().getNodeName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getFirstExternalNode().getNextExternalNode().getNodeName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 1 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getNodeName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
            if ( phylogenies[ 2 ].getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getName().equals( "Tree2" ) ) {
                return false;
            }
            if ( !phylogenies[ 2 ].isRooted() ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getFirstExternalNode().getNodeName().equals( "Scarabaeus" ) ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getFirstExternalNode().getNextExternalNode().getNodeName().equals( "Drosophila" ) ) {
                return false;
            }
            if ( !phylogenies[ 2 ].getFirstExternalNode().getNextExternalNode().getNextExternalNode().getNodeName()
                    .equals( "Aranaeus" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNHParsing() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p1 = factory.create( "(A,B1)", new NHXParser() )[ 0 ];
            if ( !p1.toNewHampshireX().equals( "(A,B1)" ) ) {
                return false;
            }
            final NHXParser nhxp = new NHXParser();
            nhxp.setTaxonomyExtraction( TAXONOMY_EXTRACTION.NO );
            nhxp.setReplaceUnderscores( true );
            final Phylogeny uc0 = factory.create( "(A__A_,_B_B)", nhxp )[ 0 ];
            if ( !uc0.getRoot().getChildNode( 0 ).getNodeName().equals( "A A " ) ) {
                return false;
            }
            if ( !uc0.getRoot().getChildNode( 1 ).getNodeName().equals( " B B" ) ) {
                return false;
            }
            final Phylogeny p1b = factory
                    .create( "   \n  \t  \b   \r \f   ; (  \n  \t  \b   \r \f; A ;  \n  \t  \b   \r \f,  \n  \t  \b   \r \f; B ;   \n  \t  \b   \r \f 1  \n  \t  \b   \r \f ;  \n  \t  \b   \r \f );;;;; \n  \t  \b   \r \f;;;  \n  \t  \b   \r \f ",
                             new NHXParser() )[ 0 ];
            if ( !p1b.toNewHampshireX().equals( "(_A_,_B_1_)" ) ) {
                return false;
            }
            final Phylogeny p2 = factory.create( new StringBuffer( "(A,B2)" ), new NHXParser() )[ 0 ];
            final Phylogeny p3 = factory.create( new char[] { '(', 'A', ',', 'B', '3', ')' }, new NHXParser() )[ 0 ];
            final Phylogeny p4 = factory.create( "(A,B4);", new NHXParser() )[ 0 ];
            final Phylogeny p5 = factory.create( new StringBuffer( "(A,B5);" ), new NHXParser() )[ 0 ];
            final Phylogeny[] p7 = factory.create( "(A,B7);(C,D7)", new NHXParser() );
            final Phylogeny[] p8 = factory.create( "(A,B8) (C,D8)", new NHXParser() );
            final Phylogeny[] p9 = factory.create( "(A,B9)\n(C,D9)", new NHXParser() );
            final Phylogeny[] p10 = factory.create( "(A,B10);(C,D10);", new NHXParser() );
            final Phylogeny[] p11 = factory.create( "(A,B11);(C,D11) (E,F11)\t(G,H11)", new NHXParser() );
            final Phylogeny[] p12 = factory.create( "(A,B12) (C,D12) (E,F12) (G,H12)", new NHXParser() );
            final Phylogeny[] p13 = factory.create( " ; (;A; , ; B ; 1  3 ; \n)\t ( \n ;"
                                                            + " C ; ,; D;13;);;;;;;(;E;,;F;13 ;) ; "
                                                            + "; ; ( \t\n\r\b; G ;, ;H ;1 3; )  ;  ;   ;",
                                                    new NHXParser() );
            if ( !p13[ 0 ].toNewHampshireX().equals( "(_A_,_B_13_)" ) ) {
                return false;
            }
            if ( !p13[ 1 ].toNewHampshireX().equals( "(_C_,_D_13_)" ) ) {
                return false;
            }
            if ( !p13[ 2 ].toNewHampshireX().equals( "(_E_,_F_13_)" ) ) {
                return false;
            }
            if ( !p13[ 3 ].toNewHampshireX().equals( "(_G_,_H_13_)" ) ) {
                return false;
            }
            final Phylogeny[] p14 = factory.create( "(A,B14)ab", new NHXParser() );
            final Phylogeny[] p15 = factory.create( "(A,B15)ab;", new NHXParser() );
            final String p16_S = "((A,B),C)";
            final Phylogeny[] p16 = factory.create( p16_S, new NHXParser() );
            if ( !p16[ 0 ].toNewHampshireX().equals( p16_S ) ) {
                return false;
            }
            final String p17_S = "(C,(A,B))";
            final Phylogeny[] p17 = factory.create( p17_S, new NHXParser() );
            if ( !p17[ 0 ].toNewHampshireX().equals( p17_S ) ) {
                return false;
            }
            final String p18_S = "((A,B),(C,D))";
            final Phylogeny[] p18 = factory.create( p18_S, new NHXParser() );
            if ( !p18[ 0 ].toNewHampshireX().equals( p18_S ) ) {
                return false;
            }
            final String p19_S = "(((A,B),C),D)";
            final Phylogeny[] p19 = factory.create( p19_S, new NHXParser() );
            if ( !p19[ 0 ].toNewHampshireX().equals( p19_S ) ) {
                return false;
            }
            final String p20_S = "(A,(B,(C,D)))";
            final Phylogeny[] p20 = factory.create( p20_S, new NHXParser() );
            if ( !p20[ 0 ].toNewHampshireX().equals( p20_S ) ) {
                return false;
            }
            final String p21_S = "(A,(B,(C,(D,E))))";
            final Phylogeny[] p21 = factory.create( p21_S, new NHXParser() );
            if ( !p21[ 0 ].toNewHampshireX().equals( p21_S ) ) {
                return false;
            }
            final String p22_S = "((((A,B),C),D),E)";
            final Phylogeny[] p22 = factory.create( p22_S, new NHXParser() );
            if ( !p22[ 0 ].toNewHampshireX().equals( p22_S ) ) {
                return false;
            }
            final String p23_S = "(A,(B,(C,(D,E)de)cde)bcde)abcde";
            final Phylogeny[] p23 = factory.create( p23_S, new NHXParser() );
            if ( !p23[ 0 ].toNewHampshireX().equals( p23_S ) ) {
                return false;
            }
            final String p24_S = "((((A,B)ab,C)abc,D)abcd,E)abcde";
            final Phylogeny[] p24 = factory.create( p24_S, new NHXParser() );
            if ( !p24[ 0 ].toNewHampshireX().equals( p24_S ) ) {
                return false;
            }
            final String p241_S1 = "(A,(B,(C,(D,E)de)cde)bcde)abcde";
            final String p241_S2 = "((((A,B)ab,C)abc,D)abcd,E)abcde";
            final Phylogeny[] p241 = factory.create( p241_S1 + p241_S2, new NHXParser() );
            if ( !p241[ 0 ].toNewHampshireX().equals( p241_S1 ) ) {
                return false;
            }
            if ( !p241[ 1 ].toNewHampshireX().equals( p241_S2 ) ) {
                return false;
            }
            final String p25_S = "((((((((((((((A,B)ab,C)abc,D)abcd,E)"
                    + "abcde,(B,(C,(D,E)de)cde)bcde)abcde,(B,((A,(B,(C,(D,"
                    + "E)de)cde)bcde)abcde,(D,E)de)cde)bcde)abcde,B)ab,C)"
                    + "abc,((((A,B)ab,C)abc,D)abcd,E)abcde)abcd,E)abcde,"
                    + "((((A,((((((((A,B)ab,C)abc,((((A,B)ab,C)abc,D)abcd,"
                    + "E)abcde)abcd,E)abcde,((((A,B)ab,C)abc,D)abcd,E)abcde)"
                    + "ab,C)abc,((((A,B)ab,C)abc,D)abcd,E)abcde)abcd,E)abcde"
                    + ")ab,C)abc,D)abcd,E)abcde)ab,C)abc,((((A,B)ab,C)abc,D)" + "abcd,E)abcde)abcd,E)abcde";
            final Phylogeny[] p25 = factory.create( p25_S, new NHXParser() );
            if ( !p25[ 0 ].toNewHampshireX().equals( p25_S ) ) {
                return false;
            }
            final String p26_S = "(A,B)ab";
            final Phylogeny[] p26 = factory.create( p26_S, new NHXParser() );
            if ( !p26[ 0 ].toNewHampshireX().equals( p26_S ) ) {
                return false;
            }
            final String p27_S = "((((A,B)ab,C)abc,D)abcd,E)abcde";
            final Phylogeny[] p27 = factory.create( new File( Test.PATH_TO_TEST_DATA + "phylogeny27.nhx" ),
                                                    new NHXParser() );
            if ( !p27[ 0 ].toNewHampshireX().equals( p27_S ) ) {
                return false;
            }
            final String p28_S1 = "((((A,B)ab,C)abc,D)abcd,E)abcde";
            final String p28_S2 = "(A,(B,(C,(D,E)de)cde)bcde)abcde";
            final String p28_S3 = "(A,B)ab";
            final String p28_S4 = "((((A,B),C),D),;E;)";
            final Phylogeny[] p28 = factory.create( new File( Test.PATH_TO_TEST_DATA + "phylogeny28.nhx" ),
                                                    new NHXParser() );
            if ( !p28[ 0 ].toNewHampshireX().equals( p28_S1 ) ) {
                return false;
            }
            if ( !p28[ 1 ].toNewHampshireX().equals( p28_S2 ) ) {
                return false;
            }
            if ( !p28[ 2 ].toNewHampshireX().equals( p28_S3 ) ) {
                return false;
            }
            if ( !p28[ 3 ].toNewHampshireX().equals( "((((A,B),C),D),_E_)" ) ) {
                return false;
            }
            final String p29_S = "((((A:0.01,B:0.684)ab:0.345,C:0.3451)abc:0.3451,D:1.5)abcd:0.134,E:0.32)abcde:0.1345";
            final Phylogeny[] p29 = factory.create( p29_S, new NHXParser() );
            if ( !p29[ 0 ].toNewHampshireX().equals( p29_S ) ) {
                return false;
            }
            final String p30_S = "((((A:0.01,B:0.02):0.93,C:0.04):0.05,D:1.4):0.06,E):0.72";
            final Phylogeny[] p30 = factory.create( p30_S, new NHXParser() );
            if ( !p30[ 0 ].toNewHampshireX().equals( p30_S ) ) {
                return false;
            }
            final String p32_S = " ;   ; 	\n  \t  \b   \f  \r  ;;;;;; ";
            final Phylogeny[] p32 = factory.create( p32_S, new NHXParser() );
            if ( ( p32.length != 1 ) || !p32[ 0 ].isEmpty() ) {
                return false;
            }
            final String p33_S = "A";
            final Phylogeny[] p33 = factory.create( p33_S, new NHXParser() );
            if ( !p33[ 0 ].toNewHampshireX().equals( p33_S ) ) {
                return false;
            }
            final String p34_S = "B;";
            final Phylogeny[] p34 = factory.create( p34_S, new NHXParser() );
            if ( !p34[ 0 ].toNewHampshireX().equals( "B" ) ) {
                return false;
            }
            final String p35_S = "B:0.2";
            final Phylogeny[] p35 = factory.create( p35_S, new NHXParser() );
            if ( !p35[ 0 ].toNewHampshireX().equals( p35_S ) ) {
                return false;
            }
            final String p36_S = "(A)";
            final Phylogeny[] p36 = factory.create( p36_S, new NHXParser() );
            if ( !p36[ 0 ].toNewHampshireX().equals( p36_S ) ) {
                return false;
            }
            final String p37_S = "((A))";
            final Phylogeny[] p37 = factory.create( p37_S, new NHXParser() );
            if ( !p37[ 0 ].toNewHampshireX().equals( p37_S ) ) {
                return false;
            }
            final String p38_S = "(((((((A:0.2):0.2):0.3):0.4):0.5):0.6):0.7):0.8";
            final Phylogeny[] p38 = factory.create( p38_S, new NHXParser() );
            if ( !p38[ 0 ].toNewHampshireX().equals( p38_S ) ) {
                return false;
            }
            final String p39_S = "(((B,((((A:0.2):0.2):0.3):0.4):0.5):0.6):0.7):0.8";
            final Phylogeny[] p39 = factory.create( p39_S, new NHXParser() );
            if ( !p39[ 0 ].toNewHampshireX().equals( p39_S ) ) {
                return false;
            }
            final String p40_S = "(A,B,C)";
            final Phylogeny[] p40 = factory.create( p40_S, new NHXParser() );
            if ( !p40[ 0 ].toNewHampshireX().equals( p40_S ) ) {
                return false;
            }
            final String p41_S = "(A,B,C,D,E,F,G,H,I,J,K)";
            final Phylogeny[] p41 = factory.create( p41_S, new NHXParser() );
            if ( !p41[ 0 ].toNewHampshireX().equals( p41_S ) ) {
                return false;
            }
            final String p42_S = "(A,B,(X,Y,Z),D,E,F,G,H,I,J,K)";
            final Phylogeny[] p42 = factory.create( p42_S, new NHXParser() );
            if ( !p42[ 0 ].toNewHampshireX().equals( p42_S ) ) {
                return false;
            }
            final String p43_S = "(A,B,C,(AA,BB,CC,(CCC,DDD,EEE,(FFFF,GGGG)x)y,DD,EE,FF,GG,HH),D,E,(EE,FF),F,G,H,(((((5)4)3)2)1),I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,(XX,(YY)),Y,Z)";
            final Phylogeny[] p43 = factory.create( p43_S, new NHXParser() );
            if ( !p43[ 0 ].toNewHampshireX().equals( p43_S ) ) {
                return false;
            }
            final String p44_S = "(((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)))";
            final Phylogeny[] p44 = factory.create( p44_S, new NHXParser() );
            if ( !p44[ 0 ].toNewHampshireX().equals( p44_S ) ) {
                return false;
            }
            final String p45_S = "((((((((((A))))))))),(((((((((B))))))))),(((((((((C))))))))))";
            final Phylogeny[] p45 = factory.create( p45_S, new NHXParser() );
            if ( !p45[ 0 ].toNewHampshireX().equals( p45_S ) ) {
                return false;
            }
            final String p46_S = "";
            final Phylogeny[] p46 = factory.create( p46_S, new NHXParser() );
            if ( ( p46.length != 1 ) || !p46[ 0 ].isEmpty() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNHXconversion() {
        try {
            final PhylogenyNode n1 = new PhylogenyNode();
            final PhylogenyNode n2 = new PhylogenyNode( "" );
            final PhylogenyNode n3 = new PhylogenyNode( "n3" );
            final PhylogenyNode n4 = new PhylogenyNode( "n4:0.01" );
            final PhylogenyNode n5 = new PhylogenyNode( "n5:0.1[&&NHX:S=Ecoli:E=1.1.1.1:D=Y:Co=Y:B=56:T=1:W=2:C=10.20.30:XN=S=tag1=value1=unit1]" );
            final PhylogenyNode n6 = new PhylogenyNode( "n6:0.000001[&&NHX:S=Ecoli:E=1.1.1.1:D=N:Co=N:B=100:T=1:W=2:C=0.0.0:XN=B=bool_tag=T]" );
            if ( !n1.toNewHampshireX().equals( "" ) ) {
                return false;
            }
            if ( !n2.toNewHampshireX().equals( "" ) ) {
                return false;
            }
            if ( !n3.toNewHampshireX().equals( "n3" ) ) {
                return false;
            }
            if ( !n4.toNewHampshireX().equals( "n4:0.01" ) ) {
                return false;
            }
            if ( !n5.toNewHampshireX()
                    .equals( "n5:0.1[&&NHX:T=1:S=Ecoli:D=Y:XN=S=tag1=value1=unit1:B=56.0:W=2.0:C=10.20.30]" ) ) {
                return false;
            }
            if ( !n6.toNewHampshireX()
                    .equals( "n6:1.0E-6[&&NHX:T=1:S=Ecoli:D=N:XN=B=bool_tag=T:B=100.0:W=2.0:C=0.0.0]" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNHXNodeParsing() {
        try {
            final PhylogenyNode n1 = new PhylogenyNode();
            final PhylogenyNode n2 = new PhylogenyNode( "" );
            final PhylogenyNode n3 = new PhylogenyNode( "n3" );
            final PhylogenyNode n4 = new PhylogenyNode( "n4:0.01" );
            final PhylogenyNode n5 = new PhylogenyNode( "n5:0.1[&&NHX:S=Ecoli:E=1.1.1.1:D=Y:B=56:T=1:On=22:SOn=33:SNn=44:W=2:C=10.20.30:XN=S=tag1=value1=unit1:XN=S=tag3=value3=unit3]" );
            if ( !n3.getNodeName().equals( "n3" ) ) {
                return false;
            }
            if ( n3.getDistanceToParent() != PhylogenyNode.DISTANCE_DEFAULT ) {
                return false;
            }
            if ( n3.isDuplication() ) {
                return false;
            }
            if ( n3.isHasAssignedEvent() ) {
                return false;
            }
            if ( PhylogenyMethods.getBranchWidthValue( n3 ) != BranchWidth.BRANCH_WIDTH_DEFAULT_VALUE ) {
                return false;
            }
            if ( !n4.getNodeName().equals( "n4" ) ) {
                return false;
            }
            if ( n4.getDistanceToParent() != 0.01 ) {
                return false;
            }
            if ( !n5.getNodeName().equals( "n5" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( n5 ) != 56 ) {
                return false;
            }
            if ( n5.getDistanceToParent() != 0.1 ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n5 ).equals( "Ecoli" ) ) {
                return false;
            }
            if ( !n5.isDuplication() ) {
                return false;
            }
            if ( !n5.isHasAssignedEvent() ) {
                return false;
            }
            if ( PhylogenyMethods.getBranchWidthValue( n5 ) != 2 ) {
                return false;
            }
            if ( n5.getNodeData().getProperties().getPropertyRefs().length != 2 ) {
                return false;
            }
            final PhylogenyNode n8 = new PhylogenyNode( "n8_ECOLI/12:0.01", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n8.getNodeName().equals( "n8_ECOLI/12" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n8 ).equals( "ECOLI" ) ) {
                return false;
            }
            final PhylogenyNode n9 = new PhylogenyNode( "n9_ECOLI/12=12:0.01", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n9.getNodeName().equals( "n9_ECOLI/12=12" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n9 ).equals( "ECOLI" ) ) {
                return false;
            }
            final PhylogenyNode n10 = new PhylogenyNode( "n10.ECOLI", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n10.getNodeName().equals( "n10.ECOLI" ) ) {
                return false;
            }
            final PhylogenyNode n20 = new PhylogenyNode( "n20_ECOLI/1-2", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n20.getNodeName().equals( "n20_ECOLI/1-2" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n20 ).equals( "ECOLI" ) ) {
                return false;
            }
            final PhylogenyNode n20x = new PhylogenyNode( "n20_ECOL1/1-2", TAXONOMY_EXTRACTION.YES );
            if ( !n20x.getNodeName().equals( "n20_ECOL1/1-2" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n20x ).equals( "ECOL1" ) ) {
                return false;
            }
            final PhylogenyNode n20xx = new PhylogenyNode( "n20_eCOL1/1-2", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n20xx.getNodeName().equals( "n20_eCOL1/1-2" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n20xx ).length() > 0 ) {
                return false;
            }
            final PhylogenyNode n20xxx = new PhylogenyNode( "n20_ecoli/1-2", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n20xxx.getNodeName().equals( "n20_ecoli/1-2" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n20xxx ).length() > 0 ) {
                return false;
            }
            final PhylogenyNode n20xxxx = new PhylogenyNode( "n20_Ecoli/1-2", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n20xxxx.getNodeName().equals( "n20_Ecoli/1-2" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n20xxxx ).length() > 0 ) {
                return false;
            }
            final PhylogenyNode n21 = new PhylogenyNode( "n21_PIG", TAXONOMY_EXTRACTION.YES );
            if ( !n21.getNodeName().equals( "n21_PIG" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n21 ).equals( "PIG" ) ) {
                return false;
            }
            final PhylogenyNode n21x = new PhylogenyNode( "n21_PIG", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n21x.getNodeName().equals( "n21_PIG" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n21x ).length() > 0 ) {
                return false;
            }
            final PhylogenyNode n22 = new PhylogenyNode( "n22/PIG", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n22.getNodeName().equals( "n22/PIG" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n22 ).length() > 0 ) {
                return false;
            }
            final PhylogenyNode n23 = new PhylogenyNode( "n23/PIG_1", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n23.getNodeName().equals( "n23/PIG_1" ) ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n23 ).length() > 0 ) {
                return false;
            }
            if ( NHXParser.LIMIT_SPECIES_NAMES_TO_FIVE_CHARS ) {
                final PhylogenyNode a = new PhylogenyNode( "n10_ECOLI/1-2", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
                if ( !a.getNodeName().equals( "n10_ECOLI/1-2" ) ) {
                    return false;
                }
                if ( !PhylogenyMethods.getSpecies( a ).equals( "ECOLI" ) ) {
                    return false;
                }
                final PhylogenyNode b = new PhylogenyNode( "n10_ECOLI1/1-2", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
                if ( !b.getNodeName().equals( "n10_ECOLI1/1-2" ) ) {
                    return false;
                }
                if ( !PhylogenyMethods.getSpecies( b ).equals( "ECOLI" ) ) {
                    return false;
                }
                final PhylogenyNode c = new PhylogenyNode( "n10_RATAF12/1000-2000", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
                if ( !c.getNodeName().equals( "n10_RATAF12/1000-2000" ) ) {
                    return false;
                }
                if ( !PhylogenyMethods.getSpecies( c ).equals( "RATAF" ) ) {
                    return false;
                }
                final PhylogenyNode d = new PhylogenyNode( "n10_RAT1/1-2", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
                if ( !d.getNodeName().equals( "n10_RAT1/1-2" ) ) {
                    return false;
                }
                if ( !PhylogenyMethods.getSpecies( d ).equals( "RAT" ) ) {
                    return false;
                }
                final PhylogenyNode e = new PhylogenyNode( "n10_RAT1", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
                if ( !e.getNodeName().equals( "n10_RAT1" ) ) {
                    return false;
                }
                if ( !ForesterUtil.isEmpty( PhylogenyMethods.getSpecies( e ) ) ) {
                    return false;
                }
            }
            final PhylogenyNode n11 = new PhylogenyNode( "n111111_ECOLI/jdj:0.4", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n11.getNodeName().equals( "n111111_ECOLI/jdj" ) ) {
                return false;
            }
            if ( n11.getDistanceToParent() != 0.4 ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n11 ).equals( "ECOLI" ) ) {
                return false;
            }
            final PhylogenyNode n12 = new PhylogenyNode( "n111111-ECOLI---/jdj:0.4",
                                                         TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n12.getNodeName().equals( "n111111-ECOLI---/jdj" ) ) {
                return false;
            }
            if ( n12.getDistanceToParent() != 0.4 ) {
                return false;
            }
            if ( PhylogenyMethods.getSpecies( n12 ).length() > 0 ) {
                return false;
            }
            final Property tvu1 = n5.getNodeData().getProperties().getProperty( "tag1" );
            final Property tvu3 = n5.getNodeData().getProperties().getProperty( "tag3" );
            if ( !tvu1.getRef().equals( "tag1" ) ) {
                return false;
            }
            if ( !tvu1.getDataType().equals( "xsd:string" ) ) {
                return false;
            }
            if ( !tvu1.getUnit().equals( "unit1" ) ) {
                return false;
            }
            if ( !tvu1.getValue().equals( "value1" ) ) {
                return false;
            }
            if ( !tvu3.getRef().equals( "tag3" ) ) {
                return false;
            }
            if ( !tvu3.getDataType().equals( "xsd:string" ) ) {
                return false;
            }
            if ( !tvu3.getUnit().equals( "unit3" ) ) {
                return false;
            }
            if ( !tvu3.getValue().equals( "value3" ) ) {
                return false;
            }
            if ( n1.getNodeName().compareTo( "" ) != 0 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( n1 ) != Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                return false;
            }
            if ( n1.getDistanceToParent() != PhylogenyNode.DISTANCE_DEFAULT ) {
                return false;
            }
            if ( n2.getNodeName().compareTo( "" ) != 0 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( n2 ) != Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                return false;
            }
            if ( n2.getDistanceToParent() != PhylogenyNode.DISTANCE_DEFAULT ) {
                return false;
            }
            final PhylogenyNode n00 = new PhylogenyNode( "n7:0.000001[&&NHX:GN=gene_name:AC=accession123:ID=node_identifier:S=Ecoli:D=N:Co=N:B=100:T=1:On=100:SOn=100:SNn=100:W=2:C=0.0.0:XN=U=url_tag=www.yahoo.com]" );
            if ( !n00.getNodeData().getNodeIdentifier().getValue().equals( "node_identifier" ) ) {
                return false;
            }
            if ( !n00.getNodeData().getSequence().getName().equals( "gene_name" ) ) {
                return false;
            }
            if ( !n00.getNodeData().getSequence().getAccession().getValue().equals( "accession123" ) ) {
                return false;
            }
            if ( !n00.getNodeData().getProperties().getProperty( "url_tag" ).getRef().equals( "url_tag" ) ) {
                return false;
            }
            if ( n00.getNodeData().getProperties().getProperty( "url_tag" ).getAppliesTo() != Property.AppliesTo.NODE ) {
                return false;
            }
            if ( !n00.getNodeData().getProperties().getProperty( "url_tag" ).getDataType().equals( "xsd:anyURI" ) ) {
                return false;
            }
            if ( !n00.getNodeData().getProperties().getProperty( "url_tag" ).getValue().equals( "www.yahoo.com" ) ) {
                return false;
            }
            if ( !n00.getNodeData().getProperties().getProperty( "url_tag" ).getUnit().equals( "" ) ) {
                return false;
            }
            final PhylogenyNode nx = new PhylogenyNode( "n5:0.1[&&NHX:S=Ecoli:GN=gene_1]" );
            if ( !nx.getNodeData().getSequence().getName().equals( "gene_1" ) ) {
                return false;
            }
            final PhylogenyNode nx2 = new PhylogenyNode( "n5:0.1[&&NHX:S=Ecoli:G=gene_2]" );
            if ( !nx2.getNodeData().getSequence().getName().equals( "gene_2" ) ) {
                return false;
            }
            final PhylogenyNode n13 = new PhylogenyNode( "blah_12345/1-2", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n13.getNodeName().equals( "blah_12345/1-2" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n13 ).equals( "" ) ) {
                return false;
            }
            final PhylogenyNode n14 = new PhylogenyNode( "blah_12X45/1-2", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n14.getNodeName().equals( "blah_12X45/1-2" ) ) {
                return false;
            }
            if ( !PhylogenyMethods.getSpecies( n14 ).equals( "12X45" ) ) {
                return false;
            }
            final PhylogenyNode n15 = new PhylogenyNode( "something_wicked[123]", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n15.getNodeName().equals( "something_wicked" ) ) {
                return false;
            }
            if ( n15.getBranchData().getNumberOfConfidences() != 1 ) {
                return false;
            }
            if ( !isEqual( n15.getBranchData().getConfidence( 0 ).getValue(), 123 ) ) {
                return false;
            }
            final PhylogenyNode n16 = new PhylogenyNode( "something_wicked2[9]", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n16.getNodeName().equals( "something_wicked2" ) ) {
                return false;
            }
            if ( n16.getBranchData().getNumberOfConfidences() != 1 ) {
                return false;
            }
            if ( !isEqual( n16.getBranchData().getConfidence( 0 ).getValue(), 9 ) ) {
                return false;
            }
            final PhylogenyNode n17 = new PhylogenyNode( "something_wicked3[a]", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !n17.getNodeName().equals( "something_wicked3" ) ) {
                return false;
            }
            if ( n17.getBranchData().getNumberOfConfidences() != 0 ) {
                return false;
            }
            final PhylogenyNode n18 = new PhylogenyNode( ":0.5[91]", TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
            if ( !isEqual( n18.getDistanceToParent(), 0.5 ) ) {
                return false;
            }
            if ( n18.getBranchData().getNumberOfConfidences() != 1 ) {
                return false;
            }
            if ( !isEqual( n18.getBranchData().getConfidence( 0 ).getValue(), 91 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNHXParsing() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p1 = factory.create( "(A     [&&NHX:S=a_species],B1[&&NHX:S=b_species])", new NHXParser() )[ 0 ];
            if ( !p1.toNewHampshireX().equals( "(A[&&NHX:S=a_species],B1[&&NHX:S=b_species])" ) ) {
                return false;
            }
            final String p2_S = "(((((((A:0.2[&&NHX:S=qwerty]):0.2[&&NHX:S=uiop]):0.3[&&NHX:S=asdf]):0.4[&&NHX:S=zxc]):0.5[&&NHX:S=a]):0.6[&&NHX:S=asd]):0.7[&&NHX:S=za]):0.8[&&NHX:S=zaq]";
            final Phylogeny[] p2 = factory.create( p2_S, new NHXParser() );
            if ( !p2[ 0 ].toNewHampshireX().equals( p2_S ) ) {
                return false;
            }
            final String p2b_S = "(((((((A:0.2[&NHX:S=qwerty]):0.2[&:S=uiop]):0.3[&NHX:S=asdf]):0.4[S=zxc]):0.5[]):0.6[&&NH:S=asd]):0.7[&&HX:S=za]):0.8[&&:S=zaq]";
            final Phylogeny[] p2b = factory.create( p2b_S, new NHXParser() );
            if ( !p2b[ 0 ].toNewHampshireX().equals( "(((((((A:0.2):0.2):0.3):0.4):0.5):0.6):0.7):0.8" ) ) {
                return false;
            }
            final Phylogeny[] p3 = factory
                    .create( "[  comment&&NHX,())))](((((((A:0.2[&&NHX:S=qwerty]):0.2[&&NHX:S=uiop]):0.3[&&NHX:S=asdf]):0.4[&&NHX:S=zxc]):0.5[&&NHX:S=a]):0.6[&&NHX:S=asd]):0.7[&&NHX:S=za]):0.8[&&NHX:S=zaq]",
                             new NHXParser() );
            if ( !p3[ 0 ].toNewHampshireX().equals( p2_S ) ) {
                return false;
            }
            final Phylogeny[] p4 = factory
                    .create( "(((((((A:0.2[&&NHX:S=qwerty]):0.2[&&NHX:S=uiop]):0.3[&&NHX:S=asdf]):0.4[&&NHX:S=zxc]):0.5[&&NHX:S=a]):0.6[&&NHX:S=asd]):0.7[&&NHX:S=za]):0.8[&&NHX:S=zaq][comment(]",
                             new NHXParser() );
            if ( !p4[ 0 ].toNewHampshireX().equals( p2_S ) ) {
                return false;
            }
            final Phylogeny[] p5 = factory
                    .create( "[]  (  [][ ][   ]  ([((( &&NHXcomment only![[[[[[]([]((((A:0.2[&&NHX:S=q[comment )))]werty][,,,,))]):0.2[&&NHX:S=uiop]):0.3[&&NHX:S=a[comment,,))]sdf])[comment(((]:0.4[&&NHX:S=zxc][comment(((][comment(((]):0.5[&&NHX:S=a]):0.6[&&NHX:S=a[comment(((]sd]):0.7[&&NHX:S=za]):0.8[&&NHX:S=zaq][comment(((]",
                             new NHXParser() );
            if ( !p5[ 0 ].toNewHampshireX().equals( p2_S ) ) {
                return false;
            }
            final String p6_S_C = "(A[][][][1][22][333][4444][55555][666666][&&NHX:S=Aspecies],B[))],C,(AA,BB,CC,(CCC,DDD,EEE,[comment](FFFF,GGGG)x)y,D[comment]D,EE,FF,GG,HH),D,E,(EE,FF),F,G,H,(((((5)4)3)2)1),I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,(XX,(YY)),Y,Z)";
            final String p6_S_WO_C = "(A[&&NHX:S=Aspecies],B,C,(AA,BB,CC,(CCC,DDD,EEE,(FFFF,GGGG)x)y,DD,EE,FF,GG,HH),D,E,(EE,FF),F,G,H,(((((5)4)3)2)1),I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,(XX,(YY)),Y,Z)";
            final Phylogeny[] p6 = factory.create( p6_S_C, new NHXParser() );
            if ( !p6[ 0 ].toNewHampshireX().equals( p6_S_WO_C ) ) {
                return false;
            }
            final String p7_S_C = "(((A [&&NHX:S=species_a], B [&&NHX:S=Vstorri] , C   , D),(A,B,C,D[comment])[],[c][]([xxx]A[comment],[comment]B[comment][comment],[comment][comment]C[comment][comment],[comment][comment]D[comment][comment])[comment][comment],[comment]   [comment](A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C[comment][comment][comment][comment][comment]    [comment],D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),[comment][comment]((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)))";
            final String p7_S_WO_C = "(((A[&&NHX:S=species_a],B[&&NHX:S=Vstorri],C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)),((A,B,C,D),(A,B,C,D),(A,B,C,D),(A,B,C,D)))";
            final Phylogeny[] p7 = factory.create( p7_S_C, new NHXParser() );
            if ( !p7[ 0 ].toNewHampshireX().equals( p7_S_WO_C ) ) {
                return false;
            }
            final String p8_S_C = "[cmt](((([]([))))))](((((A[&&NHX:S= [a comment] a])))))))[too many comments!:)])),(((((((((B[&&NHX[ a comment in a bad place]:S   =b])))))[] []   )))),(((((((((C[&&NHX:S=c])   ))[,,, ])))))))";
            final String p8_S_WO_C = "((((((((((A[&&NHX:S=a]))))))))),(((((((((B[&&NHX:S=b]))))))))),(((((((((C[&&NHX:S=c]))))))))))";
            final Phylogeny[] p8 = factory.create( p8_S_C, new NHXParser() );
            if ( !p8[ 0 ].toNewHampshireX().equals( p8_S_WO_C ) ) {
                return false;
            }
            final Phylogeny p9 = factory.create( "((A:0.2,B:0.3):0.5[91],C:0.1)root:0.1[100]", new NHXParser() )[ 0 ];
            if ( !p9.toNewHampshireX().equals( "((A:0.2,B:0.3):0.5[&&NHX:B=91.0],C:0.1)root:0.1[&&NHX:B=100.0]" ) ) {
                return false;
            }
            final Phylogeny p10 = factory
                    .create( " [79]   ( (A [co mment] :0 .2[comment],B:0.3[com])[com ment]: 0. 5 \t[ 9 1 ][ comment],C: 0.1)[comment]root:0.1[100] [comment]",
                             new NHXParser() )[ 0 ];
            if ( !p10.toNewHampshireX().equals( "((A:0.2,B:0.3):0.5[&&NHX:B=91.0],C:0.1)root:0.1[&&NHX:B=100.0]" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNHXParsingQuotes() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final NHXParser p = new NHXParser();
            p.setIgnoreQuotes( false );
            final Phylogeny[] phylogenies_0 = factory.create( new File( Test.PATH_TO_TEST_DATA + "quotes.nhx" ), p );
            if ( phylogenies_0.length != 5 ) {
                return false;
            }
            final Phylogeny phy = phylogenies_0[ 4 ];
            if ( phy.getNumberOfExternalNodes() != 7 ) {
                return false;
            }
            if ( phy.getNodes( "a name in double quotes from tree ((a,b),c)" ).size() != 1 ) {
                return false;
            }
            if ( phy.getNodes( "charles darwin 'origin of species'" ).size() != 1 ) {
                return false;
            }
            if ( !phy.getNodes( "charles darwin 'origin of species'" ).get( 0 ).getNodeData().getTaxonomy()
                    .getScientificName().equals( "hsapiens" ) ) {
                return false;
            }
            if ( phy.getNodes( "shouldbetogether single quotes" ).size() != 1 ) {
                return false;
            }
            if ( phy.getNodes( "'single quotes' inside double quotes" ).size() != 1 ) {
                return false;
            }
            if ( phy.getNodes( "double quotes inside single quotes" ).size() != 1 ) {
                return false;
            }
            if ( phy.getNodes( "noquotes" ).size() != 1 ) {
                return false;
            }
            if ( phy.getNodes( "A   (  B    C '" ).size() != 1 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testPhylogenyBranch() {
        try {
            final PhylogenyNode a1 = new PhylogenyNode( "a" );
            final PhylogenyNode b1 = new PhylogenyNode( "b" );
            final PhylogenyBranch a1b1 = new PhylogenyBranch( a1, b1 );
            final PhylogenyBranch b1a1 = new PhylogenyBranch( b1, a1 );
            if ( !a1b1.equals( a1b1 ) ) {
                return false;
            }
            if ( !a1b1.equals( b1a1 ) ) {
                return false;
            }
            if ( !b1a1.equals( a1b1 ) ) {
                return false;
            }
            final PhylogenyBranch a1_b1 = new PhylogenyBranch( a1, b1, true );
            final PhylogenyBranch b1_a1 = new PhylogenyBranch( b1, a1, true );
            final PhylogenyBranch a1_b1_ = new PhylogenyBranch( a1, b1, false );
            if ( a1_b1.equals( b1_a1 ) ) {
                return false;
            }
            if ( a1_b1.equals( a1_b1_ ) ) {
                return false;
            }
            final PhylogenyBranch b1_a1_ = new PhylogenyBranch( b1, a1, false );
            if ( !a1_b1.equals( b1_a1_ ) ) {
                return false;
            }
            if ( a1_b1_.equals( b1_a1_ ) ) {
                return false;
            }
            if ( !a1_b1_.equals( b1_a1 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testPostOrderIterator() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0 = factory.create( "((A,B)ab,(C,D)cd)r", new NHXParser() )[ 0 ];
            PhylogenyNodeIterator it0;
            for( it0 = t0.iteratorPostorder(); it0.hasNext(); ) {
                it0.next();
            }
            for( it0.reset(); it0.hasNext(); ) {
                it0.next();
            }
            final Phylogeny t1 = factory.create( "(((A,B)ab,(C,D)cd)abcd,((E,F)ef,(G,H)gh)efgh)r", new NHXParser() )[ 0 ];
            final PhylogenyNodeIterator it = t1.iteratorPostorder();
            if ( !it.next().getNodeName().equals( "A" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "B" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "ab" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "C" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "D" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "cd" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "abcd" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "E" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "F" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "ef" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "G" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "H" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "gh" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "efgh" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "r" ) ) {
                return false;
            }
            if ( it.hasNext() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testPreOrderIterator() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0 = factory.create( "((A,B)ab,(C,D)cd)r", new NHXParser() )[ 0 ];
            PhylogenyNodeIterator it0;
            for( it0 = t0.iteratorPreorder(); it0.hasNext(); ) {
                it0.next();
            }
            for( it0.reset(); it0.hasNext(); ) {
                it0.next();
            }
            PhylogenyNodeIterator it = t0.iteratorPreorder();
            if ( !it.next().getNodeName().equals( "r" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "ab" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "A" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "B" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "cd" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "C" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "D" ) ) {
                return false;
            }
            if ( it.hasNext() ) {
                return false;
            }
            final Phylogeny t1 = factory.create( "(((A,B)ab,(C,D)cd)abcd,((E,F)ef,(G,H)gh)efgh)r", new NHXParser() )[ 0 ];
            it = t1.iteratorPreorder();
            if ( !it.next().getNodeName().equals( "r" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "abcd" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "ab" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "A" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "B" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "cd" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "C" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "D" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "efgh" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "ef" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "E" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "F" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "gh" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "G" ) ) {
                return false;
            }
            if ( !it.next().getNodeName().equals( "H" ) ) {
                return false;
            }
            if ( it.hasNext() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testReIdMethods() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p = factory.create( "((1,2)A,(((X,Y,Z)a,b)3)B,(4,5,6)C)r", new NHXParser() )[ 0 ];
            final int last = p.levelOrderReID( 120 );
            if ( last != 124 ) {
                return false;
            }
            if ( p.getNode( "r" ).getNodeId() != 120 ) {
                return false;
            }
            if ( p.getNode( "A" ).getNodeId() != 121 ) {
                return false;
            }
            if ( p.getNode( "B" ).getNodeId() != 121 ) {
                return false;
            }
            if ( p.getNode( "C" ).getNodeId() != 121 ) {
                return false;
            }
            if ( p.getNode( "1" ).getNodeId() != 122 ) {
                return false;
            }
            if ( p.getNode( "2" ).getNodeId() != 122 ) {
                return false;
            }
            if ( p.getNode( "3" ).getNodeId() != 122 ) {
                return false;
            }
            if ( p.getNode( "4" ).getNodeId() != 122 ) {
                return false;
            }
            if ( p.getNode( "5" ).getNodeId() != 122 ) {
                return false;
            }
            if ( p.getNode( "6" ).getNodeId() != 122 ) {
                return false;
            }
            if ( p.getNode( "a" ).getNodeId() != 123 ) {
                return false;
            }
            if ( p.getNode( "b" ).getNodeId() != 123 ) {
                return false;
            }
            if ( p.getNode( "X" ).getNodeId() != 124 ) {
                return false;
            }
            if ( p.getNode( "Y" ).getNodeId() != 124 ) {
                return false;
            }
            if ( p.getNode( "Z" ).getNodeId() != 124 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testRerooting() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t1 = factory.create( "((A:1,B:2)AB:1[&&NHX:B=55],(C:3,D:5)CD:3[&&NHX:B=10])ABCD:0.5",
                                                 new NHXParser() )[ 0 ];
            if ( !t1.isRooted() ) {
                return false;
            }
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "CD" ) );
            t1.reRoot( t1.getNode( "A" ) );
            t1.reRoot( t1.getNode( "B" ) );
            t1.reRoot( t1.getNode( "AB" ) );
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "C" ) );
            t1.reRoot( t1.getNode( "CD" ) );
            t1.reRoot( t1.getNode( "A" ) );
            t1.reRoot( t1.getNode( "B" ) );
            t1.reRoot( t1.getNode( "AB" ) );
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "C" ) );
            t1.reRoot( t1.getNode( "A" ) );
            t1.reRoot( t1.getNode( "B" ) );
            t1.reRoot( t1.getNode( "AB" ) );
            t1.reRoot( t1.getNode( "C" ) );
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "CD" ) );
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "A" ) );
            t1.reRoot( t1.getNode( "B" ) );
            t1.reRoot( t1.getNode( "AB" ) );
            t1.reRoot( t1.getNode( "C" ) );
            t1.reRoot( t1.getNode( "D" ) );
            t1.reRoot( t1.getNode( "CD" ) );
            t1.reRoot( t1.getNode( "D" ) );
            if ( !isEqual( t1.getNode( "A" ).getDistanceToParent(), 1 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "B" ).getDistanceToParent(), 2 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "C" ).getDistanceToParent(), 3 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "D" ).getDistanceToParent(), 2.5 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "CD" ).getDistanceToParent(), 2.5 ) ) {
                return false;
            }
            if ( !isEqual( t1.getNode( "AB" ).getDistanceToParent(), 4 ) ) {
                return false;
            }
            final Phylogeny t2 = factory.create( "(((A:1,B:2)AB:10[&&NHX:B=55],C)ABC:3[&&NHX:B=33],D:5)ABCD:0.5",
                                                 new NHXParser() )[ 0 ];
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "ABC" ) );
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "C" ) );
            t2.reRoot( t2.getNode( "ABC" ) );
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "C" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "ABC" ) );
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "C" ) );
            t2.reRoot( t2.getNode( "ABC" ) );
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "C" ) );
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "C" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "ABC" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "A" ) );
            t2.reRoot( t2.getNode( "B" ) );
            t2.reRoot( t2.getNode( "AB" ) );
            t2.reRoot( t2.getNode( "C" ) );
            t2.reRoot( t2.getNode( "D" ) );
            t2.reRoot( t2.getNode( "ABC" ) );
            t2.reRoot( t2.getNode( "D" ) );
            if ( !isEqual( t2.getNode( "AB" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "ABC" ).getBranchData().getConfidence( 0 ).getValue(), 33 ) ) {
                return false;
            }
            t2.reRoot( t2.getNode( "ABC" ) );
            if ( !isEqual( t2.getNode( "AB" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "ABC" ).getBranchData().getConfidence( 0 ).getValue(), 33 ) ) {
                return false;
            }
            t2.reRoot( t2.getNode( "AB" ) );
            if ( !isEqual( t2.getNode( "AB" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "ABC" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "D" ).getBranchData().getConfidence( 0 ).getValue(), 33 ) ) {
                return false;
            }
            t2.reRoot( t2.getNode( "AB" ) );
            if ( !isEqual( t2.getNode( "AB" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "ABC" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "D" ).getBranchData().getConfidence( 0 ).getValue(), 33 ) ) {
                return false;
            }
            t2.reRoot( t2.getNode( "D" ) );
            if ( !isEqual( t2.getNode( "AB" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "ABC" ).getBranchData().getConfidence( 0 ).getValue(), 33 ) ) {
                return false;
            }
            t2.reRoot( t2.getNode( "ABC" ) );
            if ( !isEqual( t2.getNode( "AB" ).getBranchData().getConfidence( 0 ).getValue(), 55 ) ) {
                return false;
            }
            if ( !isEqual( t2.getNode( "ABC" ).getBranchData().getConfidence( 0 ).getValue(), 33 ) ) {
                return false;
            }
            final Phylogeny t3 = factory.create( "(A[&&NHX:B=10],B[&&NHX:B=20],C[&&NHX:B=30],D[&&NHX:B=40])",
                                                 new NHXParser() )[ 0 ];
            t3.reRoot( t3.getNode( "B" ) );
            if ( t3.getNode( "B" ).getBranchData().getConfidence( 0 ).getValue() != 20 ) {
                return false;
            }
            if ( t3.getNode( "A" ).getParent().getBranchData().getConfidence( 0 ).getValue() != 20 ) {
                return false;
            }
            if ( t3.getNode( "A" ).getParent().getNumberOfDescendants() != 3 ) {
                return false;
            }
            t3.reRoot( t3.getNode( "B" ) );
            if ( t3.getNode( "B" ).getBranchData().getConfidence( 0 ).getValue() != 20 ) {
                return false;
            }
            if ( t3.getNode( "A" ).getParent().getBranchData().getConfidence( 0 ).getValue() != 20 ) {
                return false;
            }
            if ( t3.getNode( "A" ).getParent().getNumberOfDescendants() != 3 ) {
                return false;
            }
            t3.reRoot( t3.getRoot() );
            if ( t3.getNode( "B" ).getBranchData().getConfidence( 0 ).getValue() != 20 ) {
                return false;
            }
            if ( t3.getNode( "A" ).getParent().getBranchData().getConfidence( 0 ).getValue() != 20 ) {
                return false;
            }
            if ( t3.getNode( "A" ).getParent().getNumberOfDescendants() != 3 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSDIse() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny species1 = factory.create( "[&&NHX:S=yeast]", new NHXParser() )[ 0 ];
            final Phylogeny gene1 = factory.create( "(A1[&&NHX:S=yeast],A2[&&NHX:S=yeast])", new NHXParser() )[ 0 ];
            gene1.setRooted( true );
            species1.setRooted( true );
            final SDI sdi = new SDIse( gene1, species1 );
            if ( !gene1.getRoot().isDuplication() ) {
                return false;
            }
            final Phylogeny species2 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=B]),[&&NHX:S=C]),[&&NHX:S=D]),([&&NHX:S=E],[&&NHX:S=F]))",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene2 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=B])ab,[&&NHX:S=C])abc,[&&NHX:S=D])abcd,([&&NHX:S=E],[&&NHX:S=F])ef)r",
                             new NHXParser() )[ 0 ];
            species2.setRooted( true );
            gene2.setRooted( true );
            final SDI sdi2 = new SDIse( gene2, species2 );
            if ( sdi2.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( !gene2.getNode( "ab" ).isSpeciation() ) {
                return false;
            }
            if ( !gene2.getNode( "ab" ).isHasAssignedEvent() ) {
                return false;
            }
            if ( !gene2.getNode( "abc" ).isSpeciation() ) {
                return false;
            }
            if ( !gene2.getNode( "abc" ).isHasAssignedEvent() ) {
                return false;
            }
            if ( !gene2.getNode( "r" ).isSpeciation() ) {
                return false;
            }
            if ( !gene2.getNode( "r" ).isHasAssignedEvent() ) {
                return false;
            }
            final Phylogeny species3 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=B]),[&&NHX:S=C]),[&&NHX:S=D]),([&&NHX:S=E],[&&NHX:S=F]))",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene3 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=A])aa,[&&NHX:S=C])abc,[&&NHX:S=D])abcd,([&&NHX:S=E],[&&NHX:S=F])ef)r",
                             new NHXParser() )[ 0 ];
            species3.setRooted( true );
            gene3.setRooted( true );
            final SDI sdi3 = new SDIse( gene3, species3 );
            if ( sdi3.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( !gene3.getNode( "aa" ).isDuplication() ) {
                return false;
            }
            if ( !gene3.getNode( "aa" ).isHasAssignedEvent() ) {
                return false;
            }
            final Phylogeny species4 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=B]),[&&NHX:S=C]),[&&NHX:S=D]),([&&NHX:S=E],[&&NHX:S=F]))",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene4 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=C])ac,[&&NHX:S=B])abc,[&&NHX:S=D])abcd,([&&NHX:S=E],[&&NHX:S=F])ef)r",
                             new NHXParser() )[ 0 ];
            species4.setRooted( true );
            gene4.setRooted( true );
            final SDI sdi4 = new SDIse( gene4, species4 );
            if ( sdi4.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( !gene4.getNode( "ac" ).isSpeciation() ) {
                return false;
            }
            if ( !gene4.getNode( "abc" ).isDuplication() ) {
                return false;
            }
            if ( gene4.getNode( "abcd" ).isDuplication() ) {
                return false;
            }
            if ( species4.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            if ( gene4.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            final Phylogeny species5 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=B]),[&&NHX:S=C]),[&&NHX:S=D]),([&&NHX:S=E],[&&NHX:S=F]))",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene5 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=D])ad,[&&NHX:S=C])adc,[&&NHX:S=B])abcd,([&&NHX:S=E],[&&NHX:S=F])ef)r",
                             new NHXParser() )[ 0 ];
            species5.setRooted( true );
            gene5.setRooted( true );
            final SDI sdi5 = new SDIse( gene5, species5 );
            if ( sdi5.getDuplicationsSum() != 2 ) {
                return false;
            }
            if ( !gene5.getNode( "ad" ).isSpeciation() ) {
                return false;
            }
            if ( !gene5.getNode( "adc" ).isDuplication() ) {
                return false;
            }
            if ( !gene5.getNode( "abcd" ).isDuplication() ) {
                return false;
            }
            if ( species5.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            if ( gene5.getNumberOfExternalNodes() != 6 ) {
                return false;
            }
            // Trees from Louxin Zhang 1997 "On a Mirkin-Muchnik-Smith
            // Conjecture for Comparing Molecular Phylogenies"
            // J. of Comput Bio. Vol. 4, No 2, pp.177-187
            final Phylogeny species6 = factory
                    .create( "(((1:[&&NHX:S=1],5:[&&NHX:S=5])1-5,((4:[&&NHX:S=4],6:[&&NHX:S=6])4-6,2:[&&NHX:S=2])4-6-2)1-5-4-6-2,"
                                     + "((9:[&&NHX:S=9],3:[&&NHX:S=3])9-3,(8:[&&NHX:S=8],7:[&&NHX:S=7])8-7)9-3-8-7)",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene6 = factory
                    .create( "(((1:0.1[&&NHX:S=1],2:0.1[&&NHX:S=2])1-2:0.1,3:0.1[&&NHX:S=3])1-2-3:0.1,"
                                     + "((4:0.1[&&NHX:S=4],(5:0.1[&&NHX:S=5],6:0.1[&&NHX:S=6])5-6:0.1)4-5-6:0.1,"
                                     + "(7:0.1[&&NHX:S=7],(8:0.1[&&NHX:S=8],9:0.1[&&NHX:S=9])8-9:0.1)7-8-9:0.1)4-5-6-7-8-9:0.1)r;",
                             new NHXParser() )[ 0 ];
            species6.setRooted( true );
            gene6.setRooted( true );
            final SDI sdi6 = new SDIse( gene6, species6 );
            if ( sdi6.getDuplicationsSum() != 3 ) {
                return false;
            }
            if ( !gene6.getNode( "r" ).isDuplication() ) {
                return false;
            }
            if ( !gene6.getNode( "4-5-6" ).isDuplication() ) {
                return false;
            }
            if ( !gene6.getNode( "7-8-9" ).isDuplication() ) {
                return false;
            }
            if ( !gene6.getNode( "1-2" ).isSpeciation() ) {
                return false;
            }
            if ( !gene6.getNode( "1-2-3" ).isSpeciation() ) {
                return false;
            }
            if ( !gene6.getNode( "5-6" ).isSpeciation() ) {
                return false;
            }
            if ( !gene6.getNode( "8-9" ).isSpeciation() ) {
                return false;
            }
            if ( !gene6.getNode( "4-5-6-7-8-9" ).isSpeciation() ) {
                return false;
            }
            sdi6.computeMappingCostL();
            if ( sdi6.computeMappingCostL() != 17 ) {
                return false;
            }
            if ( species6.getNumberOfExternalNodes() != 9 ) {
                return false;
            }
            if ( gene6.getNumberOfExternalNodes() != 9 ) {
                return false;
            }
            final Phylogeny species7 = Test.createPhylogeny( "(((((((" + "([&&NHX:S=a1],[&&NHX:S=a2]),"
                    + "([&&NHX:S=b1],[&&NHX:S=b2])" + "),[&&NHX:S=x]),(" + "([&&NHX:S=m1],[&&NHX:S=m2]),"
                    + "([&&NHX:S=n1],[&&NHX:S=n2])" + ")),(" + "([&&NHX:S=i1],[&&NHX:S=i2]),"
                    + "([&&NHX:S=j1],[&&NHX:S=j2])" + ")),(" + "([&&NHX:S=e1],[&&NHX:S=e2]),"
                    + "([&&NHX:S=f1],[&&NHX:S=f2])" + ")),[&&NHX:S=y]),[&&NHX:S=z])" );
            species7.setRooted( true );
            final Phylogeny gene7_1 = Test
                    .createPhylogeny( "((((((((a1[&&NHX:S=a1],a2[&&NHX:S=a2]),b1[&&NHX:S=b1]),x[&&NHX:S=x]),m1[&&NHX:S=m1]),i1[&&NHX:S=i1]),e1[&&NHX:S=e1]),y[&&NHX:S=y]),z[&&NHX:S=z])" );
            gene7_1.setRooted( true );
            final SDI sdi7 = new SDIse( gene7_1, species7 );
            if ( sdi7.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "a2" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "x" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "m1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "i1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "e1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "y" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_1, "a1", "z" ).isSpeciation() ) {
                return false;
            }
            final Phylogeny gene7_2 = Test
                    .createPhylogeny( "(((((((((a1[&&NHX:S=a1],a2[&&NHX:S=a2]),b1[&&NHX:S=b1]),x[&&NHX:S=x]),m1[&&NHX:S=m1]),i1[&&NHX:S=i1]),j2[&&NHX:S=j2]),e1[&&NHX:S=e1]),y[&&NHX:S=y]),z[&&NHX:S=z])" );
            gene7_2.setRooted( true );
            final SDI sdi7_2 = new SDIse( gene7_2, species7 );
            if ( sdi7_2.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "a2" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "x" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "m1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "i1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "j2" ).isDuplication() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "e1" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "y" ).isSpeciation() ) {
                return false;
            }
            if ( !Test.getEvent( gene7_2, "a1", "z" ).isSpeciation() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            return false;
        }
        return true;
    }

    private static boolean testSDIunrooted() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p0 = factory.create( "((((A,B)ab,(C1,C2)cc)abc,D)abcd,(E,F)ef)abcdef", new NHXParser() )[ 0 ];
            final List<PhylogenyBranch> l = SDIR.getBranchesInPreorder( p0 );
            final Iterator<PhylogenyBranch> iter = l.iterator();
            PhylogenyBranch br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "abcd" ) && !br.getFirstNode().getNodeName().equals( "ef" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "abcd" ) && !br.getSecondNode().getNodeName().equals( "ef" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "abcd" ) && !br.getFirstNode().getNodeName().equals( "abc" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "abcd" ) && !br.getSecondNode().getNodeName().equals( "abc" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "abc" ) && !br.getFirstNode().getNodeName().equals( "ab" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "abc" ) && !br.getSecondNode().getNodeName().equals( "ab" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "ab" ) && !br.getFirstNode().getNodeName().equals( "A" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "ab" ) && !br.getSecondNode().getNodeName().equals( "A" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "ab" ) && !br.getFirstNode().getNodeName().equals( "B" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "ab" ) && !br.getSecondNode().getNodeName().equals( "B" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "ab" ) && !br.getFirstNode().getNodeName().equals( "abc" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "ab" ) && !br.getSecondNode().getNodeName().equals( "abc" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "abc" ) && !br.getFirstNode().getNodeName().equals( "cc" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "abc" ) && !br.getSecondNode().getNodeName().equals( "cc" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "C1" ) && !br.getFirstNode().getNodeName().equals( "cc" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "C1" ) && !br.getSecondNode().getNodeName().equals( "cc" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "C2" ) && !br.getFirstNode().getNodeName().equals( "cc" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "C2" ) && !br.getSecondNode().getNodeName().equals( "cc" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "abc" ) && !br.getFirstNode().getNodeName().equals( "cc" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "abc" ) && !br.getSecondNode().getNodeName().equals( "cc" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "abc" ) && !br.getFirstNode().getNodeName().equals( "abcd" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "abc" ) && !br.getSecondNode().getNodeName().equals( "abcd" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "abcd" ) && !br.getFirstNode().getNodeName().equals( "D" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "abcd" ) && !br.getSecondNode().getNodeName().equals( "D" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "ef" ) && !br.getFirstNode().getNodeName().equals( "abcd" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "ef" ) && !br.getSecondNode().getNodeName().equals( "abcd" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "ef" ) && !br.getFirstNode().getNodeName().equals( "E" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "ef" ) && !br.getSecondNode().getNodeName().equals( "E" ) ) {
                return false;
            }
            br = iter.next();
            if ( !br.getFirstNode().getNodeName().equals( "ef" ) && !br.getFirstNode().getNodeName().equals( "F" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "ef" ) && !br.getSecondNode().getNodeName().equals( "F" ) ) {
                return false;
            }
            if ( iter.hasNext() ) {
                return false;
            }
            final Phylogeny p1 = factory.create( "(C,(A,B)ab)abc", new NHXParser() )[ 0 ];
            final List<PhylogenyBranch> l1 = SDIR.getBranchesInPreorder( p1 );
            final Iterator<PhylogenyBranch> iter1 = l1.iterator();
            br = iter1.next();
            if ( !br.getFirstNode().getNodeName().equals( "ab" ) && !br.getFirstNode().getNodeName().equals( "C" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "ab" ) && !br.getSecondNode().getNodeName().equals( "C" ) ) {
                return false;
            }
            br = iter1.next();
            if ( !br.getFirstNode().getNodeName().equals( "ab" ) && !br.getFirstNode().getNodeName().equals( "A" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "ab" ) && !br.getSecondNode().getNodeName().equals( "A" ) ) {
                return false;
            }
            br = iter1.next();
            if ( !br.getFirstNode().getNodeName().equals( "ab" ) && !br.getFirstNode().getNodeName().equals( "B" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "ab" ) && !br.getSecondNode().getNodeName().equals( "B" ) ) {
                return false;
            }
            if ( iter1.hasNext() ) {
                return false;
            }
            final Phylogeny p2 = factory.create( "((A,B)ab,C)abc", new NHXParser() )[ 0 ];
            final List<PhylogenyBranch> l2 = SDIR.getBranchesInPreorder( p2 );
            final Iterator<PhylogenyBranch> iter2 = l2.iterator();
            br = iter2.next();
            if ( !br.getFirstNode().getNodeName().equals( "ab" ) && !br.getFirstNode().getNodeName().equals( "C" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "ab" ) && !br.getSecondNode().getNodeName().equals( "C" ) ) {
                return false;
            }
            br = iter2.next();
            if ( !br.getFirstNode().getNodeName().equals( "ab" ) && !br.getFirstNode().getNodeName().equals( "A" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "ab" ) && !br.getSecondNode().getNodeName().equals( "A" ) ) {
                return false;
            }
            br = iter2.next();
            if ( !br.getFirstNode().getNodeName().equals( "ab" ) && !br.getFirstNode().getNodeName().equals( "B" ) ) {
                return false;
            }
            if ( !br.getSecondNode().getNodeName().equals( "ab" ) && !br.getSecondNode().getNodeName().equals( "B" ) ) {
                return false;
            }
            if ( iter2.hasNext() ) {
                return false;
            }
            final Phylogeny species0 = factory
                    .create( "(((([&&NHX:S=A],[&&NHX:S=B]),[&&NHX:S=C]),[&&NHX:S=D]),([&&NHX:S=E],[&&NHX:S=F]))",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene1 = factory
                    .create( "(((((A:0.6[&&NHX:S=A],B:0.1[&&NHX:S=B])ab:0.1,C:0.1[&&NHX:S=C])abc:0.3,D:1.0[&&NHX:S=D])abcd:0.2,E:0.1[&&NHX:S=E])abcde:0.2,F:0.2[&&NHX:S=F])",
                             new NHXParser() )[ 0 ];
            species0.setRooted( true );
            gene1.setRooted( true );
            final SDIR sdi_unrooted = new SDIR();
            sdi_unrooted.infer( gene1, species0, false, true, true, true, 10 );
            if ( sdi_unrooted.getCount() != 1 ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalDuplications() != 0 ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalDiffInSubTreeHeights(), 0.4 ) ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalTreeHeight(), 1.0 ) ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalMappingCost() != Integer.MAX_VALUE ) {
                return false;
            }
            final Phylogeny gene2 = factory
                    .create( "(((((A:2.6[&&NHX:S=A],B:0.1[&&NHX:S=B])ab:0.1,C:0.1[&&NHX:S=C])abc:0.3,D:1.0[&&NHX:S=D])abcd:0.2,E:0.1[&&NHX:S=E])abcde:0.2,F:0.2[&&NHX:S=F])",
                             new NHXParser() )[ 0 ];
            gene2.setRooted( true );
            sdi_unrooted.infer( gene2, species0, false, false, true, true, 10 );
            if ( sdi_unrooted.getCount() != 1 ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalDuplications() != 3 ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalDiffInSubTreeHeights(), 0.0 ) ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalTreeHeight(), 2.0 ) ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalMappingCost() != Integer.MAX_VALUE ) {
                return false;
            }
            final Phylogeny species6 = factory
                    .create( "(((1:[&&NHX:S=1],5:[&&NHX:S=5])1-5,((4:[&&NHX:S=4],6:[&&NHX:S=6])4-6,2:[&&NHX:S=2])4-6-2)1-5-4-6-2,"
                                     + "((9:[&&NHX:S=9],3:[&&NHX:S=3])9-3,(8:[&&NHX:S=8],7:[&&NHX:S=7])8-7)9-3-8-7)",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene6 = factory
                    .create( "((5:0.1[&&NHX:S=5],6:0.1[&&NHX:S=6])5-6:0.05[&&NHX:S=6],(4:0.1[&&NHX:S=4],"
                                     + "(((1:0.1[&&NHX:S=1],2:0.1[&&NHX:S=2])1-2:0.1[&&NHX:S=2],3:0.25[&&NHX:S=3])1-2-3:0.2[&&NHX:S=2],"
                                     + "(7:0.1[&&NHX:S=7],(8:0.1[&&NHX:S=8],"
                                     + "9:0.1[&&NHX:S=9])8-9:0.1[&&NHX:S=9])7-8-9:0.1[&&NHX:S=8])"
                                     + "4-5-6-7-8-9:0.1[&&NHX:S=5])4-5-6:0.05[&&NHX:S=5])",
                             new NHXParser() )[ 0 ];
            species6.setRooted( true );
            gene6.setRooted( true );
            Phylogeny[] p6 = sdi_unrooted.infer( gene6, species6, false, true, true, true, 10 );
            if ( sdi_unrooted.getCount() != 1 ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalDiffInSubTreeHeights(), 0.0 ) ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalTreeHeight(), 0.375 ) ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalDuplications() != 3 ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalMappingCost() != Integer.MAX_VALUE ) {
                return false;
            }
            if ( !p6[ 0 ].getRoot().isDuplication() ) {
                return false;
            }
            if ( !p6[ 0 ].getNode( "4-5-6" ).isDuplication() ) {
                return false;
            }
            if ( !p6[ 0 ].getNode( "7-8-9" ).isDuplication() ) {
                return false;
            }
            if ( p6[ 0 ].getNode( "1-2" ).isDuplication() ) {
                return false;
            }
            if ( p6[ 0 ].getNode( "1-2-3" ).isDuplication() ) {
                return false;
            }
            if ( p6[ 0 ].getNode( "5-6" ).isDuplication() ) {
                return false;
            }
            if ( p6[ 0 ].getNode( "8-9" ).isDuplication() ) {
                return false;
            }
            if ( p6[ 0 ].getNode( "4-5-6-7-8-9" ).isDuplication() ) {
                return false;
            }
            p6 = null;
            final Phylogeny species7 = factory
                    .create( "(((1:[&&NHX:S=1],5:[&&NHX:S=5])1-5,((4:[&&NHX:S=4],6:[&&NHX:S=6])4-6,2:[&&NHX:S=2])4-6-2)1-5-4-6-2,"
                                     + "((9:[&&NHX:S=9],3:[&&NHX:S=3])9-3,(8:[&&NHX:S=8],7:[&&NHX:S=7])8-7)9-3-8-7)",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene7 = factory
                    .create( "((5:0.1[&&NHX:S=5],6:0.1[&&NHX:S=6])5-6:0.05[&&NHX:S=6],(4:0.1[&&NHX:S=4],"
                                     + "(((1:0.1[&&NHX:S=1],2:0.1[&&NHX:S=2])1-2:0.1[&&NHX:S=2],3:0.25[&&NHX:S=3])1-2-3:0.2[&&NHX:S=2],"
                                     + "(7:0.1[&&NHX:S=7],(8:0.1[&&NHX:S=8],"
                                     + "9:0.1[&&NHX:S=9])8-9:0.1[&&NHX:S=9])7-8-9:0.1[&&NHX:S=8])"
                                     + "4-5-6-7-8-9:0.1[&&NHX:S=5])4-5-6:0.05[&&NHX:S=5])",
                             new NHXParser() )[ 0 ];
            species7.setRooted( true );
            gene7.setRooted( true );
            Phylogeny[] p7 = sdi_unrooted.infer( gene7, species7, true, true, true, true, 10 );
            if ( sdi_unrooted.getCount() != 1 ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalDiffInSubTreeHeights(), 0.0 ) ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalTreeHeight(), 0.375 ) ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalDuplications() != 3 ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalMappingCost() != 17 ) {
                return false;
            }
            if ( !p7[ 0 ].getRoot().isDuplication() ) {
                return false;
            }
            if ( !p7[ 0 ].getNode( "4-5-6" ).isDuplication() ) {
                return false;
            }
            if ( !p7[ 0 ].getNode( "7-8-9" ).isDuplication() ) {
                return false;
            }
            if ( p7[ 0 ].getNode( "1-2" ).isDuplication() ) {
                return false;
            }
            if ( p7[ 0 ].getNode( "1-2-3" ).isDuplication() ) {
                return false;
            }
            if ( p7[ 0 ].getNode( "5-6" ).isDuplication() ) {
                return false;
            }
            if ( p7[ 0 ].getNode( "8-9" ).isDuplication() ) {
                return false;
            }
            if ( p7[ 0 ].getNode( "4-5-6-7-8-9" ).isDuplication() ) {
                return false;
            }
            p7 = null;
            final Phylogeny species8 = factory
                    .create( "(((1:[&&NHX:S=1],5:[&&NHX:S=5])1-5,((4:[&&NHX:S=4],6:[&&NHX:S=6])4-6,2:[&&NHX:S=2])4-6-2)1-5-4-6-2,"
                                     + "((9:[&&NHX:S=9],3:[&&NHX:S=3])9-3,(8:[&&NHX:S=8],7:[&&NHX:S=7])8-7)9-3-8-7)",
                             new NHXParser() )[ 0 ];
            final Phylogeny gene8 = factory
                    .create( "((5:0.1[&&NHX:S=5],6:0.1[&&NHX:S=6])5-6:0.05[&&NHX:S=6],(4:0.1[&&NHX:S=4],"
                                     + "(((1:0.1[&&NHX:S=1],2:0.1[&&NHX:S=2])1-2:0.1[&&NHX:S=2],3:0.25[&&NHX:S=3])1-2-3:0.2[&&NHX:S=2],"
                                     + "(7:0.1[&&NHX:S=7],(8:0.1[&&NHX:S=8],"
                                     + "9:0.1[&&NHX:S=9])8-9:0.1[&&NHX:S=9])7-8-9:0.1[&&NHX:S=8])"
                                     + "4-5-6-7-8-9:0.1[&&NHX:S=5])4-5-6:0.05[&&NHX:S=5])",
                             new NHXParser() )[ 0 ];
            species8.setRooted( true );
            gene8.setRooted( true );
            Phylogeny[] p8 = sdi_unrooted.infer( gene8, species8, false, false, true, true, 10 );
            if ( sdi_unrooted.getCount() != 1 ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalDiffInSubTreeHeights(), 0.0 ) ) {
                return false;
            }
            if ( !Test.isEqual( sdi_unrooted.getMinimalTreeHeight(), 0.375 ) ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalDuplications() != 3 ) {
                return false;
            }
            if ( sdi_unrooted.getMinimalMappingCost() != Integer.MAX_VALUE ) {
                return false;
            }
            if ( !p8[ 0 ].getRoot().isDuplication() ) {
                return false;
            }
            if ( !p8[ 0 ].getNode( "4-5-6" ).isDuplication() ) {
                return false;
            }
            if ( !p8[ 0 ].getNode( "7-8-9" ).isDuplication() ) {
                return false;
            }
            if ( p8[ 0 ].getNode( "1-2" ).isDuplication() ) {
                return false;
            }
            if ( p8[ 0 ].getNode( "1-2-3" ).isDuplication() ) {
                return false;
            }
            if ( p8[ 0 ].getNode( "5-6" ).isDuplication() ) {
                return false;
            }
            if ( p8[ 0 ].getNode( "8-9" ).isDuplication() ) {
                return false;
            }
            if ( p8[ 0 ].getNode( "4-5-6-7-8-9" ).isDuplication() ) {
                return false;
            }
            p8 = null;
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSubtreeDeletion() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t1 = factory.create( "((A,B,C)abc,(D,E,F)def)r", new NHXParser() )[ 0 ];
            t1.deleteSubtree( t1.getNode( "A" ), false );
            if ( t1.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            t1.toNewHampshireX();
            t1.deleteSubtree( t1.getNode( "E" ), false );
            if ( t1.getNumberOfExternalNodes() != 4 ) {
                return false;
            }
            t1.toNewHampshireX();
            t1.deleteSubtree( t1.getNode( "F" ), false );
            if ( t1.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            t1.toNewHampshireX();
            t1.deleteSubtree( t1.getNode( "D" ), false );
            t1.toNewHampshireX();
            if ( t1.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "def" ), false );
            t1.toNewHampshireX();
            if ( t1.getNumberOfExternalNodes() != 2 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "B" ), false );
            t1.toNewHampshireX();
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "C" ), false );
            t1.toNewHampshireX();
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "abc" ), false );
            t1.toNewHampshireX();
            if ( t1.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
            t1.deleteSubtree( t1.getNode( "r" ), false );
            if ( t1.getNumberOfExternalNodes() != 0 ) {
                return false;
            }
            if ( !t1.isEmpty() ) {
                return false;
            }
            final Phylogeny t2 = factory.create( "(((1,2,3)A,B,C)abc,(D,E,F)def)r", new NHXParser() )[ 0 ];
            t2.deleteSubtree( t2.getNode( "A" ), false );
            t2.toNewHampshireX();
            if ( t2.getNumberOfExternalNodes() != 5 ) {
                return false;
            }
            t2.deleteSubtree( t2.getNode( "abc" ), false );
            t2.toNewHampshireX();
            if ( t2.getNumberOfExternalNodes() != 3 ) {
                return false;
            }
            t2.deleteSubtree( t2.getNode( "def" ), false );
            t2.toNewHampshireX();
            if ( t2.getNumberOfExternalNodes() != 1 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSupportCount() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny t0_1 = factory.create( "(((A,B),C),(D,E))", new NHXParser() )[ 0 ];
            final Phylogeny[] phylogenies_1 = factory.create( "(((A,B),C),(D,E)) " + "(((C,B),A),(D,E))"
                                                                      + "(((A,B),C),(D,E)) " + "(((A,B),C),(D,E))"
                                                                      + "(((A,B),C),(D,E))" + "(((C,B),A),(D,E))"
                                                                      + "(((E,B),D),(C,A))" + "(((C,B),A),(D,E))"
                                                                      + "(((A,B),C),(D,E))" + "(((A,B),C),(D,E))",
                                                              new NHXParser() );
            SupportCount.count( t0_1, phylogenies_1, true, false );
            final Phylogeny t0_2 = factory.create( "(((((A,B),C),D),E),(F,G))", new NHXParser() )[ 0 ];
            final Phylogeny[] phylogenies_2 = factory.create( "(((((A,B),C),D),E),(F,G))"
                                                                      + "(((((A,B),C),D),E),((F,G),X))"
                                                                      + "(((((A,Y),B),C),D),((F,G),E))"
                                                                      + "(((((A,B),C),D),E),(F,G))"
                                                                      + "(((((A,B),C),D),E),(F,G))"
                                                                      + "(((((A,B),C),D),E),(F,G))"
                                                                      + "(((((A,B),C),D),E),(F,G),Z)"
                                                                      + "(((((A,B),C),D),E),(F,G))"
                                                                      + "((((((A,B),C),D),E),F),G)"
                                                                      + "(((((X,Y),F,G),E),((A,B),C)),D)",
                                                              new NHXParser() );
            SupportCount.count( t0_2, phylogenies_2, true, false );
            final PhylogenyNodeIterator it = t0_2.iteratorPostorder();
            while ( it.hasNext() ) {
                final PhylogenyNode n = it.next();
                if ( !n.isExternal() && ( PhylogenyMethods.getConfidenceValue( n ) != 10 ) ) {
                    return false;
                }
            }
            final Phylogeny t0_3 = factory.create( "(((A,B)ab,C)abc,((D,E)de,F)def)", new NHXParser() )[ 0 ];
            final Phylogeny[] phylogenies_3 = factory.create( "(((A,B),C),((D,E),F))" + "(((A,C),B),((D,F),E))"
                    + "(((C,A),B),((F,D),E))" + "(((A,B),F),((D,E),C))" + "(((((A,B),C),D),E),F)", new NHXParser() );
            SupportCount.count( t0_3, phylogenies_3, true, false );
            t0_3.reRoot( t0_3.getNode( "def" ).getNodeId() );
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "ab" ) ) != 3 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "abc" ) ) != 4 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "def" ) ) != 4 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "de" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "A" ) ) != 5 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "B" ) ) != 5 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "C" ) ) != 5 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "D" ) ) != 5 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "E" ) ) != 5 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_3.getNode( "F" ) ) != 5 ) {
                return false;
            }
            final Phylogeny t0_4 = factory.create( "(((((A,B)1,C)2,D)3,E)4,F)", new NHXParser() )[ 0 ];
            final Phylogeny[] phylogenies_4 = factory.create( "((((((A,X),C),B),D),E),F) "
                    + "(((A,B,Z),C,Q),(((D,Y),E),F))", new NHXParser() );
            SupportCount.count( t0_4, phylogenies_4, true, false );
            t0_4.reRoot( t0_4.getNode( "F" ).getNodeId() );
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "1" ) ) != 1 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "2" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "3" ) ) != 1 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "4" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "A" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "B" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "C" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "D" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "E" ) ) != 2 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( t0_4.getNode( "F" ) ) != 2 ) {
                return false;
            }
            Phylogeny a = factory.create( "(((((A,B)1,C)2,D)3,E)4,F)", new NHXParser() )[ 0 ];
            final Phylogeny b1 = factory.create( "(((((B,A)1,C)2,D)3,E)4,F)", new NHXParser() )[ 0 ];
            double d = SupportCount.compare( b1, a, true, true, true );
            if ( !Test.isEqual( d, 5.0 / 5.0 ) ) {
                return false;
            }
            a = factory.create( "(((((A,B)1,C)2,D)3,E)4,F)", new NHXParser() )[ 0 ];
            final Phylogeny b2 = factory.create( "(((((C,B)1,A)2,D)3,E)4,F)", new NHXParser() )[ 0 ];
            d = SupportCount.compare( b2, a, true, true, true );
            if ( !Test.isEqual( d, 4.0 / 5.0 ) ) {
                return false;
            }
            a = factory.create( "(((((A,B)1,C)2,D)3,E)4,F)", new NHXParser() )[ 0 ];
            final Phylogeny b3 = factory.create( "(((((F,C)1,A)2,B)3,D)4,E)", new NHXParser() )[ 0 ];
            d = SupportCount.compare( b3, a, true, true, true );
            if ( !Test.isEqual( d, 2.0 / 5.0 ) ) {
                return false;
            }
            a = factory.create( "(((((A,B)1,C)2,D)3,E)4,F)r", new NHXParser() )[ 0 ];
            final Phylogeny b4 = factory.create( "(((((F,C)1,A)2,B)3,D)4,E)r", new NHXParser() )[ 0 ];
            d = SupportCount.compare( b4, a, true, true, false );
            if ( !Test.isEqual( d, 1.0 / 5.0 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSupportTransfer() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny p1 = factory.create( "(((A,B)ab:97,C)abc:57,((D,E)de:10,(F,G)fg:50,(H,I)hi:64)defghi)",
                                                 new NHXParser() )[ 0 ];
            final Phylogeny p2 = factory
                    .create( "(((A:0.1,B:0.3)ab:0.4,C)abc:0.5,((D,E)de,(F,G)fg,(H,I)hi:0.59)defghi)", new NHXParser() )[ 0 ];
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "ab" ) ) >= 0.0 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "abc" ) ) >= 0.0 ) {
                return false;
            }
            support_transfer.moveBranchLengthsToBootstrap( p1 );
            support_transfer.transferSupportValues( p1, p2 );
            if ( p2.getNode( "ab" ).getDistanceToParent() != 0.4 ) {
                return false;
            }
            if ( p2.getNode( "abc" ).getDistanceToParent() != 0.5 ) {
                return false;
            }
            if ( p2.getNode( "hi" ).getDistanceToParent() != 0.59 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "ab" ) ) != 97 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "abc" ) ) != 57 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "de" ) ) != 10 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "fg" ) ) != 50 ) {
                return false;
            }
            if ( PhylogenyMethods.getConfidenceValue( p2.getNode( "hi" ) ) != 64 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testTaxonomyAssigner() {
        try {
            String s0_str = "(((([&&NHX:S=A],[&&NHX:S=B])[&&NHX:S=AB],[&&NHX:S=C])[&&NHX:S=ABC],[&&NHX:S=D])[&&NHX:S=ABCD],[&&NHX:S=E])[&&NHX:S=ABCDE]";
            String g0_str = "((([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=A])a,[&&NHX:S=B])b,[&&NHX:S=C])c";
            Phylogeny s0 = ParserBasedPhylogenyFactory.getInstance().create( s0_str, new NHXParser() )[ 0 ];
            Phylogeny g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            s0.setRooted( true );
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "A" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "AB" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "ABC" ) ) {
                return false;
            }
            g0_str = "((([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=A])a,[&&NHX:S=A])b,[&&NHX:S=A])c";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "A" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "A" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "A" ) ) {
                return false;
            }
            g0_str = "((([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=B])a,[&&NHX:S=A])b,[&&NHX:S=A])c";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "AB" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "AB" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "AB" ) ) {
                return false;
            }
            g0_str = "((([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=B])a,[&&NHX:S=C])b,[&&NHX:S=A])c";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "AB" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "ABC" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "ABC" ) ) {
                return false;
            }
            g0_str = "((([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=B])a,[&&NHX:S=C])b,[&&NHX:S=D])c";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "AB" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "ABC" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            g0_str = "((([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=E])a,[&&NHX:S=C])b,[&&NHX:S=D])c";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCDE" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCDE" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCDE" ) ) {
                return false;
            }
            g0_str = "((([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=E])a,[&&NHX:S=A])b,[&&NHX:S=A])c";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCDE" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCDE" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCDE" ) ) {
                return false;
            }
            s0_str = "(([&&NHX:S=A],[&&NHX:S=B],[&&NHX:S=C],[&&NHX:S=D])[&&NHX:S=ABCD],"
                    + "([&&NHX:S=E],[&&NHX:S=F],[&&NHX:S=G],[&&NHX:S=H])[&&NHX:S=EFGH],"
                    + "([&&NHX:S=I],[&&NHX:S=J],[&&NHX:S=K],[&&NHX:S=L])[&&NHX:S=IJKL], "
                    + "([&&NHX:S=M],[&&NHX:S=N],[&&NHX:S=O],[&&NHX:S=P])[&&NHX:S=MNOP])[&&NHX:S=ROOT]";
            s0 = ParserBasedPhylogenyFactory.getInstance().create( s0_str, new NHXParser() )[ 0 ];
            s0.setRooted( true );
            g0_str = "(([&&NHX:S=A],[&&NHX:S=B],[&&NHX:S=C],[&&NHX:S=D])a,"
                    + "([&&NHX:S=E],[&&NHX:S=F],[&&NHX:S=G],[&&NHX:S=H])b,"
                    + "([&&NHX:S=I],[&&NHX:S=J],[&&NHX:S=K],[&&NHX:S=L])c, "
                    + "([&&NHX:S=M],[&&NHX:S=N],[&&NHX:S=O],[&&NHX:S=P])d)r";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "EFGH" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "IJKL" ) ) {
                return false;
            }
            if ( !g0.getNode( "d" ).getNodeData().getTaxonomy().getScientificName().equals( "MNOP" ) ) {
                return false;
            }
            if ( !g0.getNode( "r" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            g0_str = "(([&&NHX:S=A],[&&NHX:S=B],[&&NHX:S=A],[&&NHX:S=B])a,"
                    + "([&&NHX:S=E],[&&NHX:S=F],[&&NHX:S=F],[&&NHX:S=F])b,"
                    + "([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=I])c, "
                    + "([&&NHX:S=M],[&&NHX:S=N],[&&NHX:S=O],[&&NHX:S=O])d)r";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "EFGH" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "IJKL" ) ) {
                return false;
            }
            if ( !g0.getNode( "d" ).getNodeData().getTaxonomy().getScientificName().equals( "MNOP" ) ) {
                return false;
            }
            if ( !g0.getNode( "r" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            g0_str = "(([&&NHX:S=A],[&&NHX:S=B],[&&NHX:S=A],[&&NHX:S=B])a,"
                    + "([&&NHX:S=E],[&&NHX:S=F],[&&NHX:S=F],[&&NHX:S=F])b,"
                    + "([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L])c, "
                    + "([&&NHX:S=M],[&&NHX:S=N],[&&NHX:S=A],[&&NHX:S=O])d)r";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "EFGH" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "L" ) ) {
                return false;
            }
            if ( !g0.getNode( "d" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            if ( !g0.getNode( "r" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            g0_str = "(([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L])a,"
                    + "([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L])b,"
                    + "([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L])c, "
                    + "([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L])d)r";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "L" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "L" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "L" ) ) {
                return false;
            }
            if ( !g0.getNode( "d" ).getNodeData().getTaxonomy().getScientificName().equals( "L" ) ) {
                return false;
            }
            if ( !g0.getNode( "r" ).getNodeData().getTaxonomy().getScientificName().equals( "L" ) ) {
                return false;
            }
            g0_str = "((([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=A])a,[&&NHX:S=A])b,[&&NHX:S=A])c";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "A" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "A" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "A" ) ) {
                return false;
            }
            g0_str = "((([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=B])a,[&&NHX:S=I])b,[&&NHX:S=J])c";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            g0_str = "(((([&&NHX:S=A],[&&NHX:S=B],[&&NHX:S=C],[&&NHX:S=D])a,"
                    + "([&&NHX:S=D],[&&NHX:S=C],[&&NHX:S=B],[&&NHX:S=A])b)ab,"
                    + "([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L])c)abc, "
                    + "([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L])d)r";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "ab" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "L" ) ) {
                return false;
            }
            if ( !g0.getNode( "abc" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            if ( !g0.getNode( "d" ).getNodeData().getTaxonomy().getScientificName().equals( "L" ) ) {
                return false;
            }
            if ( !g0.getNode( "r" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            g0_str = "(((([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=C],[&&NHX:S=D])a,"
                    + "([&&NHX:S=D],[&&NHX:S=D],[&&NHX:S=B],[&&NHX:S=A])b)ab,"
                    + "([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L])c)abc, "
                    + "([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L])d)r";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "ab" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "L" ) ) {
                return false;
            }
            if ( !g0.getNode( "abc" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            if ( !g0.getNode( "d" ).getNodeData().getTaxonomy().getScientificName().equals( "L" ) ) {
                return false;
            }
            if ( !g0.getNode( "r" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            g0_str = "(((([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=C],[&&NHX:S=D])a,"
                    + "([&&NHX:S=D],[&&NHX:S=D],[&&NHX:S=B],[&&NHX:S=A])b)ab,"
                    + "([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L])c)abc, "
                    + "([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=A])d)r";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "ab" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "L" ) ) {
                return false;
            }
            if ( !g0.getNode( "abc" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            if ( !g0.getNode( "d" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            if ( !g0.getNode( "r" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            g0_str = "(((([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=C],[&&NHX:S=D])a,"
                    + "([&&NHX:S=D],[&&NHX:S=D],[&&NHX:S=B],[&&NHX:S=A])b)ab,"
                    + "([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=A])c)abc, "
                    + "([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=A])d)r";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( !g0.getNode( "a" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "b" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "ab" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "A" ) ) {
                return false;
            }
            if ( !g0.getNode( "abc" ).getNodeData().getTaxonomy().getScientificName().equals( "ABCD" ) ) {
                return false;
            }
            if ( !g0.getNode( "d" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            if ( !g0.getNode( "r" ).getNodeData().getTaxonomy().getScientificName().equals( "ROOT" ) ) {
                return false;
            }
            s0_str = "(([&&NHX:S=A],[&&NHX:S=B],[&&NHX:S=C],[&&NHX:S=D]),"
                    + "([&&NHX:S=E],[&&NHX:S=F],[&&NHX:S=G],[&&NHX:S=H]),"
                    + "([&&NHX:S=I],[&&NHX:S=J],[&&NHX:S=K],[&&NHX:S=L]), "
                    + "([&&NHX:S=M],[&&NHX:S=N],[&&NHX:S=O],[&&NHX:S=P]))";
            s0 = ParserBasedPhylogenyFactory.getInstance().create( s0_str, new NHXParser() )[ 0 ];
            s0.setRooted( true );
            g0_str = "(((([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=C],[&&NHX:S=D])a,"
                    + "([&&NHX:S=D],[&&NHX:S=D],[&&NHX:S=B],[&&NHX:S=A])b)ab,"
                    + "([&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=A],[&&NHX:S=A])c)abc, "
                    + "([&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=L],[&&NHX:S=A])d)r";
            g0 = ParserBasedPhylogenyFactory.getInstance().create( g0_str, new NHXParser() )[ 0 ];
            g0.setRooted( true );
            TaxonomyAssigner.execute( g0, s0 );
            if ( g0.getNode( "a" ).getNodeData().isHasTaxonomy() ) {
                return false;
            }
            if ( !g0.getNode( "c" ).getNodeData().getTaxonomy().getScientificName().equals( "A" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
}
