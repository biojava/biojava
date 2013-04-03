// $Id: TestSurfacing.java,v 1.83 2009/11/11 01:36:50 cmzmasek Exp $
//
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

package org.forester.surfacing;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.io.parsers.HmmPfamOutputParser;
import org.forester.io.parsers.nexus.PaupLogParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogenyinference.BasicCharacterStateMatrix;
import org.forester.phylogenyinference.CharacterStateMatrix;
import org.forester.phylogenyinference.CharacterStateMatrix.BinaryStates;
import org.forester.phylogenyinference.CharacterStateMatrix.GainLossStates;
import org.forester.surfacing.BinaryDomainCombination.DomainCombinationType;
import org.forester.test.Test;
import org.forester.util.ForesterUtil;

@SuppressWarnings( "unused")
public class TestSurfacing {

    private final static double ZERO_DIFF = 1.0E-9;

    public static boolean isEqual( final double a, final double b ) {
        return ( ( Math.abs( a - b ) ) < TestSurfacing.ZERO_DIFF );
    }

    private static StringBuffer mapToStringBuffer( final Map<PhylogenyNode, CharacterStateMatrix.BinaryStates> map ) {
        final StringBuffer sb = new StringBuffer();
        for( final PhylogenyNode key : map.keySet() ) {
            if ( !key.isExternal() ) {
                sb.append( key.getNodeName() );
                sb.append( " : " );
                sb.append( map.get( key ).toString() );
                sb.append( ForesterUtil.getLineSeparator() );
            }
        }
        return sb;
    }

    public static boolean test( final File test_dir ) {
        System.out.print( "  Domain id: " );
        if ( !TestSurfacing.testDomainId() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Protein id: " );
        if ( !TestSurfacing.testProteinId() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Species: " );
        if ( !TestSurfacing.testSpecies() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Basic domain: " );
        if ( !TestSurfacing.testBasicDomain() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Basic protein: " );
        if ( !TestSurfacing.testBasicProtein() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Combinable domains: " );
        if ( !TestSurfacing.testCombinableDomains() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Directed combinable domains: " );
        if ( !TestSurfacing.testDirectedCombinableDomains() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Genome wide specific combinable domains: " );
        if ( !TestSurfacing.testGenomeWideCombinableDomains() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Domain architecture based genome similarity calculator: " );
        if ( !TestSurfacing.testDomainArchitectureBasedGenomeSimilarityCalculator() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Hmmpfam output parser: " );
        if ( !TestSurfacing.testHmmPfamOutputParser( test_dir ) ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Hmmpfam output parser with filter: " );
        if ( !TestSurfacing.testHmmPfamOutputParserWithFilter( test_dir ) ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Combinations based pairwise similarity calculator: " );
        if ( !TestSurfacing.testCombinationsBasedPairwiseSimilarityCalculator() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Copy number based pairwise similarity calculator: " );
        if ( !TestSurfacing.testCopyNumberBasedPairwiseSimilarityCalculator() ) {
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Domain combination counting: " );
        if ( !TestSurfacing.testDomainCombinationCounting( test_dir ) ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Basic domain similarity calculator: " );
        if ( !TestSurfacing.testBasicDomainSimilarityCalculator() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Basic domain similarity calculator not ignoring species specific domains: " );
        if ( !TestSurfacing.testBasicDomainSimilarityCalculatorNotIgnoringSpeciesSpeficDomains() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Basic domain similarity calculator removal of singles: " );
        if ( !TestSurfacing.testBasicDomainSimilarityCalculatorRemovalOfSingles() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Domain sorting: " );
        if ( !TestSurfacing.testDomainSorting() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Overlap removal: " );
        if ( !TestSurfacing.testOverlapRemoval() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Engulfing overlap removal: " );
        if ( !TestSurfacing.testEngulfingOverlapRemoval() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Binary domain combination: " );
        if ( !TestSurfacing.testBinaryDomainCombination() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Parsimony: " );
        if ( !TestSurfacing.testParsimony() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Directedness: " );
        if ( !TestSurfacing.testDirectedness() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Directedness and adjacency: " );
        if ( !TestSurfacing.testDirectednessAndAdjacency() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Dollo parsimony on secodary features: " );
        if ( !TestSurfacing.testParsimonyOnSecondaryFeatures() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Paup log parser: " );
        if ( !TestSurfacing.testPaupLogParser( test_dir ) ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Binary state matrix to gain loss matrix: " );
        if ( !TestSurfacing.testBinaryStateMatrixToGainLossMatrix( test_dir ) ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        return true;
    }

    private static boolean testBasicDomain() {
        try {
            final Domain pd = new BasicDomain( "id", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            if ( !pd.getDomainId().getId().equals( "id" ) ) {
                return false;
            }
            if ( pd.getNumber() != 1 ) {
                return false;
            }
            if ( pd.getTotalCount() != 4 ) {
                return false;
            }
            if ( !pd.equals( new BasicDomain( "id", 22, 111, ( short ) 1, ( short ) 4, 0.2, -12 ) ) ) {
                return false;
            }
            final Domain a1 = new BasicDomain( "a", 1, 10, ( short ) 1, ( short ) 4, 0.1, -12 );
            final BasicDomain a1_copy = new BasicDomain( "a", 1, 10, ( short ) 1, ( short ) 4, 0.1, -12 );
            final BasicDomain a1_equal = new BasicDomain( "a", 524, 743994, ( short ) 1, ( short ) 300, 3.0005, 230 );
            final BasicDomain a2 = new BasicDomain( "a", 1, 10, ( short ) 2, ( short ) 4, 0.1, -12 );
            final BasicDomain a3 = new BasicDomain( "A", 1, 10, ( short ) 1, ( short ) 4, 0.1, -12 );
            if ( !a1.equals( a1 ) ) {
                return false;
            }
            if ( !a1.equals( a1_copy ) ) {
                return false;
            }
            if ( !a1.equals( a1_equal ) ) {
                return false;
            }
            if ( !a1.equals( a2 ) ) {
                return false;
            }
            if ( a1.equals( a3 ) ) {
                return false;
            }
            if ( a1.compareTo( a1 ) != 0 ) {
                return false;
            }
            if ( a1.compareTo( a1_copy ) != 0 ) {
                return false;
            }
            if ( a1.compareTo( a1_equal ) != 0 ) {
                return false;
            }
            if ( a1.compareTo( a2 ) != 0 ) {
                return false;
            }
            if ( a1.compareTo( a3 ) != 0 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicDomainSimilarityCalculator() {
        // mouse : ABCDE
        // rabbit: A.C.EF
        // ciona : A....FGX
        // nemve : ABCDEFG
        //
        // domain A:
        // m r c n
        // m 2/(2+3) 0 4/(4+2)
        // r 1/(1+4) 3/(3+3)
        // c 2/(2+5)
        // n
        //
        // mean = ( 2/5 + 0 + 2/3 + 1/5 + 1/2 + 2/7 ) / 6
        // min = 0.0
        // max = 2/3
        // n = 6
        //
        //
        // domain B:
        // m n
        // m 4/(4+2)
        // n
        //
        // mean = 2/3
        // min = 2/3
        // max = 2/3
        // sd = 0.0
        // n = 1
        //
        //
        // domain C:
        // m r n
        // m - 2/(2+3) 4/(4+2)
        // r - - 3/(3+3)
        // n - - -
        //
        // mean = (2/5 + 2/3 + 1/2)/3 =
        // min = 2/5
        // max = 2/3
        // sd = 0.0
        // n = 3
        try {
            final Domain A = new BasicDomain( "A", 1, 2, ( short ) 1, ( short ) 1, 0.15, -12 );
            final Domain B = new BasicDomain( "B", 1, 2, ( short ) 1, ( short ) 1, 0.2, -12 );
            final Domain C = new BasicDomain( "C", 1, 2, ( short ) 1, ( short ) 1, 0.3, -12 );
            final Domain D = new BasicDomain( "D", 1, 2, ( short ) 1, ( short ) 1, 0.5, -12 );
            final Domain E = new BasicDomain( "E", 1, 2, ( short ) 1, ( short ) 1, 0.5, -12 );
            final Domain F = new BasicDomain( "F", 1, 2, ( short ) 1, ( short ) 1, 0.01, -12 );
            final Domain G = new BasicDomain( "G", 1, 2, ( short ) 1, ( short ) 1, 0.001, -12 );
            final Domain X = new BasicDomain( "X", 1, 2, ( short ) 1, ( short ) 1, 0.0001, -12 );
            if ( !TestSurfacing.isEqual( X.getPerSequenceScore(), -12 ) ) {
                return false;
            }
            final Protein mouse_1 = new BasicProtein( "1", "mouse" );
            final Protein rabbit_1 = new BasicProtein( "1", "rabbit" );
            final Protein ciona_1 = new BasicProtein( "1", "ciona" );
            final Protein nemve_1 = new BasicProtein( "1", "nemve" );
            mouse_1.addProteinDomain( A );
            mouse_1.addProteinDomain( B );
            mouse_1.addProteinDomain( C );
            mouse_1.addProteinDomain( D );
            mouse_1.addProteinDomain( E );
            rabbit_1.addProteinDomain( A );
            rabbit_1.addProteinDomain( C );
            rabbit_1.addProteinDomain( E );
            rabbit_1.addProteinDomain( F );
            rabbit_1.addProteinDomain( F );
            rabbit_1.addProteinDomain( F );
            rabbit_1.addProteinDomain( F );
            rabbit_1.addProteinDomain( F );
            rabbit_1.addProteinDomain( F );
            ciona_1.addProteinDomain( A );
            ciona_1.addProteinDomain( A );
            ciona_1.addProteinDomain( A );
            ciona_1.addProteinDomain( A );
            ciona_1.addProteinDomain( A );
            ciona_1.addProteinDomain( F );
            ciona_1.addProteinDomain( G );
            ciona_1.addProteinDomain( X );
            nemve_1.addProteinDomain( A );
            nemve_1.addProteinDomain( B );
            nemve_1.addProteinDomain( C );
            nemve_1.addProteinDomain( D );
            nemve_1.addProteinDomain( E );
            nemve_1.addProteinDomain( F );
            nemve_1.addProteinDomain( G );
            final List<Protein> protein_list_mouse = new ArrayList<Protein>();
            final List<Protein> protein_list_rabbit = new ArrayList<Protein>();
            final List<Protein> protein_list_ciona = new ArrayList<Protein>();
            final List<Protein> protein_list_nemve = new ArrayList<Protein>();
            protein_list_mouse.add( mouse_1 );
            protein_list_rabbit.add( rabbit_1 );
            protein_list_ciona.add( ciona_1 );
            protein_list_nemve.add( nemve_1 );
            final List<GenomeWideCombinableDomains> cdc_list = new ArrayList<GenomeWideCombinableDomains>();
            cdc_list.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_mouse,
                                                                           true,
                                                                           new BasicSpecies( "mouse" ) ) );
            cdc_list.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_rabbit,
                                                                           true,
                                                                           new BasicSpecies( "rabbit" ) ) );
            cdc_list.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_ciona,
                                                                           true,
                                                                           new BasicSpecies( "ciona" ) ) );
            cdc_list.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_nemve,
                                                                           true,
                                                                           new BasicSpecies( "nemve" ) ) );
            final DomainSimilarityCalculator calc = new BasicDomainSimilarityCalculator( DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID,
                                                                                         false,
                                                                                         false );
            final SortedSet<DomainSimilarity> sims = calc
                    .calculateSimilarities( new CombinationsBasedPairwiseDomainSimilarityCalculator(),
                                            cdc_list,
                                            true,
                                            true );
            final Iterator<DomainSimilarity> sims_it = sims.iterator();
            final DomainSimilarity sa = sims_it.next();
            if ( !sa.getDomainId().getId().equals( "A" ) ) {
                return false;
            }
            if ( sa.getSpeciesData().size() != 4 ) {
                return false;
            }
            if ( !sa.getSpecies().contains( new BasicSpecies( "ciona" ) ) ) {
                return false;
            }
            if ( !sa.getSpecies().contains( new BasicSpecies( "mouse" ) ) ) {
                return false;
            }
            if ( !sa.getSpecies().contains( new BasicSpecies( "nemve" ) ) ) {
                return false;
            }
            if ( !sa.getSpecies().contains( new BasicSpecies( "rabbit" ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa.getMeanSimilarityScore(),
                                         ( 2.0 / 5 + 0 + 2.0 / 3 + 1.0 / 5 + 1.0 / 2 + 2.0 / 7 ) / 6 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa.getStandardDeviationOfSimilarityScore(), ( 0.23410788192183737 ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa.getMaximalSimilarityScore(), ( 2.0 / 3 ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa.getMinimalSimilarityScore(), ( 0.0 ) ) ) {
                return false;
            }
            if ( sa.getN() != 6 ) {
                return false;
            }
            if ( sa.getMaximalDifference() != 7 ) {
                return false;
            }
            if ( sa.getMaximalDifferenceInCounts() != 3 ) {
                return false;
            }
            final DomainSimilarity sb = sims_it.next();
            if ( !sb.getDomainId().getId().equals( "B" ) ) {
                return false;
            }
            if ( sb.getSpeciesData().size() != 2 ) {
                return false;
            }
            if ( !sb.getSpecies().contains( new BasicSpecies( "mouse" ) ) ) {
                return false;
            }
            if ( !sb.getSpecies().contains( new BasicSpecies( "nemve" ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sb.getMeanSimilarityScore(), 2.0 / 3 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sb.getStandardDeviationOfSimilarityScore(), 0.0 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sb.getMaximalSimilarityScore(), ( 2.0 / 3 ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sb.getMinimalSimilarityScore(), ( 2.0 / 3 ) ) ) {
                return false;
            }
            if ( sb.getN() != 1 ) {
                return false;
            }
            if ( sb.getMaximalDifference() != 2 ) {
                return false;
            }
            if ( sb.getMaximalDifferenceInCounts() != 2 ) {
                return false;
            }
            final DomainSimilarity sc = sims_it.next();
            if ( !sc.getDomainId().getId().equals( "C" ) ) {
                return false;
            }
            if ( sc.getSpeciesData().size() != 3 ) {
                return false;
            }
            if ( !sc.getSpecies().contains( new BasicSpecies( "mouse" ) ) ) {
                return false;
            }
            if ( !sc.getSpecies().contains( new BasicSpecies( "rabbit" ) ) ) {
                return false;
            }
            if ( !sc.getSpecies().contains( new BasicSpecies( "nemve" ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sc.getMeanSimilarityScore(), ( 2.0 / 5 + 2.0 / 3 + 1.0 / 2 ) / 3 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sc.getStandardDeviationOfSimilarityScore(), 0.13471506281091264 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sc.getMaximalSimilarityScore(), ( 2.0 / 3 ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sc.getMinimalSimilarityScore(), ( 2.0 / 5 ) ) ) {
                return false;
            }
            if ( sc.getN() != 3 ) {
                return false;
            }
            if ( sc.getMaximalDifference() != 3 ) {
                return false;
            }
            if ( sc.getMaximalDifferenceInCounts() != 3 ) {
                return false;
            }
            // mouse : ....ABCDE.....
            // rabbit: ....A.C.EFFF..
            // ciona : AAAAA......FGX
            // nemve : ....ABCDEFG...
            //
            // domain A:
            // m r c n
            // m 2/(2+3) 0 4/(4+2)
            // r - 1/(1+5) 3/(3+3)
            // c - 2/(2+6)
            // n
            //
            // mean = ( 2/5 + 0 + 2/3 + 1/6 + 1/2 + 2/8 ) / 6
            // min = 0.0
            // max = 2/3
            // n = 6
            final List<GenomeWideCombinableDomains> cdc_list2 = new ArrayList<GenomeWideCombinableDomains>();
            cdc_list2.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_mouse,
                                                                            false,
                                                                            new BasicSpecies( "mouse" ) ) );
            cdc_list2.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_rabbit,
                                                                            false,
                                                                            new BasicSpecies( "rabbit" ) ) );
            cdc_list2.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_ciona,
                                                                            false,
                                                                            new BasicSpecies( "ciona" ) ) );
            cdc_list2.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_nemve,
                                                                            false,
                                                                            new BasicSpecies( "nemve" ) ) );
            final DomainSimilarityCalculator calc2 = new BasicDomainSimilarityCalculator( DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID,
                                                                                          false,
                                                                                          false );
            final SortedSet<DomainSimilarity> sims2 = calc2
                    .calculateSimilarities( new CombinationsBasedPairwiseDomainSimilarityCalculator(),
                                            cdc_list2,
                                            false,
                                            true );
            final Iterator<DomainSimilarity> sims_it2 = sims2.iterator();
            final DomainSimilarity sa2 = sims_it2.next();
            if ( !sa2.getDomainId().getId().equals( "A" ) ) {
                return false;
            }
            if ( sa2.getSpeciesData().size() != 4 ) {
                return false;
            }
            if ( !sa2.getSpecies().contains( new BasicSpecies( "ciona" ) ) ) {
                return false;
            }
            if ( !sa2.getSpecies().contains( new BasicSpecies( "mouse" ) ) ) {
                return false;
            }
            if ( !sa2.getSpecies().contains( new BasicSpecies( "nemve" ) ) ) {
                return false;
            }
            if ( !sa2.getSpeciesData().keySet().contains( new BasicSpecies( "rabbit" ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa2.getMeanSimilarityScore(),
                                         ( 2.0 / 5 + 0 + 2.0 / 3 + 1.0 / 6 + 1.0 / 2 + 2.0 / 8 ) / 6 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa2.getStandardDeviationOfSimilarityScore(), ( 0.2404663678647683 ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa2.getMaximalSimilarityScore(), ( 2.0 / 3 ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa2.getMinimalSimilarityScore(), ( 0.0 ) ) ) {
                return false;
            }
            if ( sa2.getN() != 6 ) {
                return false;
            }
            if ( sa2.getMaximalDifference() != 8 ) {
                return false;
            }
            if ( sa2.getMaximalDifferenceInCounts() != 3 ) {
                return false;
            }
            final Protein ciona_2 = new BasicProtein( "2", "ciona" );
            ciona_2.addProteinDomain( A );
            ciona_2.addProteinDomain( A );
            ciona_2.addProteinDomain( A );
            ciona_2.addProteinDomain( B );
            ciona_2.addProteinDomain( B );
            ciona_2.addProteinDomain( B );
            ciona_2.addProteinDomain( F );
            ciona_2.addProteinDomain( F );
            ciona_2.addProteinDomain( F );
            ciona_2.addProteinDomain( F );
            ciona_2.addProteinDomain( G );
            ciona_2.addProteinDomain( X );
            final Protein ciona_3 = new BasicProtein( "3", "ciona" );
            ciona_3.addProteinDomain( A );
            ciona_3.addProteinDomain( A );
            ciona_3.addProteinDomain( A );
            ciona_3.addProteinDomain( A );
            ciona_3.addProteinDomain( B );
            ciona_3.addProteinDomain( B );
            ciona_3.addProteinDomain( X );
            ciona_3.addProteinDomain( X );
            protein_list_ciona.add( ciona_2 );
            protein_list_ciona.add( ciona_3 );
            final List<GenomeWideCombinableDomains> cdc_list3 = new ArrayList<GenomeWideCombinableDomains>();
            cdc_list3.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_mouse,
                                                                            true,
                                                                            new BasicSpecies( "mouse" ) ) );
            cdc_list3.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_rabbit,
                                                                            true,
                                                                            new BasicSpecies( "rabbit" ) ) );
            cdc_list3.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_ciona,
                                                                            true,
                                                                            new BasicSpecies( "ciona" ) ) );
            cdc_list3.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_nemve,
                                                                            true,
                                                                            new BasicSpecies( "nemve" ) ) );
            final DomainSimilarityCalculator calc3 = new BasicDomainSimilarityCalculator( DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID,
                                                                                          false,
                                                                                          false );
            final SortedSet<DomainSimilarity> sims3 = calc3
                    .calculateSimilarities( new CombinationsBasedPairwiseDomainSimilarityCalculator(),
                                            cdc_list3,
                                            false,
                                            true );
            final Iterator<DomainSimilarity> sims_it3 = sims3.iterator();
            final DomainSimilarity sa3 = sims_it3.next();
            if ( !sa3.getDomainId().getId().equals( "A" ) ) {
                return false;
            }
            final SpeciesSpecificDomainSimilariyData ssdsd = sa3.getSpeciesData().get( new BasicSpecies( "ciona" ) );
            if ( ssdsd.getCombinableDomainIdToCountsMap().size() != 4 ) {
                return false;
            }
            if ( ssdsd.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "B" ) ) != 2 ) {
                return false;
            }
            if ( ssdsd.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "F" ) ) != 2 ) {
                return false;
            }
            if ( ssdsd.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "G" ) ) != 2 ) {
                return false;
            }
            if ( ssdsd.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "X" ) ) != 3 ) {
                return false;
            }
            final List<GenomeWideCombinableDomains> cdc_list4 = new ArrayList<GenomeWideCombinableDomains>();
            cdc_list4.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_mouse,
                                                                            false,
                                                                            new BasicSpecies( "mouse" ) ) );
            cdc_list4.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_rabbit,
                                                                            false,
                                                                            new BasicSpecies( "rabbit" ) ) );
            cdc_list4.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_ciona,
                                                                            false,
                                                                            new BasicSpecies( "ciona" ) ) );
            ;
            cdc_list4.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_nemve,
                                                                            false,
                                                                            new BasicSpecies( "nemve" ) ) );
            final DomainSimilarityCalculator calc4 = new BasicDomainSimilarityCalculator( DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID,
                                                                                          true,
                                                                                          false );
            final SortedSet<DomainSimilarity> sims4 = calc4
                    .calculateSimilarities( new CombinationsBasedPairwiseDomainSimilarityCalculator(),
                                            cdc_list4,
                                            false,
                                            true );
            final Iterator<DomainSimilarity> sims_it4 = sims4.iterator();
            final DomainSimilarity sa4 = sims_it4.next();
            if ( !sa4.getDomainId().getId().equals( "A" ) ) {
                return false;
            }
            final SpeciesSpecificDomainSimilariyData ssdsd4 = sa4.getSpeciesData().get( new BasicSpecies( "ciona" ) );
            if ( ssdsd4.getCombinableDomainIdToCountsMap().size() != 5 ) {
                return false;
            }
            if ( ssdsd4.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "A" ) ) != 3 ) {
                return false;
            }
            if ( ssdsd4.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "B" ) ) != 2 ) {
                return false;
            }
            if ( ssdsd4.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "F" ) ) != 2 ) {
                return false;
            }
            if ( ssdsd4.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "G" ) ) != 2 ) {
                return false;
            }
            if ( ssdsd4.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "X" ) ) != 3 ) {
                return false;
            }
            final SortedSet<DomainSimilarity> sims4_d = calc4
                    .calculateSimilarities( new DomainCountsBasedPairwiseSimilarityCalculator(), cdc_list4, false, true );
            final Iterator<DomainSimilarity> sims_it4_d = sims4_d.iterator();
            final DomainSimilarity sa4_d = sims_it4_d.next();
            if ( !sa4_d.getDomainId().getId().equals( "A" ) ) {
                return false;
            }
            if ( sa4_d.getCombinableDomainIds( new BasicSpecies( "ciona" ) ).size() != 5 ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa4_d.getMeanSimilarityScore(), ( 1 + 1 - 11.0 / 13 + 1 - 11.0 / 13 + 1 + 1
                    + 1 - 11.0 / 13 ) / 6.0 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa4_d.getMaximalSimilarityScore(), 1.0 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa4_d.getMinimalSimilarityScore(), ( 1 - 11.0 / 13 ) ) ) {
                return false;
            }
            if ( sa4_d.getN() != 6 ) {
                return false;
            }
            final SortedSet<DomainSimilarity> sims4_p = calc4
                    .calculateSimilarities( new ProteinCountsBasedPairwiseDomainSimilarityCalculator(),
                                            cdc_list4,
                                            false,
                                            true );
            final Iterator<DomainSimilarity> sims_it4_p = sims4_p.iterator();
            final DomainSimilarity sa4_p = sims_it4_p.next();
            if ( !sa4_p.getDomainId().getId().equals( "A" ) ) {
                return false;
            }
            if ( sa4_p.getCombinableDomainIds( new BasicSpecies( "ciona" ) ).size() != 5 ) {
                return false;
            }
            if ( !sa4_p.getCombinableDomainIds( new BasicSpecies( "ciona" ) ).contains( new DomainId( "A" ) ) ) {
                return false;
            }
            if ( !sa4_p.getCombinableDomainIds( new BasicSpecies( "ciona" ) ).contains( new DomainId( "B" ) ) ) {
                return false;
            }
            if ( !sa4_p.getCombinableDomainIds( new BasicSpecies( "ciona" ) ).contains( new DomainId( "F" ) ) ) {
                return false;
            }
            if ( !sa4_p.getCombinableDomainIds( new BasicSpecies( "ciona" ) ).contains( new DomainId( "G" ) ) ) {
                return false;
            }
            if ( !sa4_p.getCombinableDomainIds( new BasicSpecies( "ciona" ) ).contains( new DomainId( "X" ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa4_p.getMeanSimilarityScore(),
                                         ( 1 + 1 - 2.0 / 4 + 1 - 2.0 / 4 + 1 + 1 + 1 - 2.0 / 4 ) / 6.0 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa4_p.getMaximalSimilarityScore(), 1 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa4_p.getMinimalSimilarityScore(), ( 1 - 2.0 / 4 ) ) ) {
                return false;
            }
            if ( sa4_p.getN() != 6 ) {
                return false;
            }
            final List<GenomeWideCombinableDomains> cdc_list5 = new ArrayList<GenomeWideCombinableDomains>();
            cdc_list5.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_mouse,
                                                                            true,
                                                                            new BasicSpecies( "mouse" ) ) );
            cdc_list5.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_rabbit,
                                                                            true,
                                                                            new BasicSpecies( "rabbit" ) ) );
            cdc_list5.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_ciona,
                                                                            true,
                                                                            new BasicSpecies( "ciona" ) ) );
            cdc_list5.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_nemve,
                                                                            true,
                                                                            new BasicSpecies( "nemve" ) ) );
            final SortedSet<DomainSimilarity> sims5_d = calc4
                    .calculateSimilarities( new DomainCountsBasedPairwiseSimilarityCalculator(), cdc_list5, false, true );
            final Iterator<DomainSimilarity> sims_it5_d = sims5_d.iterator();
            final DomainSimilarity sa5_d = sims_it5_d.next();
            if ( sa5_d.getSpecies().size() != 4 ) {
                return false;
            }
            if ( !sa5_d.getSpecies().last().equals( new BasicSpecies( "rabbit" ) ) ) {
                return false;
            }
            final SpeciesSpecificDomainSimilariyData ssdsd5 = sa5_d.getSpeciesData().get( new BasicSpecies( "ciona" ) );
            if ( ssdsd5.getCombinableDomainIdToCountsMap().size() != 4 ) {
                return false;
            }
            if ( ssdsd5.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "B" ) ) != 2 ) {
                return false;
            }
            if ( ssdsd5.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "F" ) ) != 2 ) {
                return false;
            }
            if ( ssdsd5.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "G" ) ) != 2 ) {
                return false;
            }
            if ( ssdsd5.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "X" ) ) != 3 ) {
                return false;
            }
            if ( !sa5_d.getDomainId().getId().equals( "A" ) ) {
                return false;
            }
            final Species ciona = new BasicSpecies( "ciona" );
            if ( sa5_d.getCombinableDomainIds( ciona ).size() != 4 ) {
                return false;
            }
            if ( sa5_d.getCombinableDomainIds( ciona ).contains( new DomainId( "A" ) ) ) {
                return false;
            }
            if ( !sa5_d.getCombinableDomainIds( ciona ).contains( new DomainId( "B" ) ) ) {
                return false;
            }
            if ( !sa5_d.getCombinableDomainIds( ciona ).contains( new DomainId( "F" ) ) ) {
                return false;
            }
            if ( !sa5_d.getCombinableDomainIds( ciona ).contains( new DomainId( "G" ) ) ) {
                return false;
            }
            if ( !sa5_d.getCombinableDomainIds( ciona ).contains( new DomainId( "X" ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa5_d.getMeanSimilarityScore(), ( 1 + 1 - 11.0 / 13 + 1 - 11.0 / 13 + 1 + 1
                    + 1 - 11.0 / 13 ) / 6.0 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa5_d.getMaximalSimilarityScore(), 1.0 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa5_d.getMinimalSimilarityScore(), ( 1 - 11.0 / 13 ) ) ) {
                return false;
            }
            if ( sa5_d.getN() != 6 ) {
                return false;
            }
            if ( sa5_d.getMaximalDifference() != sa5_d.getMaximalDifferenceInCounts() ) {
                return false;
            }
            if ( sa5_d.getMaximalDifference() != 11 ) {
                return false;
            }
            if ( sa5_d.getMaximalDifferenceInCounts() != 11 ) {
                return false;
            }
            final SortedSet<DomainSimilarity> sims5_p = calc4
                    .calculateSimilarities( new ProteinCountsBasedPairwiseDomainSimilarityCalculator(),
                                            cdc_list5,
                                            false,
                                            true );
            final Iterator<DomainSimilarity> sims_it5_p = sims5_p.iterator();
            final DomainSimilarity sa5_p = sims_it5_p.next();
            if ( !sa5_p.getDomainId().getId().equals( "A" ) ) {
                return false;
            }
            if ( sa5_p.getCombinableDomainIds( ciona ).size() != 4 ) {
                return false;
            }
            if ( sa5_p.getCombinableDomainIds( ciona ).contains( new DomainId( "A" ) ) ) {
                return false;
            }
            if ( !sa5_p.getCombinableDomainIds( ciona ).contains( new DomainId( "B" ) ) ) {
                return false;
            }
            if ( !sa5_p.getCombinableDomainIds( ciona ).contains( new DomainId( "F" ) ) ) {
                return false;
            }
            if ( !sa5_p.getCombinableDomainIds( ciona ).contains( new DomainId( "G" ) ) ) {
                return false;
            }
            if ( !sa5_p.getCombinableDomainIds( ciona ).contains( new DomainId( "X" ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa5_p.getMeanSimilarityScore(),
                                         ( 1 + 1 - 2.0 / 4 + 1 - 2.0 / 4 + 1 + 1 + 1 - 2.0 / 4 ) / 6.0 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa5_p.getMaximalSimilarityScore(), 1 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa5_p.getMinimalSimilarityScore(), ( 1 - 2.0 / 4 ) ) ) {
                return false;
            }
            if ( sa5_p.getN() != 6 ) {
                return false;
            }
            if ( sa5_p.getMaximalDifference() != sa5_p.getMaximalDifferenceInCounts() ) {
                return false;
            }
            if ( sa5_p.getMaximalDifference() != 2 ) {
                return false;
            }
            if ( sa5_p.getMaximalDifferenceInCounts() != 2 ) {
                return false;
            }
            final List<GenomeWideCombinableDomains> cdc_list6 = new ArrayList<GenomeWideCombinableDomains>();
            cdc_list6.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_mouse,
                                                                            false,
                                                                            new BasicSpecies( "mouse" ) ) );
            cdc_list6.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_rabbit,
                                                                            false,
                                                                            new BasicSpecies( "rabbit" ) ) );
            cdc_list6.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_ciona,
                                                                            false,
                                                                            new BasicSpecies( "ciona" ) ) );
            cdc_list6.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_nemve,
                                                                            false,
                                                                            new BasicSpecies( "nemve" ) ) );
            final SortedSet<DomainSimilarity> sims6_d = calc4
                    .calculateSimilarities( new DomainCountsBasedPairwiseSimilarityCalculator(), cdc_list6, false, true );
            final Iterator<DomainSimilarity> sims_it6_d = sims6_d.iterator();
            final DomainSimilarity sa6_d = sims_it6_d.next();
            if ( sa6_d.getSpecies().size() != 4 ) {
                return false;
            }
            if ( !sa6_d.getSpecies().last().equals( new BasicSpecies( "rabbit" ) ) ) {
                return false;
            }
            final SpeciesSpecificDomainSimilariyData ssdsd6 = sa6_d.getSpeciesData().get( new BasicSpecies( "ciona" ) );
            if ( ssdsd6.getCombinableDomainIdToCountsMap().size() != 5 ) {
                return false;
            }
            if ( ssdsd6.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "B" ) ) != 2 ) {
                return false;
            }
            if ( ssdsd6.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "F" ) ) != 2 ) {
                return false;
            }
            if ( ssdsd6.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "G" ) ) != 2 ) {
                return false;
            }
            if ( ssdsd6.getNumberOfProteinsExhibitingCombinationWith( new DomainId( "X" ) ) != 3 ) {
                return false;
            }
            if ( !sa5_d.getDomainId().getId().equals( "A" ) ) {
                return false;
            }
            final Species ciona6 = new BasicSpecies( "ciona" );
            if ( sa6_d.getCombinableDomainIds( ciona6 ).size() != 5 ) {
                return false;
            }
            if ( !sa6_d.getCombinableDomainIds( ciona6 ).contains( new DomainId( "A" ) ) ) {
                return false;
            }
            if ( !sa6_d.getCombinableDomainIds( ciona6 ).contains( new DomainId( "B" ) ) ) {
                return false;
            }
            if ( !sa6_d.getCombinableDomainIds( ciona6 ).contains( new DomainId( "F" ) ) ) {
                return false;
            }
            if ( !sa6_d.getCombinableDomainIds( ciona6 ).contains( new DomainId( "G" ) ) ) {
                return false;
            }
            if ( !sa6_d.getCombinableDomainIds( ciona6 ).contains( new DomainId( "X" ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa6_d.getMeanSimilarityScore(), ( 1 + 1 - 11.0 / 13 + 1 - 11.0 / 13 + 1 + 1
                    + 1 - 11.0 / 13 ) / 6.0 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa6_d.getMaximalSimilarityScore(), 1.0 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa6_d.getMinimalSimilarityScore(), ( 1 - 11.0 / 13 ) ) ) {
                return false;
            }
            if ( sa6_d.getN() != 6 ) {
                return false;
            }
            if ( sa6_d.getMaximalDifference() != sa6_d.getMaximalDifferenceInCounts() ) {
                return false;
            }
            if ( sa6_d.getMaximalDifference() != 11 ) {
                return false;
            }
            if ( sa6_d.getMaximalDifferenceInCounts() != 11 ) {
                return false;
            }
            final SortedSet<DomainSimilarity> sims6_p = calc4
                    .calculateSimilarities( new ProteinCountsBasedPairwiseDomainSimilarityCalculator(),
                                            cdc_list6,
                                            false,
                                            true );
            final Iterator<DomainSimilarity> sims_it6_p = sims6_p.iterator();
            final DomainSimilarity sa6_p = sims_it6_p.next();
            if ( !sa6_p.getDomainId().getId().equals( "A" ) ) {
                return false;
            }
            if ( sa6_p.getCombinableDomainIds( ciona ).size() != 5 ) {
                return false;
            }
            if ( !sa6_p.getCombinableDomainIds( ciona ).contains( new DomainId( "A" ) ) ) {
                return false;
            }
            if ( !sa6_p.getCombinableDomainIds( ciona ).contains( new DomainId( "B" ) ) ) {
                return false;
            }
            if ( !sa6_p.getCombinableDomainIds( ciona ).contains( new DomainId( "F" ) ) ) {
                return false;
            }
            if ( !sa6_p.getCombinableDomainIds( ciona ).contains( new DomainId( "G" ) ) ) {
                return false;
            }
            if ( !sa6_p.getCombinableDomainIds( ciona ).contains( new DomainId( "X" ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa6_p.getMeanSimilarityScore(),
                                         ( 1 + 1 - 2.0 / 4 + 1 - 2.0 / 4 + 1 + 1 + 1 - 2.0 / 4 ) / 6.0 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa6_p.getMaximalSimilarityScore(), 1 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa6_p.getMinimalSimilarityScore(), ( 1 - 2.0 / 4 ) ) ) {
                return false;
            }
            if ( sa6_p.getN() != 6 ) {
                return false;
            }
            if ( sa6_p.getMaximalDifference() != sa6_p.getMaximalDifferenceInCounts() ) {
                return false;
            }
            if ( sa6_p.getMaximalDifference() != 2 ) {
                return false;
            }
            if ( sa6_p.getMaximalDifferenceInCounts() != 2 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicDomainSimilarityCalculatorNotIgnoringSpeciesSpeficDomains() {
        try {
            final Domain A = new BasicDomain( "A", 1, 2, ( short ) 1, ( short ) 1, 0.15, -12 );
            final Domain B = new BasicDomain( "B", 1, 2, ( short ) 1, ( short ) 1, 0.2, -12 );
            final Domain D = new BasicDomain( "D", 1, 2, ( short ) 1, ( short ) 1, 0.5, -12 );
            final Domain E = new BasicDomain( "E", 1, 2, ( short ) 1, ( short ) 1, 0.5, -12 );
            final Domain F = new BasicDomain( "F", 1, 2, ( short ) 1, ( short ) 1, 0.01, -12 );
            final Domain G = new BasicDomain( "G", 1, 2, ( short ) 1, ( short ) 1, 0.001, -12 );
            final Domain X = new BasicDomain( "X", 1, 2, ( short ) 1, ( short ) 1, 0.0001, -12 );
            if ( !TestSurfacing.isEqual( X.getPerSequenceScore(), -12 ) ) {
                return false;
            }
            final Protein mouse_1 = new BasicProtein( "1", "mouse" );
            final Protein rabbit_1 = new BasicProtein( "1", "rabbit" );
            final Protein ciona_1 = new BasicProtein( "1", "ciona" );
            final Protein nemve_1 = new BasicProtein( "1", "nemve" );
            mouse_1.addProteinDomain( A );
            mouse_1.addProteinDomain( D );
            mouse_1.addProteinDomain( E );
            rabbit_1.addProteinDomain( B );
            rabbit_1.addProteinDomain( E );
            rabbit_1.addProteinDomain( F );
            rabbit_1.addProteinDomain( F );
            rabbit_1.addProteinDomain( F );
            rabbit_1.addProteinDomain( F );
            rabbit_1.addProteinDomain( F );
            rabbit_1.addProteinDomain( F );
            ciona_1.addProteinDomain( F );
            ciona_1.addProteinDomain( G );
            ciona_1.addProteinDomain( X );
            nemve_1.addProteinDomain( D );
            nemve_1.addProteinDomain( E );
            nemve_1.addProteinDomain( F );
            nemve_1.addProteinDomain( G );
            final List<Protein> protein_list_mouse = new ArrayList<Protein>();
            final List<Protein> protein_list_rabbit = new ArrayList<Protein>();
            final List<Protein> protein_list_ciona = new ArrayList<Protein>();
            final List<Protein> protein_list_nemve = new ArrayList<Protein>();
            protein_list_mouse.add( mouse_1 );
            protein_list_rabbit.add( rabbit_1 );
            protein_list_ciona.add( ciona_1 );
            protein_list_nemve.add( nemve_1 );
            final List<GenomeWideCombinableDomains> cdc_list = new ArrayList<GenomeWideCombinableDomains>();
            cdc_list.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_mouse,
                                                                           true,
                                                                           new BasicSpecies( "mouse" ) ) );
            cdc_list.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_rabbit,
                                                                           true,
                                                                           new BasicSpecies( "rabbit" ) ) );
            cdc_list.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_ciona,
                                                                           true,
                                                                           new BasicSpecies( "ciona" ) ) );
            cdc_list.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_nemve,
                                                                           true,
                                                                           new BasicSpecies( "nemve" ) ) );
            final DomainSimilarityCalculator calc = new BasicDomainSimilarityCalculator( DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID,
                                                                                         false,
                                                                                         false );
            final SortedSet<DomainSimilarity> sims = calc
                    .calculateSimilarities( new CombinationsBasedPairwiseDomainSimilarityCalculator(),
                                            cdc_list,
                                            true,
                                            false );
            final Iterator<DomainSimilarity> sims_it = sims.iterator();
            final DomainSimilarity sa = sims_it.next();
            if ( !sa.getDomainId().getId().equals( "A" ) ) {
                return false;
            }
            if ( sa.getSpeciesData().size() != 1 ) {
                return false;
            }
            if ( !sa.getSpecies().contains( new BasicSpecies( "mouse" ) ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa.getMeanSimilarityScore(), 1.0 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa.getStandardDeviationOfSimilarityScore(), 0.0 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa.getMaximalSimilarityScore(), 1.0 ) ) {
                return false;
            }
            if ( !TestSurfacing.isEqual( sa.getMinimalSimilarityScore(), 1.0 ) ) {
                return false;
            }
            if ( sa.getN() != 0 ) {
                return false;
            }
            if ( sa.getMaximalDifference() != 0 ) {
                return false;
            }
            if ( sa.getMaximalDifferenceInCounts() != 0 ) {
                return false;
            }
            final DomainSimilarity sb = sims_it.next();
            if ( !sb.getDomainId().getId().equals( "B" ) ) {
                return false;
            }
            if ( sb.getSpeciesData().size() != 1 ) {
                return false;
            }
            if ( !sb.getSpecies().contains( new BasicSpecies( "rabbit" ) ) ) {
                return false;
            }
            final SortedSet<DomainSimilarity> sims2 = calc
                    .calculateSimilarities( new CombinationsBasedPairwiseDomainSimilarityCalculator(),
                                            cdc_list,
                                            true,
                                            true );
            final Iterator<DomainSimilarity> sims_it2 = sims2.iterator();
            final DomainSimilarity sa2 = sims_it2.next();
            if ( !sa2.getDomainId().getId().equals( "D" ) ) {
                return false;
            }
            if ( sa2.getSpeciesData().size() != 2 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicDomainSimilarityCalculatorRemovalOfSingles() {
        try {
            final Domain A = new BasicDomain( "A", 1, 2, ( short ) 1, ( short ) 1, 0.15, -12 );
            final Domain B = new BasicDomain( "B", 1, 2, ( short ) 1, ( short ) 1, 0.2, -12 );
            final Protein mouse_1 = new BasicProtein( "1", "mouse" );
            final Protein rabbit_1 = new BasicProtein( "1", "rabbit" );
            final Protein ciona_1 = new BasicProtein( "1", "ciona" );
            final Protein nemve_1 = new BasicProtein( "1", "nemve" );
            mouse_1.addProteinDomain( A );
            rabbit_1.addProteinDomain( A );
            ciona_1.addProteinDomain( A );
            ciona_1.addProteinDomain( A );
            ciona_1.addProteinDomain( A );
            ciona_1.addProteinDomain( A );
            ciona_1.addProteinDomain( A );
            nemve_1.addProteinDomain( A );
            final List<Protein> protein_list_mouse = new ArrayList<Protein>();
            final List<Protein> protein_list_rabbit = new ArrayList<Protein>();
            final List<Protein> protein_list_ciona = new ArrayList<Protein>();
            final List<Protein> protein_list_nemve = new ArrayList<Protein>();
            protein_list_mouse.add( mouse_1 );
            protein_list_rabbit.add( rabbit_1 );
            protein_list_ciona.add( ciona_1 );
            protein_list_nemve.add( nemve_1 );
            final List<GenomeWideCombinableDomains> cdc_list = new ArrayList<GenomeWideCombinableDomains>();
            cdc_list.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_mouse,
                                                                           true,
                                                                           new BasicSpecies( "mouse" ) ) );
            cdc_list.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_rabbit,
                                                                           true,
                                                                           new BasicSpecies( "rabbit" ) ) );
            cdc_list.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_ciona,
                                                                           true,
                                                                           new BasicSpecies( "ciona" ) ) );
            cdc_list.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_nemve,
                                                                           true,
                                                                           new BasicSpecies( "nemve" ) ) );
            final DomainSimilarityCalculator calc = new BasicDomainSimilarityCalculator( DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID,
                                                                                         false,
                                                                                         false );
            final SortedSet<DomainSimilarity> sims = calc
                    .calculateSimilarities( new CombinationsBasedPairwiseDomainSimilarityCalculator(),
                                            cdc_list,
                                            false,
                                            true );
            if ( sims.size() != 1 ) {
                return false;
            }
            final Iterator<DomainSimilarity> sims_it = sims.iterator();
            final DomainSimilarity sa = sims_it.next();
            if ( !sa.getDomainId().getId().equals( "A" ) ) {
                return false;
            }
            if ( sa.getSpeciesData().size() != 4 ) {
                return false;
            }
            if ( !sa.getSpecies().contains( new BasicSpecies( "ciona" ) ) ) {
                return false;
            }
            if ( !sa.getSpecies().contains( new BasicSpecies( "mouse" ) ) ) {
                return false;
            }
            if ( !sa.getSpecies().contains( new BasicSpecies( "nemve" ) ) ) {
                return false;
            }
            if ( !sa.getSpecies().contains( new BasicSpecies( "rabbit" ) ) ) {
                return false;
            }
            final SortedSet<DomainSimilarity> sims_ns = calc
                    .calculateSimilarities( new CombinationsBasedPairwiseDomainSimilarityCalculator(),
                                            cdc_list,
                                            true,
                                            true );
            if ( sims_ns.size() != 0 ) {
                return false;
            }
            final Protein mouse_2 = new BasicProtein( "1", "mouse" );
            final Protein rabbit_2 = new BasicProtein( "1", "rabbit" );
            final Protein ciona_2 = new BasicProtein( "1", "ciona" );
            final Protein nemve_2 = new BasicProtein( "1", "nemve" );
            mouse_2.addProteinDomain( A );
            rabbit_2.addProteinDomain( A );
            ciona_2.addProteinDomain( A );
            ciona_2.addProteinDomain( A );
            ciona_2.addProteinDomain( B );
            ciona_2.addProteinDomain( A );
            ciona_2.addProteinDomain( A );
            ciona_2.addProteinDomain( A );
            nemve_2.addProteinDomain( A );
            final List<Protein> protein_list_mouse2 = new ArrayList<Protein>();
            final List<Protein> protein_list_rabbit2 = new ArrayList<Protein>();
            final List<Protein> protein_list_ciona2 = new ArrayList<Protein>();
            final List<Protein> protein_list_nemve2 = new ArrayList<Protein>();
            protein_list_mouse2.add( mouse_2 );
            protein_list_rabbit2.add( rabbit_2 );
            protein_list_ciona2.add( ciona_2 );
            protein_list_nemve2.add( nemve_2 );
            final List<GenomeWideCombinableDomains> cdc_list2 = new ArrayList<GenomeWideCombinableDomains>();
            cdc_list2.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_mouse2,
                                                                            true,
                                                                            new BasicSpecies( "mouse" ) ) );
            cdc_list2.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_rabbit2,
                                                                            true,
                                                                            new BasicSpecies( "rabbit" ) ) );
            cdc_list2.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_ciona2,
                                                                            true,
                                                                            new BasicSpecies( "ciona" ) ) );
            cdc_list2.add( BasicGenomeWideCombinableDomains.createInstance( protein_list_nemve2,
                                                                            true,
                                                                            new BasicSpecies( "nemve" ) ) );
            final SortedSet<DomainSimilarity> sims2 = calc
                    .calculateSimilarities( new CombinationsBasedPairwiseDomainSimilarityCalculator(),
                                            cdc_list2,
                                            true,
                                            true );
            if ( sims2.size() != 1 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicProtein() {
        try {
            // A0  A10  B15  A20  B25  A30  B35  B40  C50  A60  C70  D80
            final Domain A0 = new BasicDomain( "A", 0, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain A10 = new BasicDomain( "A", 10, 11, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain B15 = new BasicDomain( "B", 11, 16, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain A20 = new BasicDomain( "A", 20, 100, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain B25 = new BasicDomain( "B", 25, 26, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain A30 = new BasicDomain( "A", 30, 31, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain B35 = new BasicDomain( "B", 31, 40, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain B40 = new BasicDomain( "B", 40, 600, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain C50 = new BasicDomain( "C", 50, 59, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain A60 = new BasicDomain( "A", 60, 395, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain C70 = new BasicDomain( "C", 70, 71, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain D80 = new BasicDomain( "D", 80, 81, ( short ) 1, ( short ) 4, 0.1, -12 );
            final BasicProtein p = new BasicProtein( "p", "owl" );
            p.addProteinDomain( B15 );
            p.addProteinDomain( C50 );
            p.addProteinDomain( A60 );
            p.addProteinDomain( A30 );
            p.addProteinDomain( C70 );
            p.addProteinDomain( B35 );
            p.addProteinDomain( B40 );
            p.addProteinDomain( A0 );
            p.addProteinDomain( A10 );
            p.addProteinDomain( A20 );
            p.addProteinDomain( B25 );
            p.addProteinDomain( D80 );
            List<DomainId> domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "C" ) );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( !p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids.add( new DomainId( "X" ) );
            if ( p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "D" ) );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( !p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "D" ) );
            domains_ids.add( new DomainId( "C" ) );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( !p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "B" ) );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( !p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "B" ) );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( !p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "B" ) );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "D" ) );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( !p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "D" ) );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( !p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "D" ) );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "D" ) );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "D" ) );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( !p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "D" ) );
            domains_ids.add( new DomainId( "X" ) );
            if ( p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "X" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "D" ) );
            if ( p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( p.contains( domains_ids, true ) ) {
                return false;
            }
            domains_ids = new ArrayList<DomainId>();
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "B" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "A" ) );
            domains_ids.add( new DomainId( "C" ) );
            domains_ids.add( new DomainId( "D" ) );
            if ( !p.contains( domains_ids, false ) ) {
                return false;
            }
            if ( p.contains( domains_ids, true ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBinaryDomainCombination() {
        try {
            final BasicBinaryDomainCombination s0 = new BasicBinaryDomainCombination( "a", "a" );
            final BasicBinaryDomainCombination s1 = new BasicBinaryDomainCombination( "b", "a" );
            final BasicBinaryDomainCombination s2 = new BasicBinaryDomainCombination( "a", "b" );
            final BasicBinaryDomainCombination s3 = new BasicBinaryDomainCombination( "B", "A" );
            final BasicBinaryDomainCombination s4 = new BasicBinaryDomainCombination( "A", "B" );
            final BasicBinaryDomainCombination s5 = new BasicBinaryDomainCombination( "c", "a" );
            final BasicBinaryDomainCombination s6 = new BasicBinaryDomainCombination( "b", "c" );
            final BasicBinaryDomainCombination s7 = new BasicBinaryDomainCombination( "d", "a" );
            final BasicBinaryDomainCombination s8 = new BasicBinaryDomainCombination( "b", "d" );
            final BinaryDomainCombination s9 = BasicBinaryDomainCombination.createInstance( " z-z=a-aa " );
            if ( !s9.toString().equals( "a-aa=z-z" ) ) {
                return false;
            }
            if ( !s0.equals( s0 ) ) {
                return false;
            }
            if ( s0.equals( s1 ) ) {
                return false;
            }
            if ( s1.equals( s0 ) ) {
                return false;
            }
            if ( !s1.equals( s2 ) ) {
                return false;
            }
            if ( !s2.equals( s1 ) ) {
                return false;
            }
            if ( s2.equals( s3 ) ) {
                return false;
            }
            if ( s2.equals( s3 ) ) {
                return false;
            }
            if ( s2.equals( s4 ) ) {
                return false;
            }
            final SortedSet<BasicBinaryDomainCombination> sorted = new TreeSet<BasicBinaryDomainCombination>();
            sorted.add( s0 );
            sorted.add( s1 );
            sorted.add( s2 );
            sorted.add( s3 );
            sorted.add( s3 );
            sorted.add( s3 );
            sorted.add( s4 );
            sorted.add( s5 );
            sorted.add( s6 );
            sorted.add( s7 );
            sorted.add( s7 );
            sorted.add( s8 );
            if ( sorted.size() != 6 ) {
                return false;
            }
            final DirectedBinaryDomainCombination aa = new DirectedBinaryDomainCombination( "a", "a" );
            final DirectedBinaryDomainCombination ba = new DirectedBinaryDomainCombination( "b", "a" );
            final DirectedBinaryDomainCombination ab = new DirectedBinaryDomainCombination( "a", "b" );
            final DirectedBinaryDomainCombination bb = new DirectedBinaryDomainCombination( "b", "b" );
            if ( !aa.equals( aa ) ) {
                return false;
            }
            if ( aa.equals( bb ) ) {
                return false;
            }
            if ( ab.equals( ba ) ) {
                return false;
            }
            if ( ba.equals( ab ) ) {
                return false;
            }
            if ( !ab.equals( ab ) ) {
                return false;
            }
            if ( ab.equals( aa ) ) {
                return false;
            }
            if ( ab.equals( bb ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBinaryStateMatrixToGainLossMatrix( final File test_dir ) {
        final BinaryStates I = BinaryStates.PRESENT;
        final BinaryStates O = BinaryStates.ABSENT;
        try {
            final CharacterStateMatrix<BinaryStates> binary_states_matrix_0 = new BasicCharacterStateMatrix<BinaryStates>( 7,
                                                                                                                           6 );
            binary_states_matrix_0.setIdentifier( 0, "A" );
            binary_states_matrix_0.setIdentifier( 1, "B" );
            binary_states_matrix_0.setIdentifier( 2, "C" );
            binary_states_matrix_0.setIdentifier( 3, "D" );
            binary_states_matrix_0.setIdentifier( 4, "1" );
            binary_states_matrix_0.setIdentifier( 5, "2" );
            binary_states_matrix_0.setIdentifier( 6, "3" );
            binary_states_matrix_0.setState( 0, 0, O );
            binary_states_matrix_0.setState( 1, 0, O );
            binary_states_matrix_0.setState( 2, 0, O );
            binary_states_matrix_0.setState( 3, 0, O );
            binary_states_matrix_0.setState( 4, 0, O );
            binary_states_matrix_0.setState( 5, 0, O );
            binary_states_matrix_0.setState( 6, 0, O );
            binary_states_matrix_0.setState( 0, 1, I );
            binary_states_matrix_0.setState( 1, 1, O );
            binary_states_matrix_0.setState( 2, 1, O );
            binary_states_matrix_0.setState( 3, 1, O );
            binary_states_matrix_0.setState( 4, 1, O );
            binary_states_matrix_0.setState( 5, 1, O );
            binary_states_matrix_0.setState( 6, 1, O );
            binary_states_matrix_0.setState( 0, 2, O );
            binary_states_matrix_0.setState( 1, 2, O );
            binary_states_matrix_0.setState( 2, 2, O );
            binary_states_matrix_0.setState( 3, 2, O );
            binary_states_matrix_0.setState( 4, 2, I );
            binary_states_matrix_0.setState( 5, 2, O );
            binary_states_matrix_0.setState( 6, 2, O );
            binary_states_matrix_0.setState( 0, 3, I );
            binary_states_matrix_0.setState( 1, 3, O );
            binary_states_matrix_0.setState( 2, 3, O );
            binary_states_matrix_0.setState( 3, 3, O );
            binary_states_matrix_0.setState( 4, 3, I );
            binary_states_matrix_0.setState( 5, 3, O );
            binary_states_matrix_0.setState( 6, 3, I );
            binary_states_matrix_0.setState( 0, 4, I );
            binary_states_matrix_0.setState( 1, 4, O );
            binary_states_matrix_0.setState( 2, 4, I );
            binary_states_matrix_0.setState( 3, 4, O );
            binary_states_matrix_0.setState( 4, 4, I );
            binary_states_matrix_0.setState( 5, 4, O );
            binary_states_matrix_0.setState( 6, 4, I );
            binary_states_matrix_0.setState( 0, 5, I );
            binary_states_matrix_0.setState( 1, 5, I );
            binary_states_matrix_0.setState( 2, 5, I );
            binary_states_matrix_0.setState( 3, 5, I );
            binary_states_matrix_0.setState( 4, 5, I );
            binary_states_matrix_0.setState( 5, 5, I );
            binary_states_matrix_0.setState( 6, 5, I );
            final String[] character_labels_0 = new String[ 6 ];
            character_labels_0[ 0 ] = "first";
            character_labels_0[ 1 ] = "second";
            character_labels_0[ 2 ] = "third";
            character_labels_0[ 3 ] = "forth";
            character_labels_0[ 4 ] = "fifth";
            character_labels_0[ 5 ] = "sixth";
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final Phylogeny phylogeny_0 = factory.create( "(((A,B)1,C)2,D)3", new NHXParser() )[ 0 ];
            final DomainParsimonyCalculator dom_pars = DomainParsimonyCalculator.createInstance( phylogeny_0 );
            dom_pars.executeOnGivenBinaryStatesMatrix( binary_states_matrix_0, character_labels_0 );
            final CharacterStateMatrix<GainLossStates> gl_matrix_0 = dom_pars.getGainLossMatrix();
            // final StringWriter sw = new StringWriter();
            //  gl_matrix_0.toWriter( sw );
            // System.out.println( sw.toString() );
            if ( dom_pars.getCost() != 13 ) {
                return false;
            }
            if ( dom_pars.getTotalGains() != 5 ) {
                return false;
            }
            if ( dom_pars.getTotalLosses() != 8 ) {
                return false;
            }
            if ( dom_pars.getTotalUnchanged() != 29 ) {
                return false;
            }
            if ( gl_matrix_0.getState( "A", 1 ) != GainLossStates.GAIN ) {
                return false;
            }
            if ( gl_matrix_0.getState( "A", 4 ) != GainLossStates.UNCHANGED_PRESENT ) {
                return false;
            }
            if ( gl_matrix_0.getState( "B", 4 ) != GainLossStates.LOSS ) {
                return false;
            }
            if ( gl_matrix_0.getState( "C", 4 ) != GainLossStates.GAIN ) {
                return false;
            }
            if ( gl_matrix_0.getState( "D", 4 ) != GainLossStates.LOSS ) {
                return false;
            }
            if ( gl_matrix_0.getState( "1", 4 ) != GainLossStates.GAIN ) {
                return false;
            }
            if ( gl_matrix_0.getState( "2", 4 ) != GainLossStates.LOSS ) {
                return false;
            }
            if ( gl_matrix_0.getState( "3", 4 ) != GainLossStates.UNCHANGED_PRESENT ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testCombinableDomains() {
        try {
            final Domain key0 = new BasicDomain( "key0", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain a = new BasicDomain( "a", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain b = new BasicDomain( "b", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain c = new BasicDomain( "c", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final CombinableDomains cd0 = new BasicCombinableDomains( key0.getDomainId(), new BasicSpecies( "eel" ) );
            cd0.addCombinableDomain( a.getDomainId() );
            cd0.addCombinableDomain( b.getDomainId() );
            cd0.addCombinableDomain( b.getDomainId() );
            cd0.addCombinableDomain( c.getDomainId() );
            cd0.addCombinableDomain( c.getDomainId() );
            cd0.addCombinableDomain( c.getDomainId() );
            if ( cd0.getNumberOfCombinableDomains() != 3 ) {
                return false;
            }
            if ( cd0.getNumberOfProteinsExhibitingCombination( a.getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd0.getNumberOfProteinsExhibitingCombination( b.getDomainId() ) != 2 ) {
                return false;
            }
            if ( cd0.getNumberOfProteinsExhibitingCombination( c.getDomainId() ) != 3 ) {
                return false;
            }
            if ( cd0.getNumberOfProteinsExhibitingCombination( key0.getDomainId() ) != 0 ) {
                return false;
            }
            if ( cd0.getAllDomains().size() != 4 ) {
                return false;
            }
            if ( !cd0.getAllDomains().contains( a.getDomainId() ) ) {
                return false;
            }
            if ( !cd0.getAllDomains().contains( b.getDomainId() ) ) {
                return false;
            }
            if ( !cd0.getAllDomains().contains( c.getDomainId() ) ) {
                return false;
            }
            if ( !cd0.getAllDomains().contains( key0.getDomainId() ) ) {
                return false;
            }
            if ( cd0.toBinaryDomainCombinations().size() != 3 ) {
                return false;
            }
            final BasicBinaryDomainCombination s0 = new BasicBinaryDomainCombination( "key0", "a" );
            final BasicBinaryDomainCombination s1 = new BasicBinaryDomainCombination( "b", "key0" );
            final BasicBinaryDomainCombination s2 = new BasicBinaryDomainCombination( "key0", "c" );
            final BasicBinaryDomainCombination s3 = new BasicBinaryDomainCombination( "key0", "cc" );
            final BasicBinaryDomainCombination s4 = new BasicBinaryDomainCombination( "c", "key0" );
            if ( !cd0.toBinaryDomainCombinations().contains( s0 ) ) {
                return false;
            }
            if ( !cd0.toBinaryDomainCombinations().contains( s1 ) ) {
                return false;
            }
            if ( !cd0.toBinaryDomainCombinations().contains( s2 ) ) {
                return false;
            }
            if ( cd0.toBinaryDomainCombinations().contains( s3 ) ) {
                return false;
            }
            if ( !cd0.toBinaryDomainCombinations().contains( s4 ) ) {
                return false;
            }
            final Domain key1 = new BasicDomain( "key1", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain a1 = new BasicDomain( "a1", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain b1 = new BasicDomain( "b1", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain c1 = new BasicDomain( "c1", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final CombinableDomains cd1 = new BasicCombinableDomains( key1.getDomainId(), new BasicSpecies( "eel" ) );
            cd1.addCombinableDomain( a1.getDomainId() );
            cd1.addCombinableDomain( b1.getDomainId() );
            cd1.addCombinableDomain( c1.getDomainId() );
            cd1.addCombinableDomain( key1.getDomainId() );
            if ( cd1.getNumberOfCombinableDomains() != 4 ) {
                return false;
            }
            if ( cd1.getNumberOfProteinsExhibitingCombination( a1.getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd1.getNumberOfProteinsExhibitingCombination( b1.getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd1.getNumberOfProteinsExhibitingCombination( c1.getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd1.getNumberOfProteinsExhibitingCombination( key1.getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd1.getAllDomains().size() != 4 ) {
                return false;
            }
            if ( cd1.toBinaryDomainCombinations().size() != 4 ) {
                return false;
            }
            final BasicBinaryDomainCombination kk = new BasicBinaryDomainCombination( "key1", "key1" );
            if ( !cd1.toBinaryDomainCombinations().contains( kk ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testCombinationsBasedPairwiseSimilarityCalculator() {
        try {
            final Domain a = new BasicDomain( "A", 1, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain b = new BasicDomain( "B", 2, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain c = new BasicDomain( "C", 3, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain one_key = new BasicDomain( "bcl2", 4, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain two_key = new BasicDomain( "bcl2", 5, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final CombinableDomains one = new BasicCombinableDomains( one_key.getDomainId(), new BasicSpecies( "mouse" ) );
            final CombinableDomains two = new BasicCombinableDomains( two_key.getDomainId(),
                                                                      new BasicSpecies( "rabbit" ) );
            one.addCombinableDomain( a.getDomainId() );
            one.addCombinableDomain( a.getDomainId() );
            two.addCombinableDomain( new BasicDomain( "A", 1, 5, ( short ) 1, ( short ) 4, 0.1, -12 ).getDomainId() );
            two.addCombinableDomain( b.getDomainId() );
            two.addCombinableDomain( c.getDomainId() );
            final PairwiseDomainSimilarityCalculator calc = new CombinationsBasedPairwiseDomainSimilarityCalculator();
            final PairwiseDomainSimilarity s1 = calc.calculateSimilarity( one, two );
            if ( !TestSurfacing.isEqual( s1.getSimilarityScore(), 1.0 / ( 1 + 2 ) ) ) {
                return false;
            }
            if ( s1.getDifferenceInCounts() != ( 1 - 3 ) ) {
                return false;
            }
            if ( ( ( CombinationsBasedPairwiseDomainSimilarity ) s1 ).getNumberOfDifferentDomains() != 2 ) {
                return false;
            }
            one.addCombinableDomain( b.getDomainId() );
            one.addCombinableDomain( c.getDomainId() );
            final PairwiseDomainSimilarity s2 = calc.calculateSimilarity( one, two );
            if ( !TestSurfacing.isEqual( s2.getSimilarityScore(), 3.0 / ( 0 + 3 ) ) ) {
                return false;
            }
            if ( s2.getDifferenceInCounts() != 0 ) {
                return false;
            }
            if ( ( ( CombinationsBasedPairwiseDomainSimilarity ) s2 ).getNumberOfDifferentDomains() != 0 ) {
                return false;
            }
            final Domain d = new BasicDomain( "D", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain e = new BasicDomain( "E", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain f = new BasicDomain( "F", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            one.addCombinableDomain( d.getDomainId() );
            one.addCombinableDomain( d.getDomainId() );
            one.addCombinableDomain( e.getDomainId() );
            one.addCombinableDomain( f.getDomainId() );
            final PairwiseDomainSimilarity s3 = calc.calculateSimilarity( one, two );
            if ( !TestSurfacing.isEqual( s3.getSimilarityScore(), 3.0 / ( 3 + 3 ) ) ) {
                return false;
            }
            if ( s3.getDifferenceInCounts() != ( 6 - 3 ) ) {
                return false;
            }
            if ( ( ( CombinationsBasedPairwiseDomainSimilarity ) s3 ).getNumberOfDifferentDomains() != 3 ) {
                return false;
            }
            final Domain aaa = new BasicDomain( "aaa", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain bbb = new BasicDomain( "bbb", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain three_key = new BasicDomain( "bcl2", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain four_key = new BasicDomain( "bcl2", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final CombinableDomains three = new BasicCombinableDomains( three_key.getDomainId(),
                                                                        new BasicSpecies( "mouse" ) );
            final CombinableDomains four = new BasicCombinableDomains( four_key.getDomainId(),
                                                                       new BasicSpecies( "rabbit" ) );
            three.addCombinableDomain( aaa.getDomainId() );
            four.addCombinableDomain( bbb.getDomainId() );
            final PairwiseDomainSimilarityCalculator calc2 = new CombinationsBasedPairwiseDomainSimilarityCalculator();
            final PairwiseDomainSimilarity s4 = calc2.calculateSimilarity( three, four );
            if ( !TestSurfacing.isEqual( s4.getSimilarityScore(), 0.0 / ( 0 + 2 ) ) ) {
                return false;
            }
            final Domain aaa2 = new BasicDomain( "aaa", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            four.addCombinableDomain( aaa2.getDomainId() );
            final PairwiseDomainSimilarity s5 = calc.calculateSimilarity( three, four );
            if ( !TestSurfacing.isEqual( s5.getSimilarityScore(), 1.0 / ( 1 + 1 ) ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testCopyNumberBasedPairwiseSimilarityCalculator() {
        try {
            final Domain one_key = new BasicDomain( "bcl2", 4, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain two_key = new BasicDomain( "bcl2", 5, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final CombinableDomains one = new BasicCombinableDomains( one_key.getDomainId(), new BasicSpecies( "mouse" ) );
            final CombinableDomains two = new BasicCombinableDomains( two_key.getDomainId(),
                                                                      new BasicSpecies( "rabbit" ) );
            one.setKeyDomainCount( 2 );
            two.setKeyDomainCount( 3 );
            final PairwiseDomainSimilarityCalculator calc = new DomainCountsBasedPairwiseSimilarityCalculator();
            PairwiseDomainSimilarity s1 = calc.calculateSimilarity( one, two );
            if ( !TestSurfacing.isEqual( s1.getSimilarityScore(), 1.0 - ( 3 - 2.0 ) / ( 2 + 3 ) ) ) {
                return false;
            }
            if ( s1.getDifferenceInCounts() != ( 2 - 3 ) ) {
                return false;
            }
            one.setKeyDomainCount( 1 );
            two.setKeyDomainCount( 1 );
            s1 = calc.calculateSimilarity( one, two );
            if ( !TestSurfacing.isEqual( s1.getSimilarityScore(), 1.0 ) ) {
                return false;
            }
            if ( s1.getDifferenceInCounts() != ( 1 - 1 ) ) {
                return false;
            }
            one.setKeyDomainCount( 1 );
            two.setKeyDomainCount( 1000 );
            s1 = calc.calculateSimilarity( one, two );
            if ( !TestSurfacing.isEqual( s1.getSimilarityScore(), 1.0 - 999.0 / 1001 ) ) {
                return false;
            }
            if ( s1.getDifferenceInCounts() != ( 1 - 1000 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDirectedCombinableDomains() {
        try {
            final Domain key0 = new BasicDomain( "key0", 10, 20, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain a = new BasicDomain( "a", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain b = new BasicDomain( "b", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain c = new BasicDomain( "c", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final CombinableDomains cd0 = new DirectedCombinableDomains( key0.getDomainId(), new BasicSpecies( "eel" ) );
            cd0.addCombinableDomain( a.getDomainId() );
            cd0.addCombinableDomain( b.getDomainId() );
            cd0.addCombinableDomain( b.getDomainId() );
            cd0.addCombinableDomain( c.getDomainId() );
            cd0.addCombinableDomain( c.getDomainId() );
            cd0.addCombinableDomain( c.getDomainId() );
            if ( cd0.getNumberOfCombinableDomains() != 3 ) {
                return false;
            }
            if ( cd0.getNumberOfProteinsExhibitingCombination( a.getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd0.getNumberOfProteinsExhibitingCombination( b.getDomainId() ) != 2 ) {
                return false;
            }
            if ( cd0.getNumberOfProteinsExhibitingCombination( c.getDomainId() ) != 3 ) {
                return false;
            }
            if ( cd0.getNumberOfProteinsExhibitingCombination( key0.getDomainId() ) != 0 ) {
                return false;
            }
            if ( cd0.getAllDomains().size() != 4 ) {
                return false;
            }
            if ( !cd0.getAllDomains().contains( a.getDomainId() ) ) {
                return false;
            }
            if ( !cd0.getAllDomains().contains( b.getDomainId() ) ) {
                return false;
            }
            if ( !cd0.getAllDomains().contains( c.getDomainId() ) ) {
                return false;
            }
            if ( !cd0.getAllDomains().contains( key0.getDomainId() ) ) {
                return false;
            }
            if ( cd0.toBinaryDomainCombinations().size() != 3 ) {
                return false;
            }
            final BinaryDomainCombination s0 = new DirectedBinaryDomainCombination( "key0", "a" );
            final BinaryDomainCombination s1 = new DirectedBinaryDomainCombination( "b", "key0" );
            final BinaryDomainCombination s2 = new DirectedBinaryDomainCombination( "key0", "c" );
            final BinaryDomainCombination s3 = new DirectedBinaryDomainCombination( "key0", "cc" );
            final BinaryDomainCombination s4 = new DirectedBinaryDomainCombination( "a", "b" );
            final BinaryDomainCombination s5 = new DirectedBinaryDomainCombination( "b", "a" );
            final BinaryDomainCombination s6 = new DirectedBinaryDomainCombination( "key0", "b" );
            final BinaryDomainCombination s7 = new DirectedBinaryDomainCombination( "a", "key0" );
            final BinaryDomainCombination s8 = new DirectedBinaryDomainCombination( "c", "key0" );
            if ( !cd0.toBinaryDomainCombinations().contains( s0 ) ) {
                return false;
            }
            if ( cd0.toBinaryDomainCombinations().contains( s1 ) ) {
                return false;
            }
            if ( !cd0.toBinaryDomainCombinations().contains( s2 ) ) {
                return false;
            }
            if ( cd0.toBinaryDomainCombinations().contains( s3 ) ) {
                return false;
            }
            if ( cd0.toBinaryDomainCombinations().contains( s4 ) ) {
                return false;
            }
            if ( cd0.toBinaryDomainCombinations().contains( s5 ) ) {
                return false;
            }
            if ( !cd0.toBinaryDomainCombinations().contains( s6 ) ) {
                return false;
            }
            if ( cd0.toBinaryDomainCombinations().contains( s7 ) ) {
                return false;
            }
            if ( cd0.toBinaryDomainCombinations().contains( s8 ) ) {
                return false;
            }
            final Domain key1 = new BasicDomain( "key1", 1, 2, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain a1 = new BasicDomain( "a1", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain b1 = new BasicDomain( "b1", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain c1 = new BasicDomain( "c1", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final CombinableDomains cd1 = new DirectedCombinableDomains( key1.getDomainId(), new BasicSpecies( "eel" ) );
            cd1.addCombinableDomain( a1.getDomainId() );
            cd1.addCombinableDomain( b1.getDomainId() );
            cd1.addCombinableDomain( c1.getDomainId() );
            cd1.addCombinableDomain( key1.getDomainId() );
            if ( cd1.getNumberOfCombinableDomains() != 4 ) {
                return false;
            }
            if ( cd1.getNumberOfProteinsExhibitingCombination( a1.getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd1.getNumberOfProteinsExhibitingCombination( b1.getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd1.getNumberOfProteinsExhibitingCombination( c1.getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd1.getNumberOfProteinsExhibitingCombination( key1.getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd1.getAllDomains().size() != 4 ) {
                return false;
            }
            if ( cd1.toBinaryDomainCombinations().size() != 4 ) {
                return false;
            }
            final BinaryDomainCombination kk = new DirectedBinaryDomainCombination( "key1", "key1" );
            if ( !cd1.toBinaryDomainCombinations().contains( kk ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDirectedness() {
        try {
            final BinaryStates X = BinaryStates.PRESENT;
            final BinaryStates O = BinaryStates.ABSENT;
            final GainLossStates G = GainLossStates.GAIN;
            final GainLossStates L = GainLossStates.LOSS;
            final GainLossStates A = GainLossStates.UNCHANGED_ABSENT;
            final GainLossStates P = GainLossStates.UNCHANGED_PRESENT;
            final Protein one_1 = new BasicProtein( "one", "1" );
            final Protein two_1 = new BasicProtein( "two", "1" );
            final Protein three_1 = new BasicProtein( "three", "1" );
            final Protein four_1 = new BasicProtein( "four", "1" );
            final Protein five_1 = new BasicProtein( "five", "1" );
            one_1.addProteinDomain( new BasicDomain( "B", 12, 14, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            one_1.addProteinDomain( new BasicDomain( "C", 13, 14, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            one_1.addProteinDomain( new BasicDomain( "A", 11, 12, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            one_1.addProteinDomain( new BasicDomain( "X", 100, 110, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            one_1.addProteinDomain( new BasicDomain( "Y", 200, 210, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            two_1.addProteinDomain( new BasicDomain( "A", 10, 20, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            two_1.addProteinDomain( new BasicDomain( "B", 30, 40, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            two_1.addProteinDomain( new BasicDomain( "Y", 1, 2, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            two_1.addProteinDomain( new BasicDomain( "X", 10, 11, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            three_1.addProteinDomain( new BasicDomain( "P", 10, 11, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            three_1.addProteinDomain( new BasicDomain( "M", 1, 2, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            three_1.addProteinDomain( new BasicDomain( "M", 5, 6, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            three_1.addProteinDomain( new BasicDomain( "N", 7, 8, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            three_1.addProteinDomain( new BasicDomain( "N", 3, 4, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            four_1.addProteinDomain( new BasicDomain( "XX", 10, 20, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            five_1.addProteinDomain( new BasicDomain( "YY", 30, 40, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            final List<Protein> list_1 = new ArrayList<Protein>();
            list_1.add( one_1 );
            list_1.add( two_1 );
            list_1.add( three_1 );
            list_1.add( four_1 );
            list_1.add( five_1 );
            final GenomeWideCombinableDomains gwcd_1 = BasicGenomeWideCombinableDomains
                    .createInstance( list_1, false, new BasicSpecies( "1" ), DomainCombinationType.DIRECTED );
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "A", "B" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "B", "A" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "A", "A" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "A", "C" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "C", "A" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "B", "C" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "C", "X" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "C", "Y" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "A", "X" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "A", "Y" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "Y", "A" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "X", "A" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "C", "B" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "X", "Y" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "Y", "X" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "A", "Y" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "A", "X" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "Y", "C" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "M", "N" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "N", "M" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "N", "P" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "M", "P" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "P", "N" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "P", "M" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "XX", "YY" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "YY", "XX" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations().contains( new DirectedBinaryDomainCombination( "B", "B" ) ) ) {
                return false;
            }
            //            final List<GenomeWideCombinableDomains> gwcd_list = new ArrayList<GenomeWideCombinableDomains>();
            //            gwcd_list.add( gwcd_1 );
            //            gwcd_list.add( gwcd_2 );
            //            final CharacterStateMatrix<BinaryStates> matrix_d = DomainParsimonyCalculator
            //                    .createMatrixOfDomainPresenceOrAbsence( gwcd_list );
            //            final CharacterStateMatrix<BinaryStates> matrix_bc = DomainParsimonyCalculator
            //                    .createMatrixOfBinaryDomainCombinationPresenceOrAbsence( gwcd_list );
            //            if ( matrix_d.getState( 0, 0 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 0, 0 ) != X ) {
            //                return false;
            //            }
            //        
            //
            //            final BasicCharacterStateMatrix<BinaryStates> dm = new BasicCharacterStateMatrix<BinaryStates>( new BinaryStates[][] {
            //                    { X, X, X, X, X, X }, { X, X, X, X, X, X } } );
            //            if ( !matrix_d.equals( dm ) ) {
            //                return false;
            //            }
            //            final BasicCharacterStateMatrix<BinaryStates> bcm = new BasicCharacterStateMatrix<BinaryStates>( new BinaryStates[][] {
            //                    { X, O, X, X, X, X, O, X, X, O, X, X }, { X, X, X, O, O, O, O, X, O, O, X, X } } );
            //            if ( !matrix_d.equals( dm ) ) {
            //                return false;
            //            }
            //``````````````````````````
            //            final List<GenomeWideCombinableDomains> gwcd_list = new ArrayList<GenomeWideCombinableDomains>();
            //            gwcd_list.add( one );
            //            gwcd_list.add( two );
            //            gwcd_list.add( three );
            //            gwcd_list.add( four );
            //            final CharacterStateMatrix<BinaryStates> matrix_d = DomainParsimony
            //                    .createMatrixOfDomainPresenceOrAbsence( gwcd_list );
            //            final CharacterStateMatrix<BinaryStates> matrix_bc = DomainParsimony
            //                    .createMatrixOfBinaryDomainCombinationPresenceOrAbsence( gwcd_list );
            //            //         System.out.println( "d:"  );
            //            //         System.out.println(matrix_d.toStringBuffer().toString()  );
            //            //         System.out.println( "bc:"  );
            //            //        System.out.println(matrix_bc.toStringBuffer().toString()  );
            //            // 1 a b c e f g h l m
            //            // 2 a b c e f g i n o
            //            // 3 a b d e f g j p q
            //            // 4 a b d p r
            //            if ( matrix_d.getState( 0, 0 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_d.getState( 0, 1 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_d.getState( 0, 2 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_d.getState( 0, 3 ) != O ) {
            //                return false;
            //            }
            //            if ( matrix_d.getState( 0, 4 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_d.getState( 0, 5 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_d.getState( 0, 6 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_d.getState( 0, 7 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_d.getState( 0, 8 ) != O ) {
            //                return false;
            //            }
            //            // 1 a-a a-b a-c e-f e-g e-h f-g f-h g-h l-m
            //            // 2 a-b a-c e-f e-g e-i f-g f-i g-i n-o
            //            // 3 a-b a-d e-f e-g e-j f-g f-j g-j p-q
            //            // 4 a-b a-d p-r
            //            if ( matrix_bc.getState( 0, 0 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 0, 1 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 0, 2 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 0, 3 ) != O ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 0, 4 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 1, 0 ) != O ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 1, 1 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 1, 2 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 1, 3 ) != O ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 1, 4 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 2, 0 ) != O ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 2, 1 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 2, 2 ) != O ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 2, 3 ) != X ) {
            //                return false;
            //            }
            //            if ( matrix_bc.getState( 2, 4 ) != X ) {
            //                return false;
            //            }
            //            final PhylogenyFactory factory0 = ParserBasedPhylogenyFactory.getInstance();
            //            final String p0_str = "((one,two)1-2,(three,four)3-4)root";
            //            final Phylogeny p0 = factory0.create( p0_str, new NHXParser() )[ 0 ];
            //            final DomainParsimony dp0 = DomainParsimony.createInstance( p0, gwcd_list );
            //            dp0.executeDolloParsimonyOnDomainPresence();
            //            final CharacterStateMatrix<GainLossStates> gl_matrix_d = dp0.getGainLossMatrix();
            //            final CharacterStateMatrix<BinaryStates> is_matrix_d = dp0.getInternalStatesMatrix();
            //            dp0.executeDolloParsimonyOnBinaryDomainCombintionPresence();
            //            final CharacterStateMatrix<GainLossStates> gl_matrix_bc = dp0.getGainLossMatrix();
            //            final CharacterStateMatrix<BinaryStates> is_matrix_bc = dp0.getInternalStatesMatrix();
            //            if ( is_matrix_d.getState( "root", "A" ) != X ) {
            //                return false;
            //            }
            //            if ( is_matrix_d.getState( "root", "B" ) != X ) {
            //                return false;
            //            }
            //            if ( is_matrix_d.getState( "root", "C" ) != O ) {
            //                return false;
            //            }
            //            if ( is_matrix_d.getState( "root", "D" ) != O ) {
            //                return false;
            //            }
            //            if ( is_matrix_d.getState( "root", "E" ) != X ) {
            //                return false;
            //            }
            //            if ( is_matrix_bc.getState( "root", "A=A" ) != O ) {
            //                return false;
            //            }
            //            if ( is_matrix_bc.getState( "root", "A=B" ) != X ) {
            //                return false;
            //            }
            //            if ( is_matrix_bc.getState( "root", "A=C" ) != O ) {
            //                return false;
            //            }
            //            if ( is_matrix_bc.getState( "root", "A=D" ) != O ) {
            //                return false;
            //            }
            //            if ( is_matrix_bc.getState( "root", "G=H" ) != O ) {
            //                return false;
            //            }
            //            if ( is_matrix_bc.getState( "1-2", "G=H" ) != O ) {
            //                return false;
            //            }
            //            if ( is_matrix_bc.getState( "root", "E=F" ) != X ) {
            //                return false;
            //            }
            //            if ( gl_matrix_bc.getState( "root", "E=F" ) != P ) {
            //                return false;
            //            }
            //            if ( gl_matrix_bc.getState( "root", "A=A" ) != A ) {
            //                return false;
            //            }
            //            if ( gl_matrix_bc.getState( "one", "A=A" ) != G ) {
            //                return false;
            //            }
            //            if ( gl_matrix_bc.getState( "root", "A=B" ) != P ) {
            //                return false;
            //            }
            //            if ( gl_matrix_bc.getState( "3-4", "A=D" ) != G ) {
            //                return false;
            //            }
            //            if ( gl_matrix_bc.getState( "four", "E=F" ) != L ) {
            //                return false;
            //            }
            //            if ( gl_matrix_d.getState( "3-4", "P" ) != G ) {
            //                return false;
            //            }
            //            final Protein ab_1 = new BasicProtein( "ab", "one" );
            //            ab_1.addProteinDomain( a );
            //            ab_1.addProteinDomain( b );
            //            final Protein ac_1 = new BasicProtein( "ac", "one" );
            //            ac_1.addProteinDomain( a );
            //            ac_1.addProteinDomain( c );
            //            final Protein de_1 = new BasicProtein( "de", "one" );
            //            de_1.addProteinDomain( d );
            //            de_1.addProteinDomain( e );
            //            final Protein ac_2 = new BasicProtein( "ac", "two" );
            //            ac_2.addProteinDomain( a );
            //            ac_2.addProteinDomain( c );
            //            final Protein ab_3 = new BasicProtein( "ab", "three" );
            //            ab_3.addProteinDomain( a );
            //            ab_3.addProteinDomain( b );
            //            final Protein de_4 = new BasicProtein( "de", "four" );
            //            de_4.addProteinDomain( d );
            //            de_4.addProteinDomain( e );
            //            final Protein ab_6 = new BasicProtein( "ab", "six" );
            //            ab_6.addProteinDomain( a );
            //            ab_6.addProteinDomain( b );
            //            final List<Protein> spec_one = new ArrayList<Protein>();
            //            final List<Protein> spec_two = new ArrayList<Protein>();
            //            final List<Protein> spec_three = new ArrayList<Protein>();
            //            final List<Protein> spec_four = new ArrayList<Protein>();
            //            final List<Protein> spec_five = new ArrayList<Protein>();
            //            final List<Protein> spec_six = new ArrayList<Protein>();
            //            final List<Protein> spec_seven = new ArrayList<Protein>();
            //            spec_one.add( ab_1 );
            //            spec_one.add( ac_1 );
            //            spec_one.add( de_1 );
            //            spec_two.add( ac_2 );
            //            spec_three.add( ab_3 );
            //            spec_four.add( de_4 );
            //            spec_six.add( ab_6 );
            //            final GenomeWideCombinableDomains one_gwcd = BasicGenomeWideCombinableDomains
            //                    .createInstance( spec_one, false, new BasicSpecies( "one" ), false );
            //            final GenomeWideCombinableDomains two_gwcd = BasicGenomeWideCombinableDomains
            //                    .createInstance( spec_two, false, new BasicSpecies( "two" ), false );
            //            final GenomeWideCombinableDomains three_gwcd = BasicGenomeWideCombinableDomains
            //                    .createInstance( spec_three, false, new BasicSpecies( "three" ), false );
            //            final GenomeWideCombinableDomains four_gwcd = BasicGenomeWideCombinableDomains
            //                    .createInstance( spec_four, false, new BasicSpecies( "four" ), false );
            //            final GenomeWideCombinableDomains five_gwcd = BasicGenomeWideCombinableDomains
            //                    .createInstance( spec_five, false, new BasicSpecies( "five" ), false );
            //            final GenomeWideCombinableDomains six_gwcd = BasicGenomeWideCombinableDomains
            //                    .createInstance( spec_six, false, new BasicSpecies( "six" ), false );
            //            final GenomeWideCombinableDomains seven_gwcd = BasicGenomeWideCombinableDomains
            //                    .createInstance( spec_seven, false, new BasicSpecies( "seven" ), false
            //                                    );
            //            final List<GenomeWideCombinableDomains> gwcd_list1 = new ArrayList<GenomeWideCombinableDomains>();
            //            gwcd_list1.add( one_gwcd );
            //            gwcd_list1.add( two_gwcd );
            //            gwcd_list1.add( three_gwcd );
            //            gwcd_list1.add( four_gwcd );
            //            gwcd_list1.add( five_gwcd );
            //            gwcd_list1.add( six_gwcd );
            //            gwcd_list1.add( seven_gwcd );
            //            final PhylogenyFactory factory1 = ParserBasedPhylogenyFactory.getInstance();
            //            final String p1_str = "(((((one,two)12,three)123,(four,five)45)12345,six)123456,seven)root";
            //            final Phylogeny p1 = factory1.create( p1_str, new NHXParser() )[ 0 ];
            //            final DomainParsimony dp1 = DomainParsimony.createInstance( p1, gwcd_list1 );
            //            dp1.executeDolloParsimonyOnDomainPresence();
            //            final CharacterStateMatrix<GainLossStates> gl_dollo_d = dp1.getGainLossMatrix();
            //            final CharacterStateMatrix<BinaryStates> i_dollo_d = dp1.getInternalStatesMatrix();
            //            if ( dp1.getCost() != 14 ) {
            //                return false;
            //            }
            //            if ( dp1.getTotalGains() != 5 ) {
            //                return false;
            //            }
            //            if ( dp1.getTotalLosses() != 9 ) {
            //                return false;
            //            }
            //            if ( dp1.getTotalUnchanged() != 51 ) {
            //                return false;
            //            }
            //            if ( dp1.getNetGainsOnNode( "45" ) != -2 ) {
            //                return false;
            //            }
            //            if ( dp1.getSumOfGainsOnNode( "45" ) != 0 ) {
            //                return false;
            //            }
            //            if ( dp1.getSumOfLossesOnNode( "45" ) != 2 ) {
            //                return false;
            //            }
            //            if ( dp1.getSumOfUnchangedOnNode( "45" ) != 3 ) {
            //                return false;
            //            }
            //            if ( dp1.getSumOfUnchangedPresentOnNode( "45" ) != 2 ) {
            //                return false;
            //            }
            //            if ( dp1.getSumOfUnchangedAbsentOnNode( "45" ) != 1 ) {
            //                return false;
            //            }
            //            if ( dp1.getUnitsGainedOnNode( "45" ).contains( "A" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsLostOnNode( "45" ).contains( "A" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsLostOnNode( "45" ).contains( "B" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsGainedOnNode( "12345" ).contains( "D" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsOnNode( "12" ).contains( "A" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsOnNode( "12" ).contains( "B" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsOnNode( "12" ).contains( "C" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsOnNode( "12" ).contains( "D" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsOnNode( "12" ).contains( "E" ) ) {
            //                return false;
            //            }
            //            if ( dp1.getNetGainsOnNode( "123456" ) != 2 ) {
            //                return false;
            //            }
            //            if ( dp1.getSumOfGainsOnNode( "123456" ) != 2 ) {
            //                return false;
            //            }
            //            dp1.executeDolloParsimonyOnBinaryDomainCombintionPresence();
            //            final CharacterStateMatrix<GainLossStates> gl_dollo_bc = dp1.getGainLossMatrix();
            //            final CharacterStateMatrix<BinaryStates> i_dollo_bc = dp1.getInternalStatesMatrix();
            //            if ( dp1.getCost() != 8 ) {
            //                return false;
            //            }
            //            if ( dp1.getTotalGains() != 3 ) {
            //                return false;
            //            }
            //            if ( dp1.getTotalLosses() != 5 ) {
            //                return false;
            //            }
            //            if ( dp1.getTotalUnchanged() != 31 ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsLostOnNode( "45" ).contains( "A=B" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsGainedOnNode( "12345" ).contains( "D=E" ) ) {
            //                return false;
            //            }
            //            dp1.executeFitchParsimonyOnDomainPresence();
            //            final CharacterStateMatrix<GainLossStates> gl_fitch_d = dp1.getGainLossMatrix();
            //            final CharacterStateMatrix<BinaryStates> i_fitch_d = dp1.getInternalStatesMatrix();
            //            if ( dp1.getCost() != 10 ) {
            //                return false;
            //            }
            //            if ( dp1.getTotalGains() != 7 ) {
            //                return false;
            //            }
            //            if ( dp1.getTotalLosses() != 3 ) {
            //                return false;
            //            }
            //            if ( dp1.getTotalUnchanged() != 55 ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsGainedOnNode( "four" ).contains( "E" ) ) {
            //                return false;
            //            }
            //            dp1.executeFitchParsimonyOnBinaryDomainCombintion();
            //            final CharacterStateMatrix<GainLossStates> gl_fitch_bc = dp1.getGainLossMatrix();
            //            final CharacterStateMatrix<BinaryStates> i_fitch_bc = dp1.getInternalStatesMatrix();
            //            if ( dp1.getCost() != 6 ) {
            //                return false;
            //            }
            //            if ( dp1.getTotalGains() != 4 ) {
            //                return false;
            //            }
            //            if ( dp1.getTotalLosses() != 2 ) {
            //                return false;
            //            }
            //            if ( dp1.getTotalUnchanged() != 33 ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsLostOnNode( "45" ).contains( "A=B" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsGainedOnNode( "four" ).contains( "D=E" ) ) {
            //                return false;
            //            }
            //            if ( dp1.getNetGainsOnNode( "two" ) != -1 ) {
            //                return false;
            //            }
            //            if ( dp1.getNetGainsOnNode( "123" ) != 0 ) {
            //                return false;
            //            }
            //            if ( dp1.getSumOfUnchangedPresentOnNode( "123" ) != 1 ) {
            //                return false;
            //            }
            //            if ( dp1.getSumOfUnchangedAbsentOnNode( "123" ) != 2 ) {
            //                return false;
            //            }
            //            if ( dp1.getSumOfUnchangedOnNode( "123" ) != 3 ) {
            //                return false;
            //            }
            //            if ( dp1.getSumOfUnchangedOnNode( "two" ) != 2 ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsUnchangedAbsentOnNode( "two" ).contains( "D=E" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsUnchangedPresentOnNode( "two" ).contains( "A=C" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsUnchangedAbsentOnNode( "123" ).contains( "A=C" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsUnchangedPresentOnNode( "123" ).contains( "A=B" ) ) {
            //                return false;
            //            }
            //            if ( !dp1.getUnitsUnchangedAbsentOnNode( "123" ).contains( "D=E" ) ) {
            //                return false;
            //            }
            //            CharacterStateMatrix<BinaryStates> bsm = null;
            //            CharacterStateMatrix<GainLossStates> glm = null;
            //            bsm = new BasicCharacterStateMatrix<BinaryStates>( new BinaryStates[][] { { X, X, X, X, X },
            //                    { X, X, O, X, X }, { O, O, O, X, X }, { X, X, O, X, X }, { X, X, O, O, O }, { O, O, O, O, O } } );
            //            if ( !bsm.equals( i_dollo_d ) ) {
            //                return false;
            //            }
            //            bsm = new BasicCharacterStateMatrix<BinaryStates>( new BinaryStates[][] { { X, X, X, O, O },
            //                    { X, X, O, O, O }, { O, O, O, O, O }, { X, X, O, O, O }, { X, X, O, O, O }, { O, O, O, O, O } } );
            //            if ( !bsm.equals( i_fitch_d ) ) {
            //                return false;
            //            }
            //            glm = new BasicCharacterStateMatrix<GainLossStates>( new GainLossStates[][] { { P, P, P, P, P },
            //                    { P, L, P, L, L }, { P, P, G, P, P }, { P, P, A, L, L }, { P, P, A, P, P }, { A, A, A, P, P },
            //                    { A, A, A, L, L }, { L, L, A, P, P }, { P, P, A, G, G }, { P, P, A, A, A }, { G, G, A, A, A },
            //                    { A, A, A, A, A }, { A, A, A, A, A } } );
            //            if ( !glm.equals( gl_dollo_d ) ) {
            //                return false;
            //            }
            //            glm = new BasicCharacterStateMatrix<GainLossStates>( new GainLossStates[][] { { P, P, P, G, G },
            //                    { P, L, P, A, A }, { P, P, G, A, A }, { P, P, A, A, A }, { P, P, A, A, A }, { A, A, A, G, G },
            //                    { A, A, A, A, A }, { L, L, A, A, A }, { P, P, A, A, A }, { P, P, A, A, A }, { G, G, A, A, A },
            //                    { A, A, A, A, A }, { A, A, A, A, A } } );
            //            if ( !glm.equals( gl_fitch_d ) ) {
            //                return false;
            //            }
            //            bsm = new BasicCharacterStateMatrix<BinaryStates>( new BinaryStates[][] { { X, X, X }, { X, O, X },
            //                    { O, O, X }, { X, O, X }, { X, O, O }, { O, O, O } } );
            //            if ( !bsm.equals( i_dollo_bc ) ) {
            //                return false;
            //            }
            //            bsm = new BasicCharacterStateMatrix<BinaryStates>( new BinaryStates[][] { { X, X, O }, { X, O, O },
            //                    { O, O, O }, { X, O, O }, { X, O, O }, { O, O, O } } );
            //            if ( !bsm.equals( i_fitch_bc ) ) {
            //                return false;
            //            }
            //            glm = new BasicCharacterStateMatrix<GainLossStates>( new GainLossStates[][] { { P, P, P }, { L, P, L },
            //                    { P, G, P }, { P, A, L }, { P, A, P }, { A, A, P }, { A, A, L }, { L, A, P }, { P, A, G },
            //                    { P, A, A }, { G, A, A }, { A, A, A }, { A, A, A } } );
            //            if ( !glm.equals( gl_dollo_bc ) ) {
            //                return false;
            //            }
            //            glm = new BasicCharacterStateMatrix<GainLossStates>( new GainLossStates[][] { { P, P, G }, { L, P, A },
            //                    { P, G, A }, { P, A, A }, { P, A, A }, { A, A, G }, { A, A, A }, { L, A, A }, { P, A, A },
            //                    { P, A, A }, { G, A, A }, { A, A, A }, { A, A, A } } );
            //            if ( !glm.equals( gl_fitch_bc ) ) {
            //                return false;
            //            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDirectednessAndAdjacency() {
        try {
            final Protein one_1 = new BasicProtein( "one", "1" );
            final Protein two_1 = new BasicProtein( "two", "1" );
            final Protein three_1 = new BasicProtein( "three", "1" );
            final Protein four_1 = new BasicProtein( "four", "1" );
            final Protein five_1 = new BasicProtein( "five", "1" );
            one_1.addProteinDomain( new BasicDomain( "B", 12, 14, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            one_1.addProteinDomain( new BasicDomain( "C", 13, 14, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            one_1.addProteinDomain( new BasicDomain( "A", 11, 12, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            one_1.addProteinDomain( new BasicDomain( "X", 100, 110, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            one_1.addProteinDomain( new BasicDomain( "Y", 200, 210, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            two_1.addProteinDomain( new BasicDomain( "A", 10, 20, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            two_1.addProteinDomain( new BasicDomain( "B", 30, 40, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            two_1.addProteinDomain( new BasicDomain( "Y", 1, 2, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            two_1.addProteinDomain( new BasicDomain( "X", 10, 11, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            three_1.addProteinDomain( new BasicDomain( "P", 10, 11, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            three_1.addProteinDomain( new BasicDomain( "M", 1, 2, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            three_1.addProteinDomain( new BasicDomain( "M", 5, 6, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            three_1.addProteinDomain( new BasicDomain( "N", 7, 8, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            three_1.addProteinDomain( new BasicDomain( "N", 3, 4, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            four_1.addProteinDomain( new BasicDomain( "XX", 10, 20, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            five_1.addProteinDomain( new BasicDomain( "YY", 30, 40, ( short ) 1, ( short ) 4, 0.1, -12 ) );
            final List<Protein> list_1 = new ArrayList<Protein>();
            list_1.add( one_1 );
            list_1.add( two_1 );
            list_1.add( three_1 );
            list_1.add( four_1 );
            list_1.add( five_1 );
            final GenomeWideCombinableDomains gwcd_1 = BasicGenomeWideCombinableDomains
                    .createInstance( list_1, false, new BasicSpecies( "1" ), DomainCombinationType.DIRECTED_ADJACTANT );
            if ( !gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "A", "B" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "B", "A" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "A", "A" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "A", "C" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "C", "A" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "B", "C" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "C", "X" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "C", "Y" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "X", "Y" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "A", "X" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "A", "Y" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "Y", "A" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "X", "A" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "C", "B" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "X", "Y" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "Y", "X" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "A", "Y" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "A", "X" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "Y", "C" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "M", "N" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "N", "M" ) ) ) {
                return false;
            }
            if ( !gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "N", "P" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "M", "P" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "P", "N" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "P", "M" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "XX", "YY" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "YY", "XX" ) ) ) {
                return false;
            }
            if ( gwcd_1.toBinaryDomainCombinations()
                    .contains( new AdjactantDirectedBinaryDomainCombination( "B", "B" ) ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDomainArchitectureBasedGenomeSimilarityCalculator() {
        try {
            final Domain a = new BasicDomain( "a", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain b = new BasicDomain( "b", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain c = new BasicDomain( "c", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain d = new BasicDomain( "d", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain e = new BasicDomain( "e", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain f = new BasicDomain( "f", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain g = new BasicDomain( "g", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain h = new BasicDomain( "h", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain i = new BasicDomain( "i", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain j = new BasicDomain( "j", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain k = new BasicDomain( "k", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain l = new BasicDomain( "l", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain m = new BasicDomain( "m", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain n = new BasicDomain( "n", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Protein eel_0 = new BasicProtein( "0", "eel" );
            final Protein eel_1 = new BasicProtein( "1", "eel" );
            final Protein eel_2 = new BasicProtein( "2", "eel" );
            final Protein eel_3 = new BasicProtein( "3", "eel" );
            final Protein eel_4 = new BasicProtein( "4", "eel" );
            final Protein eel_5 = new BasicProtein( "5", "eel" );
            final Protein eel_6 = new BasicProtein( "6", "eel" );
            final Protein rat_0 = new BasicProtein( "0", "rat" );
            final Protein rat_1 = new BasicProtein( "1", "rat" );
            final Protein rat_2 = new BasicProtein( "2", "rat" );
            final Protein rat_3 = new BasicProtein( "3", "rat" );
            final Protein rat_4 = new BasicProtein( "4", "rat" );
            final Protein rat_5 = new BasicProtein( "5", "rat" );
            final Protein rat_6 = new BasicProtein( "6", "rat" );
            final Protein rat_7 = new BasicProtein( "7", "rat" );
            eel_1.addProteinDomain( a );
            eel_2.addProteinDomain( a );
            eel_2.addProteinDomain( b );
            eel_3.addProteinDomain( a );
            eel_3.addProteinDomain( a );
            eel_3.addProteinDomain( b );
            eel_4.addProteinDomain( a );
            eel_4.addProteinDomain( b );
            eel_4.addProteinDomain( c );
            eel_4.addProteinDomain( d );
            eel_4.addProteinDomain( e );
            eel_5.addProteinDomain( e );
            eel_5.addProteinDomain( e );
            eel_5.addProteinDomain( f );
            eel_5.addProteinDomain( f );
            eel_5.addProteinDomain( f );
            eel_5.addProteinDomain( f );
            eel_6.addProteinDomain( g );
            eel_6.addProteinDomain( h );
            rat_1.addProteinDomain( a );
            rat_2.addProteinDomain( a );
            rat_2.addProteinDomain( b );
            rat_3.addProteinDomain( a );
            rat_3.addProteinDomain( a );
            rat_3.addProteinDomain( b );
            rat_4.addProteinDomain( a );
            rat_4.addProteinDomain( b );
            rat_4.addProteinDomain( c );
            rat_4.addProteinDomain( i );
            rat_4.addProteinDomain( l );
            rat_5.addProteinDomain( i );
            rat_5.addProteinDomain( f );
            rat_5.addProteinDomain( f );
            rat_6.addProteinDomain( j );
            rat_6.addProteinDomain( k );
            rat_7.addProteinDomain( m );
            rat_7.addProteinDomain( n );
            final List<Protein> protein_list_eel = new ArrayList<Protein>();
            protein_list_eel.add( eel_0 );
            protein_list_eel.add( eel_1 );
            protein_list_eel.add( eel_2 );
            protein_list_eel.add( eel_3 );
            protein_list_eel.add( eel_4 );
            protein_list_eel.add( eel_5 );
            protein_list_eel.add( eel_6 );
            final List<Protein> protein_list_rat = new ArrayList<Protein>();
            protein_list_rat.add( rat_0 );
            protein_list_rat.add( rat_1 );
            protein_list_rat.add( rat_2 );
            protein_list_rat.add( rat_3 );
            protein_list_rat.add( rat_4 );
            protein_list_rat.add( rat_5 );
            protein_list_rat.add( rat_6 );
            protein_list_rat.add( rat_7 );
            final GenomeWideCombinableDomains eel_not_ignore = BasicGenomeWideCombinableDomains
                    .createInstance( protein_list_eel, false, new BasicSpecies( "eel" ) );
            final GenomeWideCombinableDomains eel_ignore = BasicGenomeWideCombinableDomains
                    .createInstance( protein_list_eel, true, new BasicSpecies( "eel" ) );
            final GenomeWideCombinableDomains rat_not_ignore = BasicGenomeWideCombinableDomains
                    .createInstance( protein_list_rat, false, new BasicSpecies( "rat" ) );
            final GenomeWideCombinableDomains rat_ignore = BasicGenomeWideCombinableDomains
                    .createInstance( protein_list_rat, true, new BasicSpecies( "rat" ) );
            final DomainArchitectureBasedGenomeSimilarityCalculator calc_ni = new DomainArchitectureBasedGenomeSimilarityCalculator( eel_not_ignore,
                                                                                                                                     rat_not_ignore );
            final DomainArchitectureBasedGenomeSimilarityCalculator calc_i = new DomainArchitectureBasedGenomeSimilarityCalculator( eel_ignore,
                                                                                                                                    rat_ignore );
            if ( calc_ni.getAllDomains().size() != 14 ) {
                return false;
            }
            if ( calc_i.getAllDomains().size() != 14 ) {
                return false;
            }
            if ( calc_ni.getDomainsSpecificToGenome0().size() != 4 ) {
                return false;
            }
            if ( calc_i.getDomainsSpecificToGenome0().size() != 4 ) {
                return false;
            }
            if ( calc_ni.getDomainsSpecificToGenome1().size() != 6 ) {
                return false;
            }
            if ( calc_i.getDomainsSpecificToGenome1().size() != 6 ) {
                return false;
            }
            if ( calc_i.getSharedDomains().size() != 4 ) {
                return false;
            }
            if ( calc_ni.getSharedDomains().size() != 4 ) {
                return false;
            }
            if ( !calc_ni.getDomainsSpecificToGenome0().contains( d.getDomainId() ) ) {
                return false;
            }
            if ( !calc_ni.getDomainsSpecificToGenome0().contains( e.getDomainId() ) ) {
                return false;
            }
            if ( !calc_ni.getDomainsSpecificToGenome0().contains( g.getDomainId() ) ) {
                return false;
            }
            if ( !calc_ni.getDomainsSpecificToGenome0().contains( h.getDomainId() ) ) {
                return false;
            }
            if ( calc_ni.getDomainsSpecificToGenome0().contains( a.getDomainId() ) ) {
                return false;
            }
            if ( calc_ni.getDomainsSpecificToGenome0().contains( i.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getDomainsSpecificToGenome0().contains( d.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getDomainsSpecificToGenome0().contains( e.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getDomainsSpecificToGenome0().contains( g.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getDomainsSpecificToGenome0().contains( h.getDomainId() ) ) {
                return false;
            }
            if ( calc_i.getDomainsSpecificToGenome0().contains( a.getDomainId() ) ) {
                return false;
            }
            if ( calc_i.getDomainsSpecificToGenome0().contains( i.getDomainId() ) ) {
                return false;
            }
            if ( !calc_ni.getDomainsSpecificToGenome1().contains( i.getDomainId() ) ) {
                return false;
            }
            if ( !calc_ni.getDomainsSpecificToGenome1().contains( l.getDomainId() ) ) {
                return false;
            }
            if ( !calc_ni.getDomainsSpecificToGenome1().contains( j.getDomainId() ) ) {
                return false;
            }
            if ( !calc_ni.getDomainsSpecificToGenome1().contains( k.getDomainId() ) ) {
                return false;
            }
            if ( !calc_ni.getDomainsSpecificToGenome1().contains( m.getDomainId() ) ) {
                return false;
            }
            if ( !calc_ni.getDomainsSpecificToGenome1().contains( n.getDomainId() ) ) {
                return false;
            }
            if ( calc_ni.getDomainsSpecificToGenome1().contains( a.getDomainId() ) ) {
                return false;
            }
            if ( calc_ni.getDomainsSpecificToGenome1().contains( b.getDomainId() ) ) {
                return false;
            }
            if ( calc_ni.getDomainsSpecificToGenome1().contains( d.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getDomainsSpecificToGenome1().contains( i.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getDomainsSpecificToGenome1().contains( l.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getDomainsSpecificToGenome1().contains( j.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getDomainsSpecificToGenome1().contains( k.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getDomainsSpecificToGenome1().contains( m.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getDomainsSpecificToGenome1().contains( n.getDomainId() ) ) {
                return false;
            }
            if ( calc_i.getDomainsSpecificToGenome1().contains( a.getDomainId() ) ) {
                return false;
            }
            if ( calc_i.getDomainsSpecificToGenome1().contains( b.getDomainId() ) ) {
                return false;
            }
            if ( calc_i.getDomainsSpecificToGenome1().contains( d.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getSharedDomains().contains( a.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getSharedDomains().contains( b.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getSharedDomains().contains( c.getDomainId() ) ) {
                return false;
            }
            if ( !calc_i.getSharedDomains().contains( f.getDomainId() ) ) {
                return false;
            }
            final Set<DomainId> all = calc_ni.getAllDomains();
            if ( !all.contains( a.getDomainId() ) ) {
                return false;
            }
            if ( !all.contains( b.getDomainId() ) ) {
                return false;
            }
            if ( !all.contains( c.getDomainId() ) ) {
                return false;
            }
            if ( !all.contains( d.getDomainId() ) ) {
                return false;
            }
            if ( !all.contains( e.getDomainId() ) ) {
                return false;
            }
            if ( !all.contains( f.getDomainId() ) ) {
                return false;
            }
            if ( !all.contains( g.getDomainId() ) ) {
                return false;
            }
            if ( !all.contains( h.getDomainId() ) ) {
                return false;
            }
            if ( !all.contains( i.getDomainId() ) ) {
                return false;
            }
            if ( !all.contains( l.getDomainId() ) ) {
                return false;
            }
            if ( !all.contains( j.getDomainId() ) ) {
                return false;
            }
            if ( !all.contains( k.getDomainId() ) ) {
                return false;
            }
            if ( !all.contains( m.getDomainId() ) ) {
                return false;
            }
            if ( !all.contains( n.getDomainId() ) ) {
                return false;
            }
            final Set<BinaryDomainCombination> s_0_ni = calc_ni.getBinaryDomainCombinationsSpecificToGenome0();
            final Set<BinaryDomainCombination> s_0_i = calc_i.getBinaryDomainCombinationsSpecificToGenome0();
            final Set<BinaryDomainCombination> s_1_ni = calc_ni.getBinaryDomainCombinationsSpecificToGenome1();
            final Set<BinaryDomainCombination> s_1_i = calc_i.getBinaryDomainCombinationsSpecificToGenome1();
            final Set<BinaryDomainCombination> a_ni = calc_ni.getAllBinaryDomainCombinations();
            final Set<BinaryDomainCombination> a_i = calc_i.getAllBinaryDomainCombinations();
            final Set<BinaryDomainCombination> shared_ni = calc_ni.getSharedBinaryDomainCombinations();
            final Set<BinaryDomainCombination> shared_i = calc_i.getSharedBinaryDomainCombinations();
            if ( a_ni.size() != 25 ) {
                return false;
            }
            if ( a_i.size() != 22 ) {
                return false;
            }
            if ( s_0_ni.size() != 10 ) {
                return false;
            }
            if ( s_0_i.size() != 9 ) {
                return false;
            }
            if ( s_1_ni.size() != 10 ) {
                return false;
            }
            if ( s_1_i.size() != 10 ) {
                return false;
            }
            if ( shared_ni.size() != 5 ) {
                return false;
            }
            if ( shared_i.size() != 3 ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "a", "a" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "b", "a" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "a", "c" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "a", "d" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "a", "e" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "b", "c" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "b", "d" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "b", "e" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "c", "d" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "c", "e" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "d", "e" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "e", "f" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "g", "h" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "f", "f" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "e", "e" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "a", "i" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "a", "l" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "b", "i" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "b", "l" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "c", "i" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "c", "l" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "i", "l" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "i", "f" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "m", "n" ) ) ) {
                return false;
            }
            if ( !a_ni.contains( new BasicBinaryDomainCombination( "j", "k" ) ) ) {
                return false;
            }
            if ( a_ni.contains( new BasicBinaryDomainCombination( "a", "g" ) ) ) {
                return false;
            }
            if ( a_ni.contains( new BasicBinaryDomainCombination( "a", "m" ) ) ) {
                return false;
            }
            if ( a_i.contains( new BasicBinaryDomainCombination( "a", "a" ) ) ) {
                return false;
            }
            if ( a_i.contains( new BasicBinaryDomainCombination( "f", "f" ) ) ) {
                return false;
            }
            if ( a_i.contains( new BasicBinaryDomainCombination( "e", "e" ) ) ) {
                return false;
            }
            if ( !shared_ni.contains( new BasicBinaryDomainCombination( "a", "a" ) ) ) {
                return false;
            }
            if ( !shared_ni.contains( new BasicBinaryDomainCombination( "a", "b" ) ) ) {
                return false;
            }
            if ( !shared_ni.contains( new BasicBinaryDomainCombination( "a", "c" ) ) ) {
                return false;
            }
            if ( !shared_ni.contains( new BasicBinaryDomainCombination( "b", "c" ) ) ) {
                return false;
            }
            if ( !shared_ni.contains( new BasicBinaryDomainCombination( "f", "f" ) ) ) {
                return false;
            }
            if ( shared_ni.contains( new BasicBinaryDomainCombination( "m", "n" ) ) ) {
                return false;
            }
            if ( shared_i.contains( new BasicBinaryDomainCombination( "a", "a" ) ) ) {
                return false;
            }
            if ( !shared_i.contains( new BasicBinaryDomainCombination( "a", "b" ) ) ) {
                return false;
            }
            if ( !shared_i.contains( new BasicBinaryDomainCombination( "a", "c" ) ) ) {
                return false;
            }
            if ( !shared_i.contains( new BasicBinaryDomainCombination( "b", "c" ) ) ) {
                return false;
            }
            if ( shared_i.contains( new BasicBinaryDomainCombination( "f", "f" ) ) ) {
                return false;
            }
            if ( shared_i.contains( new BasicBinaryDomainCombination( "m", "n" ) ) ) {
                return false;
            }
            if ( !s_0_ni.contains( new BasicBinaryDomainCombination( "a", "d" ) ) ) {
                return false;
            }
            if ( !s_0_ni.contains( new BasicBinaryDomainCombination( "a", "e" ) ) ) {
                return false;
            }
            if ( !s_0_ni.contains( new BasicBinaryDomainCombination( "b", "d" ) ) ) {
                return false;
            }
            if ( !s_0_ni.contains( new BasicBinaryDomainCombination( "b", "e" ) ) ) {
                return false;
            }
            if ( !s_0_ni.contains( new BasicBinaryDomainCombination( "c", "d" ) ) ) {
                return false;
            }
            if ( !s_0_ni.contains( new BasicBinaryDomainCombination( "c", "e" ) ) ) {
                return false;
            }
            if ( !s_0_ni.contains( new BasicBinaryDomainCombination( "d", "e" ) ) ) {
                return false;
            }
            if ( !s_0_ni.contains( new BasicBinaryDomainCombination( "e", "f" ) ) ) {
                return false;
            }
            if ( !s_0_ni.contains( new BasicBinaryDomainCombination( "g", "h" ) ) ) {
                return false;
            }
            if ( !s_0_ni.contains( new BasicBinaryDomainCombination( "e", "e" ) ) ) {
                return false;
            }
            if ( !s_0_i.contains( new BasicBinaryDomainCombination( "a", "d" ) ) ) {
                return false;
            }
            if ( !s_0_i.contains( new BasicBinaryDomainCombination( "a", "e" ) ) ) {
                return false;
            }
            if ( !s_0_i.contains( new BasicBinaryDomainCombination( "b", "d" ) ) ) {
                return false;
            }
            if ( !s_0_i.contains( new BasicBinaryDomainCombination( "b", "e" ) ) ) {
                return false;
            }
            if ( !s_0_i.contains( new BasicBinaryDomainCombination( "c", "d" ) ) ) {
                return false;
            }
            if ( !s_0_i.contains( new BasicBinaryDomainCombination( "c", "e" ) ) ) {
                return false;
            }
            if ( !s_0_i.contains( new BasicBinaryDomainCombination( "d", "e" ) ) ) {
                return false;
            }
            if ( !s_0_i.contains( new BasicBinaryDomainCombination( "e", "f" ) ) ) {
                return false;
            }
            if ( !s_0_i.contains( new BasicBinaryDomainCombination( "g", "h" ) ) ) {
                return false;
            }
            if ( s_0_i.contains( new BasicBinaryDomainCombination( "e", "e" ) ) ) {
                return false;
            }
            if ( !s_1_ni.contains( new BasicBinaryDomainCombination( "a", "i" ) ) ) {
                return false;
            }
            if ( !s_1_ni.contains( new BasicBinaryDomainCombination( "a", "l" ) ) ) {
                return false;
            }
            if ( !s_1_ni.contains( new BasicBinaryDomainCombination( "b", "i" ) ) ) {
                return false;
            }
            if ( !s_1_ni.contains( new BasicBinaryDomainCombination( "b", "l" ) ) ) {
                return false;
            }
            if ( !s_1_ni.contains( new BasicBinaryDomainCombination( "c", "i" ) ) ) {
                return false;
            }
            if ( !s_1_ni.contains( new BasicBinaryDomainCombination( "c", "l" ) ) ) {
                return false;
            }
            if ( !s_1_ni.contains( new BasicBinaryDomainCombination( "l", "i" ) ) ) {
                return false;
            }
            if ( !s_1_ni.contains( new BasicBinaryDomainCombination( "i", "f" ) ) ) {
                return false;
            }
            if ( !s_1_ni.contains( new BasicBinaryDomainCombination( "m", "n" ) ) ) {
                return false;
            }
            if ( !s_1_ni.contains( new BasicBinaryDomainCombination( "j", "k" ) ) ) {
                return false;
            }
            if ( s_1_ni.contains( new BasicBinaryDomainCombination( "a", "b" ) ) ) {
                return false;
            }
            if ( !s_1_i.contains( new BasicBinaryDomainCombination( "a", "i" ) ) ) {
                return false;
            }
            if ( !s_1_i.contains( new BasicBinaryDomainCombination( "a", "l" ) ) ) {
                return false;
            }
            if ( !s_1_i.contains( new BasicBinaryDomainCombination( "b", "i" ) ) ) {
                return false;
            }
            if ( !s_1_i.contains( new BasicBinaryDomainCombination( "b", "l" ) ) ) {
                return false;
            }
            if ( !s_1_i.contains( new BasicBinaryDomainCombination( "c", "i" ) ) ) {
                return false;
            }
            if ( !s_1_i.contains( new BasicBinaryDomainCombination( "c", "l" ) ) ) {
                return false;
            }
            if ( !s_1_i.contains( new BasicBinaryDomainCombination( "l", "i" ) ) ) {
                return false;
            }
            if ( !s_1_i.contains( new BasicBinaryDomainCombination( "i", "f" ) ) ) {
                return false;
            }
            if ( !s_1_i.contains( new BasicBinaryDomainCombination( "m", "n" ) ) ) {
                return false;
            }
            if ( !s_1_i.contains( new BasicBinaryDomainCombination( "j", "k" ) ) ) {
                return false;
            }
            if ( s_1_i.contains( new BasicBinaryDomainCombination( "a", "b" ) ) ) {
                return false;
            }
            if ( !isEqual( calc_ni.calculateSharedBinaryDomainCombinationBasedGenomeSimilarityScore(),
                           1.0 - ( 25.0 - 5.0 ) / 25.0 ) ) {
                return false;
            }
            if ( !isEqual( calc_i.calculateSharedBinaryDomainCombinationBasedGenomeSimilarityScore(),
                           1.0 - ( 22.0 - 3.0 ) / 22.0 ) ) {
                return false;
            }
            if ( !isEqual( calc_ni.calculateSharedDomainsBasedGenomeSimilarityScore(), 1.0 - ( 14.0 - 4.0 ) / 14.0 ) ) {
                return false;
            }
            if ( !isEqual( calc_i.calculateSharedDomainsBasedGenomeSimilarityScore(), 1.0 - ( 14.0 - 4.0 ) / 14.0 ) ) {
                return false;
            }
            final Domain u = new BasicDomain( "u", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain v = new BasicDomain( "v", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain w = new BasicDomain( "w", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain x = new BasicDomain( "x", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain y = new BasicDomain( "y", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain z = new BasicDomain( "z", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Protein a_0 = new BasicProtein( "0", "a" );
            final Protein a_1 = new BasicProtein( "1", "a" );
            final Protein a_2 = new BasicProtein( "2", "a" );
            final Protein b_0 = new BasicProtein( "0", "b" );
            final Protein b_1 = new BasicProtein( "1", "b" );
            a_0.addProteinDomain( u );
            a_0.addProteinDomain( v );
            a_0.addProteinDomain( w );
            a_1.addProteinDomain( w );
            a_1.addProteinDomain( x );
            a_2.addProteinDomain( y );
            a_2.addProteinDomain( z );
            b_0.addProteinDomain( u );
            b_0.addProteinDomain( w );
            b_1.addProteinDomain( y );
            b_1.addProteinDomain( z );
            final List<Protein> protein_list_a = new ArrayList<Protein>();
            protein_list_a.add( a_0 );
            protein_list_a.add( a_1 );
            protein_list_a.add( a_2 );
            final List<Protein> protein_list_b = new ArrayList<Protein>();
            protein_list_b.add( b_0 );
            protein_list_b.add( b_1 );
            final GenomeWideCombinableDomains ca = BasicGenomeWideCombinableDomains
                    .createInstance( protein_list_a, false, new BasicSpecies( "a" ) );
            final GenomeWideCombinableDomains cb = BasicGenomeWideCombinableDomains
                    .createInstance( protein_list_b, true, new BasicSpecies( "b" ) );
            final DomainArchitectureBasedGenomeSimilarityCalculator calc_u = new DomainArchitectureBasedGenomeSimilarityCalculator( ca,
                                                                                                                                    cb );
            calc_u.setAllowDomainsToBeIgnored( true );
            if ( calc_u.getAllDomains().size() != 6 ) {
                return false;
            }
            if ( calc_u.getDomainsSpecificToGenome0().size() != 2 ) {
                return false;
            }
            if ( calc_u.getDomainsSpecificToGenome1().size() != 0 ) {
                return false;
            }
            if ( !calc_u.getDomainsSpecificToGenome0().contains( v.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getDomainsSpecificToGenome0().contains( x.getDomainId() ) ) {
                return false;
            }
            if ( calc_u.getSharedDomains().size() != 4 ) {
                return false;
            }
            if ( !calc_u.getSharedDomains().contains( u.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getSharedDomains().contains( w.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getSharedDomains().contains( y.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getSharedDomains().contains( z.getDomainId() ) ) {
                return false;
            }
            if ( calc_u.getAllDomains().size() != 6 ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( u.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( w.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( y.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( z.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( v.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( x.getDomainId() ) ) {
                return false;
            }
            if ( calc_u.getBinaryDomainCombinationsSpecificToGenome0().size() != 3 ) {
                return false;
            }
            if ( calc_u.getBinaryDomainCombinationsSpecificToGenome1().size() != 0 ) {
                return false;
            }
            if ( calc_u.getSharedBinaryDomainCombinations().size() != 2 ) {
                return false;
            }
            if ( calc_u.getAllBinaryDomainCombinations().size() != 5 ) {
                return false;
            }
            if ( !calc_u.getBinaryDomainCombinationsSpecificToGenome0()
                    .contains( new BasicBinaryDomainCombination( "v", "u" ) ) ) {
                return false;
            }
            if ( !calc_u.getBinaryDomainCombinationsSpecificToGenome0()
                    .contains( new BasicBinaryDomainCombination( "w", "v" ) ) ) {
                return false;
            }
            if ( !calc_u.getBinaryDomainCombinationsSpecificToGenome0()
                    .contains( new BasicBinaryDomainCombination( "w", "x" ) ) ) {
                return false;
            }
            if ( !calc_u.getSharedBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "w", "u" ) ) ) {
                return false;
            }
            if ( !calc_u.getSharedBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "z", "y" ) ) ) {
                return false;
            }
            if ( !calc_u.getAllBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "v", "u" ) ) ) {
                return false;
            }
            if ( !calc_u.getAllBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "w", "v" ) ) ) {
                return false;
            }
            if ( !calc_u.getAllBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "w", "x" ) ) ) {
                return false;
            }
            if ( !calc_u.getAllBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "w", "u" ) ) ) {
                return false;
            }
            if ( !calc_u.getAllBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "z", "y" ) ) ) {
                return false;
            }
            calc_u.setAllowDomainsToBeIgnored( true );
            calc_u.addDomainIdToIgnore( u.getDomainId() );
            calc_u.addDomainIdToIgnore( new DomainId( "other" ) );
            calc_u.addDomainIdToIgnore( new DomainId( "other_too" ) );
            if ( calc_u.getAllDomains().size() != 5 ) {
                return false;
            }
            if ( calc_u.getDomainsSpecificToGenome0().size() != 2 ) {
                return false;
            }
            if ( calc_u.getDomainsSpecificToGenome1().size() != 0 ) {
                return false;
            }
            if ( !calc_u.getDomainsSpecificToGenome0().contains( v.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getDomainsSpecificToGenome0().contains( x.getDomainId() ) ) {
                return false;
            }
            if ( calc_u.getSharedDomains().size() != 3 ) {
                return false;
            }
            if ( calc_u.getSharedDomains().contains( u.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getSharedDomains().contains( w.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getSharedDomains().contains( y.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getSharedDomains().contains( z.getDomainId() ) ) {
                return false;
            }
            if ( calc_u.getAllDomains().size() != 5 ) {
                return false;
            }
            if ( calc_u.getAllDomains().contains( u.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( w.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( y.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( z.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( v.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( x.getDomainId() ) ) {
                return false;
            }
            if ( calc_u.getBinaryDomainCombinationsSpecificToGenome0().size() != 2 ) {
                return false;
            }
            if ( calc_u.getBinaryDomainCombinationsSpecificToGenome1().size() != 0 ) {
                return false;
            }
            if ( calc_u.getSharedBinaryDomainCombinations().size() != 1 ) {
                return false;
            }
            if ( calc_u.getAllBinaryDomainCombinations().size() != 3 ) {
                return false;
            }
            if ( calc_u.getBinaryDomainCombinationsSpecificToGenome0()
                    .contains( new BasicBinaryDomainCombination( "v", "u" ) ) ) {
                return false;
            }
            if ( !calc_u.getBinaryDomainCombinationsSpecificToGenome0()
                    .contains( new BasicBinaryDomainCombination( "w", "v" ) ) ) {
                return false;
            }
            if ( !calc_u.getBinaryDomainCombinationsSpecificToGenome0()
                    .contains( new BasicBinaryDomainCombination( "w", "x" ) ) ) {
                return false;
            }
            if ( calc_u.getSharedBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "w", "u" ) ) ) {
                return false;
            }
            if ( !calc_u.getSharedBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "z", "y" ) ) ) {
                return false;
            }
            if ( calc_u.getAllBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "v", "u" ) ) ) {
                return false;
            }
            if ( !calc_u.getAllBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "w", "v" ) ) ) {
                return false;
            }
            if ( !calc_u.getAllBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "w", "x" ) ) ) {
                return false;
            }
            if ( calc_u.getAllBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "w", "u" ) ) ) {
                return false;
            }
            if ( !calc_u.getAllBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "z", "y" ) ) ) {
                return false;
            }
            calc_u.setAllowDomainsToBeIgnored( false );
            if ( calc_u.getAllDomains().size() != 6 ) {
                return false;
            }
            //------------
            calc_u.setAllowDomainsToBeIgnored( true );
            calc_u.deleteAllDomainIdsToIgnore();
            calc_u.addDomainIdToIgnore( new DomainId( "v" ) );
            calc_u.addDomainIdToIgnore( new DomainId( "w" ) );
            calc_u.addDomainIdToIgnore( new DomainId( "other" ) );
            calc_u.addDomainIdToIgnore( new DomainId( "other_too" ) );
            if ( calc_u.getAllDomains().size() != 4 ) {
                return false;
            }
            if ( calc_u.getDomainsSpecificToGenome0().size() != 1 ) {
                return false;
            }
            if ( calc_u.getDomainsSpecificToGenome1().size() != 0 ) {
                return false;
            }
            if ( calc_u.getDomainsSpecificToGenome0().contains( v.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getDomainsSpecificToGenome0().contains( x.getDomainId() ) ) {
                return false;
            }
            if ( calc_u.getSharedDomains().size() != 3 ) {
                return false;
            }
            if ( !calc_u.getSharedDomains().contains( u.getDomainId() ) ) {
                return false;
            }
            if ( calc_u.getSharedDomains().contains( w.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getSharedDomains().contains( y.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getSharedDomains().contains( z.getDomainId() ) ) {
                return false;
            }
            if ( calc_u.getAllDomains().size() != 4 ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( u.getDomainId() ) ) {
                return false;
            }
            if ( calc_u.getAllDomains().contains( w.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( y.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( z.getDomainId() ) ) {
                return false;
            }
            if ( calc_u.getAllDomains().contains( v.getDomainId() ) ) {
                return false;
            }
            if ( !calc_u.getAllDomains().contains( x.getDomainId() ) ) {
                return false;
            }
            if ( calc_u.getBinaryDomainCombinationsSpecificToGenome0().size() != 0 ) {
                return false;
            }
            if ( calc_u.getBinaryDomainCombinationsSpecificToGenome1().size() != 0 ) {
                return false;
            }
            if ( calc_u.getSharedBinaryDomainCombinations().size() != 1 ) {
                return false;
            }
            if ( calc_u.getAllBinaryDomainCombinations().size() != 1 ) {
                return false;
            }
            if ( !calc_u.getSharedBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "y", "z" ) ) ) {
                return false;
            }
            if ( !calc_u.getAllBinaryDomainCombinations().contains( new BasicBinaryDomainCombination( "z", "y" ) ) ) {
                return false;
            }
            if ( !isEqual( calc_u.calculateSharedBinaryDomainCombinationBasedGenomeSimilarityScore(),
                           1.0 - ( 1.0 - 1.0 ) / 1.0 ) ) {
                return false;
            }
            if ( !isEqual( calc_u.calculateSharedDomainsBasedGenomeSimilarityScore(), 1.0 - ( 4.0 - 3.0 ) / 4.0 ) ) {
                return false;
            }
            calc_u.setAllowDomainsToBeIgnored( false );
            if ( !isEqual( calc_u.calculateSharedBinaryDomainCombinationBasedGenomeSimilarityScore(),
                           1.0 - ( 5.0 - 2.0 ) / 5.0 ) ) {
                return false;
            }
            if ( !isEqual( calc_u.calculateSharedDomainsBasedGenomeSimilarityScore(), 1.0 - ( 6.0 - 4.0 ) / 6.0 ) ) {
                return false;
            }
            calc_u.setAllowDomainsToBeIgnored( true );
            if ( !isEqual( calc_u.calculateSharedBinaryDomainCombinationBasedGenomeSimilarityScore(),
                           1.0 - ( 1.0 - 1.0 ) / 1.0 ) ) {
                return false;
            }
            if ( !isEqual( calc_u.calculateSharedDomainsBasedGenomeSimilarityScore(), 1.0 - ( 4.0 - 3.0 ) / 4.0 ) ) {
                return false;
            }
            calc_u.deleteAllDomainIdsToIgnore();
            if ( !isEqual( calc_u.calculateSharedBinaryDomainCombinationBasedGenomeSimilarityScore(),
                           1.0 - ( 5.0 - 2.0 ) / 5.0 ) ) {
                return false;
            }
            if ( !isEqual( calc_u.calculateSharedDomainsBasedGenomeSimilarityScore(), 1.0 - ( 6.0 - 4.0 ) / 6.0 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDomainCombinationCounting( final File test_dir ) {
        try {
            final HmmPfamOutputParser parser = new HmmPfamOutputParser( new File( test_dir
                    + ForesterUtil.getFileSeparator() + "hmmpfam_output2" ), "human", "ls" );
            parser.setEValueMaximum( 0.2 );
            parser.setIgnoreDufs( true );
            parser.setReturnType( HmmPfamOutputParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            final List<Protein> domain_collections = parser.parse();
            final BasicGenomeWideCombinableDomains cdcc = BasicGenomeWideCombinableDomains
                    .createInstance( domain_collections, false, new BasicSpecies( "human" ) );
            CombinableDomains cd = cdcc.get( new DomainId( "A" ) );
            if ( cd.getKeyDomainCount() != 9 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 7 ) {
                return false;
            }
            if ( cd.getNumberOfCombinableDomains() != 11 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "A" ).getDomainId() ) != 2 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "B" ).getDomainId() ) != 6 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "C" ).getDomainId() ) != 4 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "D" ).getDomainId() ) != 3 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "E" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "U" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "V" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "W" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "X" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "Y" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "Z" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "NN" ).getDomainId() ) != 0 ) {
                return false;
            }
            if ( cd.getKeyDomainCount() != 9 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "B" ) );
            if ( cd.getKeyDomainCount() != 12 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 7 ) {
                return false;
            }
            if ( cd.getNumberOfCombinableDomains() != 11 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "A" ).getDomainId() ) != 6 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "B" ).getDomainId() ) != 4 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "C" ).getDomainId() ) != 4 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "D" ).getDomainId() ) != 3 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "E" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "U" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "V" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "W" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "X" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "Y" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "Z" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "NN" ).getDomainId() ) != 0 ) {
                return false;
            }
            if ( cd.getKeyDomainCount() != 12 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "C" ) );
            if ( cd.getKeyDomainCount() != 10 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 7 ) {
                return false;
            }
            if ( cd.getNumberOfCombinableDomains() != 11 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "A" ).getDomainId() ) != 4 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "B" ).getDomainId() ) != 4 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "C" ).getDomainId() ) != 2 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "D" ).getDomainId() ) != 3 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "E" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "U" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "V" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "W" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "X" ).getDomainId() ) != 2 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "Y" ).getDomainId() ) != 2 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "Z" ).getDomainId() ) != 2 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "NN" ).getDomainId() ) != 0 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "D" ) );
            if ( cd.getKeyDomainCount() != 15 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 6 ) {
                return false;
            }
            if ( cd.getNumberOfCombinableDomains() != 11 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "E" ) );
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            if ( cd.getKeyDomainCount() != 1 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 1 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "U" ) );
            if ( cd.getNumberOfCombinableDomains() != 11 ) {
                return false;
            }
            if ( cd.getKeyDomainCount() != 6 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 3 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "V" ) );
            if ( cd.getNumberOfCombinableDomains() != 11 ) {
                return false;
            }
            if ( cd.getKeyDomainCount() != 3 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 2 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "W" ) );
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            if ( cd.getKeyDomainCount() != 2 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 2 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "X" ) );
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            if ( cd.getKeyDomainCount() != 2 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "Y" ) );
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "Z" ) );
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "NN" ) );
            if ( cd.getKeyDomainCount() != 1 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 1 ) {
                return false;
            }
            if ( cd.getNumberOfCombinableDomains() != 0 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "NN" ).getDomainId() ) != 0 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "MM" ) );
            if ( cd.getNumberOfCombinableDomains() != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "MM" ).getDomainId() ) != 0 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "OO" ).getDomainId() ) != 1 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "OO" ) );
            if ( cd.getNumberOfCombinableDomains() != 2 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "OO" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "MM" ).getDomainId() ) != 1 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "QQ" ) );
            if ( cd.getNumberOfCombinableDomains() != 1 ) {
                return false;
            }
            if ( cd.getKeyDomainCount() != 17 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 4 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "QQ" ).getDomainId() ) != 3 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "PP" ) );
            if ( cd.getNumberOfCombinableDomains() != 0 ) {
                return false;
            }
            if ( cd.getKeyDomainCount() != 2 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 2 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "singlet" ) );
            if ( cd.getKeyDomainCount() != 1 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 1 ) {
                return false;
            }
            if ( cd.getNumberOfCombinableDomains() != 0 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "singlet" ).getDomainId() ) != 0 ) {
                return false;
            }
            cd = cdcc.get( new DomainId( "three" ) );
            if ( cd.getKeyDomainCount() != 3 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 1 ) {
                return false;
            }
            if ( cd.getNumberOfCombinableDomains() != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "three" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "so_far_so_bad" ) ) != 0 ) {
                return false;
            }
            // Ignore combinations with same:
            final BasicGenomeWideCombinableDomains cdcc2 = BasicGenomeWideCombinableDomains
                    .createInstance( domain_collections,
                                     true,
                                     new BasicSpecies( "human" ),
                                     null,
                                     DomainCombinationType.BASIC );
            cd = cdcc2.get( new DomainId( "A" ) );
            if ( cd.getKeyDomainCount() != 9 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 7 ) {
                return false;
            }
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "A" ).getDomainId() ) != 0 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "B" ).getDomainId() ) != 6 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "C" ).getDomainId() ) != 4 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "D" ).getDomainId() ) != 3 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new SimpleDomain( "E" ).getDomainId() ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "U" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "V" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "W" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "X" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "Y" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "Z" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "NN" ) ) != 0 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "B" ) );
            if ( cd.getKeyDomainCount() != 12 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 7 ) {
                return false;
            }
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "A" ) ) != 6 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "B" ) ) != 0 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "C" ) ) != 4 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "D" ) ) != 3 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "E" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "U" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "V" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "W" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "X" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "Y" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "Z" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "NN" ) ) != 0 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "C" ) );
            if ( cd.getKeyDomainCount() != 10 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 7 ) {
                return false;
            }
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "A" ) ) != 4 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "B" ) ) != 4 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "C" ) ) != 0 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "D" ) ) != 3 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "E" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "U" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "V" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "W" ) ) != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "X" ) ) != 2 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "Y" ) ) != 2 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "Z" ) ) != 2 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "NN" ) ) != 0 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "D" ) );
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "E" ) );
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            if ( cd.getKeyDomainCount() != 1 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "U" ) );
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "V" ) );
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "W" ) );
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "X" ) );
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "Y" ) );
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "Z" ) );
            if ( cd.getNumberOfCombinableDomains() != 10 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "NN" ) );
            if ( cd.getNumberOfCombinableDomains() != 0 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "NN" ) ) != 0 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "MM" ) );
            if ( cd.getNumberOfCombinableDomains() != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "MM" ) ) != 0 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "OO" ) ) != 1 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "OO" ) );
            if ( cd.getNumberOfCombinableDomains() != 1 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "OO" ) ) != 0 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "MM" ) ) != 1 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "QQ" ) );
            if ( cd.getNumberOfCombinableDomains() != 0 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "QQ" ) ) != 0 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "singlet" ) );
            if ( cd.getKeyDomainCount() != 1 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 1 ) {
                return false;
            }
            if ( cd.getNumberOfCombinableDomains() != 0 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "singlet" ) ) != 0 ) {
                return false;
            }
            cd = cdcc2.get( new DomainId( "three" ) );
            if ( cd.getKeyDomainCount() != 3 ) {
                return false;
            }
            if ( cd.getKeyDomainProteinsCount() != 1 ) {
                return false;
            }
            if ( cd.getNumberOfCombinableDomains() != 0 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "three" ) ) != 0 ) {
                return false;
            }
            if ( cd.getNumberOfProteinsExhibitingCombination( new DomainId( "so_far_so_bad" ) ) != 0 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDomainId() {
        try {
            final DomainId id1 = new DomainId( "a" );
            final DomainId id2 = new DomainId( "a" );
            final DomainId id3 = new DomainId( "A" );
            final DomainId id4 = new DomainId( "b" );
            if ( !id1.equals( id1 ) ) {
                return false;
            }
            if ( id1.getId().equals( "x" ) ) {
                return false;
            }
            if ( id1.getId().equals( null ) ) {
                return false;
            }
            if ( !id1.equals( id2 ) ) {
                return false;
            }
            if ( id1.equals( id3 ) ) {
                return false;
            }
            if ( id1.hashCode() != id1.hashCode() ) {
                return false;
            }
            if ( id1.hashCode() != id2.hashCode() ) {
                return false;
            }
            if ( id1.hashCode() == id3.hashCode() ) {
                return false;
            }
            if ( id1.compareTo( id1 ) != 0 ) {
                return false;
            }
            if ( id1.compareTo( id2 ) != 0 ) {
                return false;
            }
            if ( id1.compareTo( id3 ) != 0 ) {
                return false;
            }
            if ( id1.compareTo( id4 ) >= 0 ) {
                return false;
            }
            if ( id4.compareTo( id1 ) <= 0 ) {
                return false;
            }
            if ( !id4.getId().equals( "b" ) ) {
                return false;
            }
            final DomainId id5 = new DomainId( " C " );
            if ( !id5.getId().equals( "C" ) ) {
                return false;
            }
            if ( id5.equals( id1 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDomainSorting() {
        try {
            final Domain A = new BasicDomain( "A", ( short ) 1, ( short ) 2, ( short ) 1, ( short ) 1, 0.1, -12 );
            final Domain B = new BasicDomain( "B", ( short ) 1, ( short ) 2, ( short ) 1, ( short ) 1, 0.1, -12 );
            final Domain C = new BasicDomain( "C", ( short ) 1, ( short ) 2, ( short ) 1, ( short ) 1, 0.2, -12 );
            final Domain D = new BasicDomain( "D", ( short ) 1, ( short ) 2, ( short ) 1, ( short ) 1, 0.3, -12 );
            final Domain E = new BasicDomain( "E", ( short ) 1, ( short ) 2, ( short ) 1, ( short ) 1, 0.4, -12 );
            final Domain F = new BasicDomain( "F", ( short ) 1, ( short ) 2, ( short ) 1, ( short ) 1, 0.5, -12 );
            final Domain G = new BasicDomain( "G", ( short ) 1, ( short ) 2, ( short ) 1, ( short ) 1, 0.6, -12 );
            final Domain H1 = new BasicDomain( "H", ( short ) 100, ( short ) 200, ( short ) 1, ( short ) 5, 0.7, -12 );
            final Domain H2 = new BasicDomain( "H", ( short ) 300, ( short ) 400, ( short ) 2, ( short ) 5, 0.7, -12 );
            final Domain H3 = new BasicDomain( "H", ( short ) 500, ( short ) 600, ( short ) 3, ( short ) 5, 0.7, -12 );
            final Domain H4 = new BasicDomain( "H", ( short ) 700, ( short ) 800, ( short ) 4, ( short ) 5, 0.7, -12 );
            final Domain H5 = new BasicDomain( "H", ( short ) 700, ( short ) 800, ( short ) 5, ( short ) 5, 0.7, -12 );
            final Domain H6 = new BasicDomain( "H",
                                               ( short ) 1199,
                                               ( short ) 1299,
                                               ( short ) 6,
                                               ( short ) 6,
                                               0.7,
                                               -0.111 );
            final Domain H7 = new BasicDomain( "H7", ( short ) 700, ( short ) 800, ( short ) 5, ( short ) 5, 0.7, -12 );
            final Domain H8 = new BasicDomain( "H7", ( short ) 700, ( short ) 800, ( short ) 5, ( short ) 200, 0.7, -12 );
            final Protein protein = new BasicProtein( "00", "bat" );
            protein.addProteinDomain( H5 );
            protein.addProteinDomain( H2 );
            protein.addProteinDomain( H7 );
            protein.addProteinDomain( H6 );
            protein.addProteinDomain( A );
            protein.addProteinDomain( G );
            protein.addProteinDomain( H4 );
            protein.addProteinDomain( D );
            protein.addProteinDomain( H1 );
            protein.addProteinDomain( C );
            protein.addProteinDomain( E );
            protein.addProteinDomain( F );
            protein.addProteinDomain( B );
            protein.addProteinDomain( H3 );
            protein.addProteinDomain( H7 );
            protein.addProteinDomain( H7 );
            protein.addProteinDomain( H8 );
            final List<Domain> sorted = SurfacingUtil.sortDomainsWithAscendingConfidenceValues( protein );
            if ( sorted.size() != 17 ) {
                return false;
            }
            if ( !sorted.get( 0 ).getDomainId().getId().equals( "A" ) ) {
                return false;
            }
            if ( sorted.get( 0 ).getNumber() != 1 ) {
                return false;
            }
            if ( !sorted.get( 1 ).getDomainId().getId().equals( "B" ) ) {
                return false;
            }
            if ( sorted.get( 1 ).getNumber() != 1 ) {
                return false;
            }
            if ( !sorted.get( 2 ).getDomainId().getId().equals( "C" ) ) {
                return false;
            }
            if ( sorted.get( 2 ).getNumber() != 1 ) {
                return false;
            }
            if ( !sorted.get( 3 ).getDomainId().getId().equals( "D" ) ) {
                return false;
            }
            if ( sorted.get( 3 ).getNumber() != 1 ) {
                return false;
            }
            if ( !sorted.get( 4 ).getDomainId().getId().equals( "E" ) ) {
                return false;
            }
            if ( sorted.get( 4 ).getNumber() != 1 ) {
                return false;
            }
            if ( !sorted.get( 5 ).getDomainId().getId().equals( "F" ) ) {
                return false;
            }
            if ( sorted.get( 5 ).getNumber() != 1 ) {
                return false;
            }
            if ( !sorted.get( 6 ).getDomainId().getId().equals( "G" ) ) {
                return false;
            }
            if ( sorted.get( 6 ).getNumber() != 1 ) {
                return false;
            }
            if ( !sorted.get( 7 ).getDomainId().getId().equals( "H" ) ) {
                return false;
            }
            if ( sorted.get( 7 ).getNumber() != 5 ) {
                return false;
            }
            if ( !sorted.get( 8 ).getDomainId().getId().equals( "H" ) ) {
                return false;
            }
            if ( sorted.get( 8 ).getNumber() != 2 ) {
                return false;
            }
            if ( !sorted.get( 9 ).getDomainId().getId().equals( "H" ) ) {
                return false;
            }
            if ( sorted.get( 9 ).getNumber() != 6 ) {
                return false;
            }
            if ( !sorted.get( 10 ).getDomainId().getId().equals( "H" ) ) {
                return false;
            }
            if ( sorted.get( 10 ).getNumber() != 4 ) {
                return false;
            }
            if ( !sorted.get( 11 ).getDomainId().getId().equals( "H" ) ) {
                return false;
            }
            if ( sorted.get( 11 ).getNumber() != 1 ) {
                return false;
            }
            if ( sorted.get( 11 ).getTotalCount() != 5 ) {
                return false;
            }
            if ( !sorted.get( 12 ).getDomainId().getId().equals( "H" ) ) {
                return false;
            }
            if ( sorted.get( 12 ).getNumber() != 3 ) {
                return false;
            }
            if ( !sorted.get( 13 ).getDomainId().getId().equals( "H7" ) ) {
                return false;
            }
            if ( sorted.get( 13 ).getNumber() != 5 ) {
                return false;
            }
            if ( !sorted.get( 14 ).getDomainId().getId().equals( "H7" ) ) {
                return false;
            }
            if ( sorted.get( 14 ).getNumber() != 5 ) {
                return false;
            }
            if ( !sorted.get( 15 ).getDomainId().getId().equals( "H7" ) ) {
                return false;
            }
            if ( sorted.get( 15 ).getNumber() != 5 ) {
                return false;
            }
            // To check if sorting is stable [as claimed by Sun for
            // Collections.sort( List )]
            if ( !sorted.get( 16 ).getDomainId().getId().equals( "H7" ) ) {
                return false;
            }
            if ( sorted.get( 16 ).getNumber() != 5 ) {
                return false;
            }
            if ( sorted.get( 16 ).getTotalCount() != 200 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testEngulfingOverlapRemoval() {
        try {
            final Domain d0 = new BasicDomain( "d0", 0, 8, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d1 = new BasicDomain( "d1", 0, 1, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d2 = new BasicDomain( "d2", 0, 2, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d3 = new BasicDomain( "d3", 7, 8, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d4 = new BasicDomain( "d4", 7, 9, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d5 = new BasicDomain( "d4", 0, 9, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d6 = new BasicDomain( "d4", 4, 5, ( short ) 1, ( short ) 1, 0.1, 1 );
            final List<Boolean> covered = new ArrayList<Boolean>();
            covered.add( true ); // 0
            covered.add( false ); // 1
            covered.add( true ); // 2
            covered.add( false ); // 3
            covered.add( true ); // 4
            covered.add( true ); // 5
            covered.add( false ); // 6
            covered.add( true ); // 7
            covered.add( true ); // 8
            if ( SurfacingUtil.isEngulfed( d0, covered ) ) {
                return false;
            }
            if ( SurfacingUtil.isEngulfed( d1, covered ) ) {
                return false;
            }
            if ( SurfacingUtil.isEngulfed( d2, covered ) ) {
                return false;
            }
            if ( !SurfacingUtil.isEngulfed( d3, covered ) ) {
                return false;
            }
            if ( SurfacingUtil.isEngulfed( d4, covered ) ) {
                return false;
            }
            if ( SurfacingUtil.isEngulfed( d5, covered ) ) {
                return false;
            }
            if ( !SurfacingUtil.isEngulfed( d6, covered ) ) {
                return false;
            }
            final Domain a = new BasicDomain( "a", 0, 10, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain b = new BasicDomain( "b", 8, 20, ( short ) 1, ( short ) 1, 0.2, 1 );
            final Domain c = new BasicDomain( "c", 15, 16, ( short ) 1, ( short ) 1, 0.3, 1 );
            final Protein abc = new BasicProtein( "abc", "nemve" );
            abc.addProteinDomain( a );
            abc.addProteinDomain( b );
            abc.addProteinDomain( c );
            final Protein abc_r1 = SurfacingUtil.removeOverlappingDomains( 3, false, abc );
            final Protein abc_r2 = SurfacingUtil.removeOverlappingDomains( 3, true, abc );
            if ( abc.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            if ( abc_r1.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            if ( abc_r2.getNumberOfProteinDomains() != 2 ) {
                return false;
            }
            if ( !abc_r2.getProteinDomain( 0 ).getDomainId().getId().equals( "a" ) ) {
                return false;
            }
            if ( !abc_r2.getProteinDomain( 1 ).getDomainId().getId().equals( "b" ) ) {
                return false;
            }
            final Domain d = new BasicDomain( "d", 0, 10, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain e = new BasicDomain( "e", 8, 20, ( short ) 1, ( short ) 1, 0.3, 1 );
            final Domain f = new BasicDomain( "f", 15, 16, ( short ) 1, ( short ) 1, 0.2, 1 );
            final Protein def = new BasicProtein( "def", "nemve" );
            def.addProteinDomain( d );
            def.addProteinDomain( e );
            def.addProteinDomain( f );
            final Protein def_r1 = SurfacingUtil.removeOverlappingDomains( 5, false, def );
            final Protein def_r2 = SurfacingUtil.removeOverlappingDomains( 5, true, def );
            if ( def.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            if ( def_r1.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            if ( def_r2.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            if ( !def_r2.getProteinDomain( 0 ).getDomainId().getId().equals( "d" ) ) {
                return false;
            }
            if ( !def_r2.getProteinDomain( 1 ).getDomainId().getId().equals( "f" ) ) {
                return false;
            }
            if ( !def_r2.getProteinDomain( 2 ).getDomainId().getId().equals( "e" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testGenomeWideCombinableDomains() {
        try {
            final Domain a = new BasicDomain( "a", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain b = new BasicDomain( "b", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain c = new BasicDomain( "c", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain d = new BasicDomain( "d", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain e = new BasicDomain( "e", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain f = new BasicDomain( "f", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain g = new BasicDomain( "g", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain h = new BasicDomain( "h", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain x = new BasicDomain( "x", 23, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Protein eel_0 = new BasicProtein( "0", "eel" );
            final Protein eel_1 = new BasicProtein( "1", "eel" );
            final Protein eel_2 = new BasicProtein( "2", "eel" );
            final Protein eel_3 = new BasicProtein( "3", "eel" );
            final Protein eel_4 = new BasicProtein( "4", "eel" );
            final Protein eel_5 = new BasicProtein( "5", "eel" );
            final Protein eel_6 = new BasicProtein( "6", "eel" );
            eel_1.addProteinDomain( a );
            eel_2.addProteinDomain( a );
            eel_2.addProteinDomain( b );
            eel_3.addProteinDomain( a );
            eel_3.addProteinDomain( a );
            eel_3.addProteinDomain( b );
            eel_4.addProteinDomain( a );
            eel_4.addProteinDomain( b );
            eel_4.addProteinDomain( c );
            eel_4.addProteinDomain( d );
            eel_4.addProteinDomain( e );
            eel_5.addProteinDomain( e );
            eel_5.addProteinDomain( e );
            eel_5.addProteinDomain( f );
            eel_5.addProteinDomain( f );
            eel_5.addProteinDomain( f );
            eel_5.addProteinDomain( f );
            eel_6.addProteinDomain( g );
            eel_6.addProteinDomain( h );
            final List<Protein> protein_list_eel = new ArrayList<Protein>();
            protein_list_eel.add( eel_0 );
            protein_list_eel.add( eel_1 );
            protein_list_eel.add( eel_2 );
            protein_list_eel.add( eel_3 );
            protein_list_eel.add( eel_4 );
            protein_list_eel.add( eel_5 );
            protein_list_eel.add( eel_6 );
            final BasicGenomeWideCombinableDomains eel_not_ignore = BasicGenomeWideCombinableDomains
                    .createInstance( protein_list_eel, false, new BasicSpecies( "eel" ) );
            final BasicGenomeWideCombinableDomains eel_ignore = BasicGenomeWideCombinableDomains
                    .createInstance( protein_list_eel, true, new BasicSpecies( "eel" ) );
            if ( !eel_not_ignore.contains( new DomainId( "a" ) ) ) {
                return false;
            }
            if ( !eel_not_ignore.contains( new DomainId( "b" ) ) ) {
                return false;
            }
            if ( !eel_not_ignore.contains( new DomainId( "c" ) ) ) {
                return false;
            }
            if ( !eel_not_ignore.contains( new DomainId( "d" ) ) ) {
                return false;
            }
            if ( !eel_not_ignore.contains( new DomainId( "e" ) ) ) {
                return false;
            }
            if ( !eel_not_ignore.contains( new DomainId( "f" ) ) ) {
                return false;
            }
            if ( !eel_not_ignore.contains( new DomainId( "g" ) ) ) {
                return false;
            }
            if ( !eel_not_ignore.contains( new DomainId( "h" ) ) ) {
                return false;
            }
            if ( eel_not_ignore.contains( new DomainId( "x" ) ) ) {
                return false;
            }
            if ( !eel_ignore.contains( new DomainId( "a" ) ) ) {
                return false;
            }
            if ( !eel_ignore.contains( new DomainId( "b" ) ) ) {
                return false;
            }
            if ( !eel_ignore.contains( new DomainId( "c" ) ) ) {
                return false;
            }
            if ( !eel_ignore.contains( new DomainId( "d" ) ) ) {
                return false;
            }
            if ( !eel_ignore.contains( new DomainId( "e" ) ) ) {
                return false;
            }
            if ( !eel_ignore.contains( new DomainId( "f" ) ) ) {
                return false;
            }
            if ( !eel_ignore.contains( new DomainId( "g" ) ) ) {
                return false;
            }
            if ( !eel_ignore.contains( new DomainId( "h" ) ) ) {
                return false;
            }
            if ( eel_ignore.contains( new DomainId( "x" ) ) ) {
                return false;
            }
            if ( eel_not_ignore.getSize() != 8 ) {
                return false;
            }
            if ( eel_ignore.getSize() != 8 ) {
                return false;
            }
            if ( eel_not_ignore.get( new DomainId( "a" ) ).getCombinableDomainsIds().size() != 5 ) {
                return false;
            }
            if ( eel_not_ignore.get( new DomainId( "b" ) ).getCombinableDomainsIds().size() != 4 ) {
                return false;
            }
            if ( eel_not_ignore.get( new DomainId( "c" ) ).getCombinableDomainsIds().size() != 4 ) {
                return false;
            }
            if ( eel_not_ignore.get( new DomainId( "d" ) ).getCombinableDomainsIds().size() != 4 ) {
                return false;
            }
            if ( eel_not_ignore.get( new DomainId( "e" ) ).getCombinableDomainsIds().size() != 6 ) {
                return false;
            }
            if ( eel_not_ignore.get( new DomainId( "f" ) ).getCombinableDomainsIds().size() != 2 ) {
                return false;
            }
            if ( eel_not_ignore.get( new DomainId( "g" ) ).getCombinableDomainsIds().size() != 1 ) {
                return false;
            }
            if ( eel_not_ignore.get( new DomainId( "h" ) ).getCombinableDomainsIds().size() != 1 ) {
                return false;
            }
            if ( eel_ignore.get( new DomainId( "a" ) ).getCombinableDomainsIds().size() != 4 ) {
                return false;
            }
            if ( eel_ignore.get( new DomainId( "b" ) ).getCombinableDomainsIds().size() != 4 ) {
                return false;
            }
            if ( eel_ignore.get( new DomainId( "c" ) ).getCombinableDomainsIds().size() != 4 ) {
                return false;
            }
            if ( eel_ignore.get( new DomainId( "d" ) ).getCombinableDomainsIds().size() != 4 ) {
                return false;
            }
            if ( eel_ignore.get( new DomainId( "e" ) ).getCombinableDomainsIds().size() != 5 ) {
                return false;
            }
            if ( eel_ignore.get( new DomainId( "f" ) ).getCombinableDomainsIds().size() != 1 ) {
                return false;
            }
            if ( eel_ignore.get( new DomainId( "g" ) ).getCombinableDomainsIds().size() != 1 ) {
                return false;
            }
            if ( eel_ignore.get( new DomainId( "h" ) ).getCombinableDomainsIds().size() != 1 ) {
                return false;
            }
            if ( eel_not_ignore.getAllDomainIds().size() != 8 ) {
                return false;
            }
            if ( !eel_not_ignore.getAllDomainIds().contains( a.getDomainId() ) ) {
                return false;
            }
            if ( !eel_not_ignore.getAllDomainIds().contains( b.getDomainId() ) ) {
                return false;
            }
            if ( !eel_not_ignore.getAllDomainIds().contains( c.getDomainId() ) ) {
                return false;
            }
            if ( !eel_not_ignore.getAllDomainIds().contains( d.getDomainId() ) ) {
                return false;
            }
            if ( !eel_not_ignore.getAllDomainIds().contains( e.getDomainId() ) ) {
                return false;
            }
            if ( !eel_not_ignore.getAllDomainIds().contains( f.getDomainId() ) ) {
                return false;
            }
            if ( !eel_not_ignore.getAllDomainIds().contains( g.getDomainId() ) ) {
                return false;
            }
            if ( !eel_not_ignore.getAllDomainIds().contains( h.getDomainId() ) ) {
                return false;
            }
            if ( eel_not_ignore.getAllDomainIds().contains( x.getDomainId() ) ) {
                return false;
            }
            if ( eel_ignore.getAllDomainIds().size() != 8 ) {
                return false;
            }
            if ( !eel_ignore.getAllDomainIds().contains( a.getDomainId() ) ) {
                return false;
            }
            if ( !eel_ignore.getAllDomainIds().contains( b.getDomainId() ) ) {
                return false;
            }
            if ( !eel_ignore.getAllDomainIds().contains( c.getDomainId() ) ) {
                return false;
            }
            if ( !eel_ignore.getAllDomainIds().contains( d.getDomainId() ) ) {
                return false;
            }
            if ( !eel_ignore.getAllDomainIds().contains( e.getDomainId() ) ) {
                return false;
            }
            if ( !eel_ignore.getAllDomainIds().contains( f.getDomainId() ) ) {
                return false;
            }
            if ( !eel_ignore.getAllDomainIds().contains( g.getDomainId() ) ) {
                return false;
            }
            if ( !eel_ignore.getAllDomainIds().contains( h.getDomainId() ) ) {
                return false;
            }
            if ( eel_ignore.getAllDomainIds().contains( x.getDomainId() ) ) {
                return false;
            }
            final SortedSet<BinaryDomainCombination> bc0 = eel_not_ignore.toBinaryDomainCombinations();
            if ( bc0.size() != 15 ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "a", "a" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "a", "b" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "b", "a" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "a", "c" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "a", "d" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "a", "e" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "b", "c" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "b", "d" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "b", "e" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "c", "d" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "c", "e" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "d", "e" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "e", "f" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "e", "e" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "f", "f" ) ) ) {
                return false;
            }
            if ( !bc0.contains( new BasicBinaryDomainCombination( "g", "h" ) ) ) {
                return false;
            }
            if ( bc0.contains( new BasicBinaryDomainCombination( "f", "a" ) ) ) {
                return false;
            }
            if ( bc0.contains( new BasicBinaryDomainCombination( "f", "b" ) ) ) {
                return false;
            }
            if ( bc0.contains( new BasicBinaryDomainCombination( "a", "h" ) ) ) {
                return false;
            }
            if ( bc0.contains( new BasicBinaryDomainCombination( "a", "g" ) ) ) {
                return false;
            }
            final SortedSet<BinaryDomainCombination> bc1 = eel_ignore.toBinaryDomainCombinations();
            if ( bc1.size() != 12 ) {
                return false;
            }
            if ( bc1.contains( new BasicBinaryDomainCombination( "a", "a" ) ) ) {
                return false;
            }
            if ( !bc1.contains( new BasicBinaryDomainCombination( "a", "b" ) ) ) {
                return false;
            }
            if ( !bc1.contains( new BasicBinaryDomainCombination( "b", "a" ) ) ) {
                return false;
            }
            if ( !bc1.contains( new BasicBinaryDomainCombination( "a", "c" ) ) ) {
                return false;
            }
            if ( !bc1.contains( new BasicBinaryDomainCombination( "a", "d" ) ) ) {
                return false;
            }
            if ( !bc1.contains( new BasicBinaryDomainCombination( "a", "e" ) ) ) {
                return false;
            }
            if ( !bc1.contains( new BasicBinaryDomainCombination( "b", "c" ) ) ) {
                return false;
            }
            if ( !bc1.contains( new BasicBinaryDomainCombination( "b", "d" ) ) ) {
                return false;
            }
            if ( !bc1.contains( new BasicBinaryDomainCombination( "b", "e" ) ) ) {
                return false;
            }
            if ( !bc1.contains( new BasicBinaryDomainCombination( "c", "d" ) ) ) {
                return false;
            }
            if ( !bc1.contains( new BasicBinaryDomainCombination( "c", "e" ) ) ) {
                return false;
            }
            if ( !bc1.contains( new BasicBinaryDomainCombination( "d", "e" ) ) ) {
                return false;
            }
            if ( !bc1.contains( new BasicBinaryDomainCombination( "e", "f" ) ) ) {
                return false;
            }
            if ( !bc1.contains( new BasicBinaryDomainCombination( "g", "h" ) ) ) {
                return false;
            }
            if ( bc1.contains( new BasicBinaryDomainCombination( "e", "e" ) ) ) {
                return false;
            }
            if ( bc1.contains( new BasicBinaryDomainCombination( "f", "f" ) ) ) {
                return false;
            }
            if ( bc1.contains( new BasicBinaryDomainCombination( "f", "a" ) ) ) {
                return false;
            }
            if ( bc1.contains( new BasicBinaryDomainCombination( "f", "b" ) ) ) {
                return false;
            }
            if ( bc1.contains( new BasicBinaryDomainCombination( "a", "g" ) ) ) {
                return false;
            }
            if ( bc1.contains( new BasicBinaryDomainCombination( "b", "g" ) ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testHmmPfamOutputParser( final File test_dir ) {
        try {
            final HmmPfamOutputParser parser = new HmmPfamOutputParser( new File( test_dir
                    + ForesterUtil.getFileSeparator() + "hmmpfam_output" ), "human", "ls" );
            parser.setEValueMaximum( 0.2 );
            parser.setIgnoreDufs( true );
            parser.setReturnType( HmmPfamOutputParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            List<?> domain_collections = null;
            domain_collections = parser.parse();
            if ( parser.getDomainsEncountered() != 4 ) {
                return false;
            }
            if ( parser.getDomainsIgnoredDueToDuf() != 0 ) {
                return false;
            }
            if ( parser.getDomainsIgnoredDueToEval() != 1 ) {
                return false;
            }
            if ( parser.getDomainsIgnoredDueToOverlap() != 0 ) {
                return false;
            }
            if ( parser.getDomainsStored() != 3 ) {
                return false;
            }
            if ( domain_collections.size() != 1 ) {
                return false;
            }
            final Protein pdc = ( Protein ) domain_collections.get( 0 );
            if ( !pdc.getProteinId().equals( new ProteinId( "ENSP00000285681" ) ) ) {
                return false;
            }
            if ( !pdc.getSpecies().getSpeciesId().equals( "human" ) ) {
                return false;
            }
            if ( pdc.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            if ( !pdc.getAccession().equals( "acc_ENSP00000285681" ) ) {
                return false;
            }
            if ( !pdc
                    .getDescription()
                    .equals( "pep:known chromosome:NCBI36:21:16024215:16174248:1 gene:ENSG00000155313 transcript:ENST00000285681" ) ) {
                return false;
            }
            final List<Domain> uba = pdc.getProteinDomains( new DomainId( "UBA" ) );
            final List<Domain> uim = pdc.getProteinDomains( new DomainId( "UIM" ) );
            final List<Domain> uch = pdc.getProteinDomains( new DomainId( "UCH" ) );
            if ( uba.size() != 1 ) {
                return false;
            }
            if ( uim.size() != 2 ) {
                return false;
            }
            if ( uch.size() != 0 ) {
                return false;
            }
            final BasicDomain uim_domain = ( BasicDomain ) uim.get( 1 );
            if ( !uim_domain.getDomainId().equals( new DomainId( "UIM" ) ) ) {
                return false;
            }
            if ( uim_domain.getTotalCount() != 2 ) {
                return false;
            }
            final BasicDomain uba_domain = ( BasicDomain ) uba.get( 0 );
            if ( !uba_domain.getDomainId().equals( new DomainId( "UBA" ) ) ) {
                return false;
            }
            if ( uba_domain.getNumber() != 1 ) {
                return false;
            }
            if ( uba_domain.getTotalCount() != 1 ) {
                return false;
            }
            if ( uba_domain.getFrom() != 16 ) {
                return false;
            }
            if ( uba_domain.getTo() != 57 ) {
                return false;
            }
            if ( !Test.isEqual( uba_domain.getPerSequenceEvalue(), 0.00084 ) ) {
                return false;
            }
            if ( !Test.isEqual( uba_domain.getPerSequenceScore(), 23.2 ) ) {
                return false;
            }
            final HmmPfamOutputParser parser2 = new HmmPfamOutputParser( new File( test_dir
                    + ForesterUtil.getFileSeparator() + "hmmpfam_output_short" ), "human", "ls" );
            parser2.setEValueMaximum( 0.2 );
            parser2.setIgnoreDufs( true );
            parser2.setReturnType( HmmPfamOutputParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            List<Protein> domain_collections2 = null;
            domain_collections2 = parser2.parse();
            if ( parser2.getDomainsEncountered() != 4 ) {
                return false;
            }
            if ( parser.getDomainsIgnoredDueToDuf() != 0 ) {
                return false;
            }
            if ( parser.getDomainsIgnoredDueToEval() != 1 ) {
                return false;
            }
            if ( parser.getDomainsIgnoredDueToOverlap() != 0 ) {
                return false;
            }
            if ( parser2.getDomainsStored() != 3 ) {
                return false;
            }
            if ( domain_collections2.size() != 1 ) {
                return false;
            }
            final Protein pdc2 = domain_collections2.get( 0 );
            if ( !pdc2.getProteinId().getId().equals( "ENSP00000285681" ) ) {
                return false;
            }
            if ( !pdc2.getSpecies().getSpeciesId().equals( "human" ) ) {
                return false;
            }
            if ( !pdc2.getName().equals( "" ) ) {
                return false;
            }
            if ( !pdc2.getAccession().equals( "223" ) ) {
                return false;
            }
            if ( !pdc2
                    .getDescription()
                    .equals( "pep:known chromosome:NCBI36:21:16024215:16174248:1 gene:ENSG00000155313 transcript:ENST00000285681" ) ) {
                return false;
            }
            if ( pdc2.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            final List<Domain> uba2 = pdc2.getProteinDomains( new DomainId( "UBA" ) );
            final List<Domain> uim2 = pdc2.getProteinDomains( new DomainId( "UIM" ) );
            final List<Domain> uch2 = pdc2.getProteinDomains( new DomainId( "UCH" ) );
            if ( uba2.size() != 1 ) {
                return false;
            }
            if ( uim2.size() != 2 ) {
                return false;
            }
            if ( uch2.size() != 0 ) {
                return false;
            }
            final BasicDomain uim_domain2 = ( BasicDomain ) uim2.get( 1 );
            if ( !uim_domain2.getDomainId().getId().equals( "UIM" ) ) {
                return false;
            }
            if ( uim_domain2.getTotalCount() != 2 ) {
                return false;
            }
            final BasicDomain uba_domain2 = ( BasicDomain ) uba2.get( 0 );
            if ( !uba_domain2.getDomainId().getId().equals( "UBA" ) ) {
                return false;
            }
            if ( uba_domain2.getNumber() != 1 ) {
                return false;
            }
            if ( uba_domain2.getTotalCount() != 1 ) {
                return false;
            }
            if ( uba_domain2.getFrom() != 16 ) {
                return false;
            }
            if ( uba_domain2.getTo() != 57 ) {
                return false;
            }
            if ( !Test.isEqual( uba_domain2.getPerSequenceEvalue(), 0.00084 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testHmmPfamOutputParserWithFilter( final File test_dir ) {
        try {
            HmmPfamOutputParser parser = new HmmPfamOutputParser( new File( test_dir + ForesterUtil.getFileSeparator()
                    + "hmmpfam_output3" ), "human", "ls" );
            parser.setEValueMaximum( 0.2 );
            parser.setIgnoreDufs( true );
            parser.setReturnType( HmmPfamOutputParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            List<Protein> proteins = null;
            proteins = parser.parse();
            if ( parser.getProteinsIgnoredDueToFilter() != 0 ) {
                return false;
            }
            if ( proteins.size() != 4 ) {
                return false;
            }
            //
            Set<DomainId> filter = new TreeSet<DomainId>();
            filter.add( new DomainId( "beauty" ) );
            filter.add( new DomainId( "strange" ) );
            parser = new HmmPfamOutputParser( new File( test_dir + ForesterUtil.getFileSeparator() + "hmmpfam_output3" ),
                                              "human",
                                              "ls",
                                              filter,
                                              HmmPfamOutputParser.FilterType.NEGATIVE_PROTEIN );
            parser.setEValueMaximum( 0.2 );
            parser.setIgnoreDufs( true );
            parser.setReturnType( HmmPfamOutputParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            proteins = null;
            proteins = parser.parse();
            if ( parser.getProteinsIgnoredDueToFilter() != 0 ) {
                return false;
            }
            if ( proteins.size() != 4 ) {
                return false;
            }
            //
            filter = new TreeSet<DomainId>();
            filter.add( new DomainId( "beauty" ) );
            filter.add( new DomainId( "strange" ) );
            parser = new HmmPfamOutputParser( new File( test_dir + ForesterUtil.getFileSeparator() + "hmmpfam_output3" ),
                                              "human",
                                              "ls",
                                              filter,
                                              HmmPfamOutputParser.FilterType.POSITIVE_PROTEIN );
            parser.setEValueMaximum( 0.2 );
            parser.setIgnoreDufs( true );
            parser.setReturnType( HmmPfamOutputParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            proteins = null;
            proteins = parser.parse();
            if ( parser.getProteinsIgnoredDueToFilter() != 4 ) {
                return false;
            }
            if ( proteins.size() != 0 ) {
                return false;
            }
            //
            filter = new TreeSet<DomainId>();
            filter.add( new DomainId( "UIM" ) );
            filter.add( new DomainId( "A" ) );
            filter.add( new DomainId( "C" ) );
            parser = new HmmPfamOutputParser( new File( test_dir + ForesterUtil.getFileSeparator() + "hmmpfam_output3" ),
                                              "human",
                                              "ls",
                                              filter,
                                              HmmPfamOutputParser.FilterType.POSITIVE_PROTEIN );
            parser.setEValueMaximum( 0.2 );
            parser.setIgnoreDufs( true );
            parser.setReturnType( HmmPfamOutputParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            proteins = null;
            proteins = parser.parse();
            if ( parser.getProteinsIgnoredDueToFilter() != 0 ) {
                return false;
            }
            if ( proteins.size() != 4 ) {
                return false;
            }
            //
            filter = new TreeSet<DomainId>();
            filter.add( new DomainId( "UIM" ) );
            filter.add( new DomainId( "A" ) );
            filter.add( new DomainId( "C" ) );
            filter.add( new DomainId( "X" ) );
            parser = new HmmPfamOutputParser( new File( test_dir + ForesterUtil.getFileSeparator() + "hmmpfam_output3" ),
                                              "human",
                                              "ls",
                                              filter,
                                              HmmPfamOutputParser.FilterType.NEGATIVE_DOMAIN );
            parser.setEValueMaximum( 0.2 );
            parser.setIgnoreDufs( true );
            parser.setReturnType( HmmPfamOutputParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            proteins = null;
            proteins = parser.parse();
            if ( parser.getDomainsIgnoredDueToNegativeDomainFilter() != 7 ) {
                return false;
            }
            if ( proteins.size() != 3 ) {
                return false;
            }
            //
            filter = new TreeSet<DomainId>();
            filter.add( new DomainId( "UIM" ) );
            filter.add( new DomainId( "A" ) );
            filter.add( new DomainId( "C" ) );
            parser = new HmmPfamOutputParser( new File( test_dir + ForesterUtil.getFileSeparator() + "hmmpfam_output3" ),
                                              "human",
                                              "ls",
                                              filter,
                                              HmmPfamOutputParser.FilterType.NEGATIVE_PROTEIN );
            parser.setEValueMaximum( 0.2 );
            parser.setIgnoreDufs( true );
            parser.setReturnType( HmmPfamOutputParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            proteins = null;
            proteins = parser.parse();
            if ( parser.getProteinsIgnoredDueToFilter() != 4 ) {
                return false;
            }
            if ( proteins.size() != 0 ) {
                return false;
            }
            //
            filter = new TreeSet<DomainId>();
            filter.add( new DomainId( "UIM" ) );
            parser = new HmmPfamOutputParser( new File( test_dir + ForesterUtil.getFileSeparator() + "hmmpfam_output3" ),
                                              "human",
                                              "ls",
                                              filter,
                                              HmmPfamOutputParser.FilterType.NEGATIVE_PROTEIN );
            parser.setEValueMaximum( 0.2 );
            parser.setIgnoreDufs( true );
            parser.setReturnType( HmmPfamOutputParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            proteins = null;
            proteins = parser.parse();
            if ( parser.getProteinsIgnoredDueToFilter() != 1 ) {
                return false;
            }
            if ( parser.getProteinsStored() != 3 ) {
                return false;
            }
            if ( proteins.size() != 3 ) {
                return false;
            }
            //
            filter = new TreeSet<DomainId>();
            filter.add( new DomainId( "UIM" ) );
            parser = new HmmPfamOutputParser( new File( test_dir + ForesterUtil.getFileSeparator() + "hmmpfam_output3" ),
                                              "human",
                                              "ls",
                                              filter,
                                              HmmPfamOutputParser.FilterType.POSITIVE_PROTEIN );
            parser.setEValueMaximum( 0.2 );
            parser.setIgnoreDufs( true );
            parser.setReturnType( HmmPfamOutputParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            proteins = null;
            proteins = parser.parse();
            if ( parser.getProteinsIgnoredDueToFilter() != 3 ) {
                return false;
            }
            if ( parser.getProteinsStored() != 1 ) {
                return false;
            }
            if ( proteins.size() != 1 ) {
                return false;
            }
            //
            filter = new TreeSet<DomainId>();
            filter.add( new DomainId( "A" ) );
            filter.add( new DomainId( "C" ) );
            parser = new HmmPfamOutputParser( new File( test_dir + ForesterUtil.getFileSeparator() + "hmmpfam_output3" ),
                                              "human",
                                              "ls",
                                              filter,
                                              HmmPfamOutputParser.FilterType.POSITIVE_PROTEIN );
            parser.setEValueMaximum( 0.2 );
            parser.setIgnoreDufs( true );
            parser.setReturnType( HmmPfamOutputParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            proteins = null;
            proteins = parser.parse();
            if ( parser.getDomainsEncountered() != 11 ) {
                return false;
            }
            if ( parser.getProteinsEncountered() != 4 ) {
                return false;
            }
            if ( parser.getProteinsIgnoredDueToFilter() != 1 ) {
                return false;
            }
            if ( parser.getProteinsStored() != 3 ) {
                return false;
            }
            if ( proteins.size() != 3 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testOverlapRemoval() {
        try {
            final Domain d0 = new BasicDomain( "d0", ( short ) 2, ( short ) 5, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d1 = new BasicDomain( "d1", ( short ) 7, ( short ) 10, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d2 = new BasicDomain( "d2", ( short ) 0, ( short ) 20, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d3 = new BasicDomain( "d3", ( short ) 9, ( short ) 10, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Domain d4 = new BasicDomain( "d4", ( short ) 7, ( short ) 8, ( short ) 1, ( short ) 1, 0.1, 1 );
            final List<Boolean> covered = new ArrayList<Boolean>();
            covered.add( true ); // 0
            covered.add( false ); // 1
            covered.add( true ); // 2
            covered.add( false ); // 3
            covered.add( true ); // 4
            covered.add( true ); // 5
            covered.add( false ); // 6
            covered.add( true ); // 7
            covered.add( true ); // 8
            if ( SurfacingUtil.calculateOverlap( d0, covered ) != 3 ) {
                return false;
            }
            if ( SurfacingUtil.calculateOverlap( d1, covered ) != 2 ) {
                return false;
            }
            if ( SurfacingUtil.calculateOverlap( d2, covered ) != 6 ) {
                return false;
            }
            if ( SurfacingUtil.calculateOverlap( d3, covered ) != 0 ) {
                return false;
            }
            if ( SurfacingUtil.calculateOverlap( d4, covered ) != 2 ) {
                return false;
            }
            final Domain a = new BasicDomain( "a", ( short ) 2, ( short ) 5, ( short ) 1, ( short ) 1, 0.01, 1 );
            final Domain b = new BasicDomain( "b", ( short ) 2, ( short ) 10, ( short ) 1, ( short ) 1, 0.1, 1 );
            final Protein ab = new BasicProtein( "ab", "varanus" );
            ab.addProteinDomain( a );
            ab.addProteinDomain( b );
            final Protein ab_s0 = SurfacingUtil.removeOverlappingDomains( 3, false, ab );
            if ( ab.getNumberOfProteinDomains() != 2 ) {
                return false;
            }
            if ( ab_s0.getNumberOfProteinDomains() != 1 ) {
                return false;
            }
            if ( !ab_s0.getProteinDomain( 0 ).getDomainId().getId().equals( "a" ) ) {
                return false;
            }
            final Protein ab_s1 = SurfacingUtil.removeOverlappingDomains( 4, false, ab );
            if ( ab.getNumberOfProteinDomains() != 2 ) {
                return false;
            }
            if ( ab_s1.getNumberOfProteinDomains() != 2 ) {
                return false;
            }
            final Domain c = new BasicDomain( "c", ( short ) 20000, ( short ) 20500, ( short ) 1, ( short ) 1, 10, 1 );
            final Domain d = new BasicDomain( "d",
                                              ( short ) 10000,
                                              ( short ) 10500,
                                              ( short ) 1,
                                              ( short ) 1,
                                              0.0000001,
                                              1 );
            final Domain e = new BasicDomain( "e", ( short ) 5000, ( short ) 5500, ( short ) 1, ( short ) 1, 0.0001, 1 );
            final Protein cde = new BasicProtein( "cde", "varanus" );
            cde.addProteinDomain( c );
            cde.addProteinDomain( d );
            cde.addProteinDomain( e );
            final Protein cde_s0 = SurfacingUtil.removeOverlappingDomains( 0, false, cde );
            if ( cde.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            if ( cde_s0.getNumberOfProteinDomains() != 3 ) {
                return false;
            }
            final Domain f = new BasicDomain( "f", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 10, 1 );
            final Domain g = new BasicDomain( "g", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 0.01, 1 );
            final Domain h = new BasicDomain( "h", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 0.0001, 1 );
            final Domain i = new BasicDomain( "i", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 0.5, 1 );
            final Domain i2 = new BasicDomain( "i", ( short ) 5, ( short ) 30, ( short ) 1, ( short ) 1, 0.5, 10 );
            final Protein fghi = new BasicProtein( "fghi", "varanus" );
            fghi.addProteinDomain( f );
            fghi.addProteinDomain( g );
            fghi.addProteinDomain( h );
            fghi.addProteinDomain( i );
            fghi.addProteinDomain( i );
            fghi.addProteinDomain( i );
            fghi.addProteinDomain( i2 );
            final Protein fghi_s0 = SurfacingUtil.removeOverlappingDomains( 10, false, fghi );
            if ( fghi.getNumberOfProteinDomains() != 7 ) {
                return false;
            }
            if ( fghi_s0.getNumberOfProteinDomains() != 1 ) {
                return false;
            }
            if ( !fghi_s0.getProteinDomain( 0 ).getDomainId().getId().equals( "h" ) ) {
                return false;
            }
            final Protein fghi_s1 = SurfacingUtil.removeOverlappingDomains( 11, false, fghi );
            if ( fghi.getNumberOfProteinDomains() != 7 ) {
                return false;
            }
            if ( fghi_s1.getNumberOfProteinDomains() != 7 ) {
                return false;
            }
            final Domain j = new BasicDomain( "j", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 10, 1 );
            final Domain k = new BasicDomain( "k", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 0.01, 1 );
            final Domain l = new BasicDomain( "l", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 1, 0.0001, 1 );
            final Domain m = new BasicDomain( "m", ( short ) 10, ( short ) 20, ( short ) 1, ( short ) 4, 0.5, 1 );
            final Domain m0 = new BasicDomain( "m", ( short ) 10, ( short ) 20, ( short ) 2, ( short ) 4, 0.5, 1 );
            final Domain m1 = new BasicDomain( "m", ( short ) 10, ( short ) 20, ( short ) 3, ( short ) 4, 0.5, 1 );
            final Domain m2 = new BasicDomain( "m", ( short ) 5, ( short ) 30, ( short ) 4, ( short ) 4, 0.5, 10 );
            final Protein jklm = new BasicProtein( "jklm", "varanus" );
            jklm.addProteinDomain( j );
            jklm.addProteinDomain( k );
            jklm.addProteinDomain( l );
            jklm.addProteinDomain( m );
            jklm.addProteinDomain( m0 );
            jklm.addProteinDomain( m1 );
            jklm.addProteinDomain( m2 );
            final Protein jklm_s0 = SurfacingUtil.removeOverlappingDomains( 10, false, jklm );
            if ( jklm.getNumberOfProteinDomains() != 7 ) {
                return false;
            }
            if ( jklm_s0.getNumberOfProteinDomains() != 1 ) {
                return false;
            }
            if ( !jklm_s0.getProteinDomain( 0 ).getDomainId().getId().equals( "l" ) ) {
                return false;
            }
            final Protein jklm_s1 = SurfacingUtil.removeOverlappingDomains( 11, false, jklm );
            if ( jklm.getNumberOfProteinDomains() != 7 ) {
                return false;
            }
            if ( jklm_s1.getNumberOfProteinDomains() != 7 ) {
                return false;
            }
            final Domain only = new BasicDomain( "only", ( short ) 5, ( short ) 30, ( short ) 4, ( short ) 4, 0.5, 10 );
            final Protein od = new BasicProtein( "od", "varanus" );
            od.addProteinDomain( only );
            final Protein od_s0 = SurfacingUtil.removeOverlappingDomains( 0, false, od );
            if ( od.getNumberOfProteinDomains() != 1 ) {
                return false;
            }
            if ( od_s0.getNumberOfProteinDomains() != 1 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testParsimony() {
        try {
            final BinaryStates X = BinaryStates.PRESENT;
            final BinaryStates O = BinaryStates.ABSENT;
            final GainLossStates G = GainLossStates.GAIN;
            final GainLossStates L = GainLossStates.LOSS;
            final GainLossStates A = GainLossStates.UNCHANGED_ABSENT;
            final GainLossStates P = GainLossStates.UNCHANGED_PRESENT;
            final Domain a = new BasicDomain( "A", 1, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain b = new BasicDomain( "B", 2, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain c = new BasicDomain( "C", 3, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain d = new BasicDomain( "D", 1, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain e = new BasicDomain( "E", 2, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain f = new BasicDomain( "F", 3, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain g = new BasicDomain( "G", 1, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain h = new BasicDomain( "H", 2, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain i = new BasicDomain( "I", 3, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain j = new BasicDomain( "J", 1, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain l = new BasicDomain( "L", 3, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain m = new BasicDomain( "M", 1, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain n = new BasicDomain( "N", 2, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain o = new BasicDomain( "O", 3, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain p = new BasicDomain( "P", 1, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain q = new BasicDomain( "Q", 2, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain r = new BasicDomain( "R", 3, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            // 1 a-a a-b a-c e-f-g-h l-m
            // 2 a-b a-c e-f-g-i n-o
            // 3 a-b a-d e-f-g-j p-q
            // 4 a-b a-d p-r
            // 1 a-a a-b a-c e-f e-g e-h f-g f-h g-h l-m
            // 2 a-b a-c e-f e-g e-i f-g f-i g-i n-o
            // 3 a-b a-d e-f e-g e-j f-g f-j g-j p-q
            // 4 a-b a-d p-r
            // 1 a b c e f g h l m
            // 2 a b c e f g i n o
            // 3 a b d e f g j p q
            // 4 a b d p r
            final Protein aa1 = new BasicProtein( "aa1", "one" );
            aa1.addProteinDomain( a );
            aa1.addProteinDomain( a );
            final Protein ab1 = new BasicProtein( "ab1", "one" );
            ab1.addProteinDomain( a );
            ab1.addProteinDomain( b );
            final Protein ac1 = new BasicProtein( "ac1", "one" );
            ac1.addProteinDomain( a );
            ac1.addProteinDomain( c );
            final Protein efgh1 = new BasicProtein( "efgh1", "one" );
            efgh1.addProteinDomain( e );
            efgh1.addProteinDomain( f );
            efgh1.addProteinDomain( g );
            efgh1.addProteinDomain( h );
            final Protein lm1 = new BasicProtein( "lm1", "one" );
            lm1.addProteinDomain( l );
            lm1.addProteinDomain( m );
            final Protein ab2 = new BasicProtein( "ab2", "two" );
            ab2.addProteinDomain( a );
            ab2.addProteinDomain( b );
            final Protein ac2 = new BasicProtein( "ac2", "two" );
            ac2.addProteinDomain( a );
            ac2.addProteinDomain( c );
            final Protein efgi2 = new BasicProtein( "efgi2", "two" );
            efgi2.addProteinDomain( e );
            efgi2.addProteinDomain( f );
            efgi2.addProteinDomain( g );
            efgi2.addProteinDomain( i );
            final Protein no2 = new BasicProtein( "no2", "two" );
            no2.addProteinDomain( n );
            no2.addProteinDomain( o );
            final Protein ab3 = new BasicProtein( "ab3", "three" );
            ab3.addProteinDomain( a );
            ab3.addProteinDomain( b );
            final Protein ad3 = new BasicProtein( "ad3", "three" );
            ad3.addProteinDomain( a );
            ad3.addProteinDomain( d );
            final Protein efgj3 = new BasicProtein( "efgj3", "three" );
            efgj3.addProteinDomain( e );
            efgj3.addProteinDomain( f );
            efgj3.addProteinDomain( g );
            efgj3.addProteinDomain( j );
            final Protein pq3 = new BasicProtein( "pq3", "three" );
            pq3.addProteinDomain( p );
            pq3.addProteinDomain( q );
            final Protein ab4 = new BasicProtein( "ab4", "four" );
            ab4.addProteinDomain( a );
            ab4.addProteinDomain( b );
            final Protein ad4 = new BasicProtein( "ad4", "four" );
            ad4.addProteinDomain( a );
            ad4.addProteinDomain( d );
            final Protein pr4 = new BasicProtein( "pr4", "four" );
            pr4.addProteinDomain( p );
            pr4.addProteinDomain( r );
            final List<Protein> one_list = new ArrayList<Protein>();
            one_list.add( aa1 );
            one_list.add( ab1 );
            one_list.add( ac1 );
            one_list.add( efgh1 );
            one_list.add( lm1 );
            final List<Protein> two_list = new ArrayList<Protein>();
            two_list.add( ab2 );
            two_list.add( ac2 );
            two_list.add( efgi2 );
            two_list.add( no2 );
            final List<Protein> three_list = new ArrayList<Protein>();
            three_list.add( ab3 );
            three_list.add( ad3 );
            three_list.add( efgj3 );
            three_list.add( pq3 );
            final List<Protein> four_list = new ArrayList<Protein>();
            four_list.add( ab4 );
            four_list.add( ad4 );
            four_list.add( pr4 );
            final GenomeWideCombinableDomains one = BasicGenomeWideCombinableDomains
                    .createInstance( one_list, false, new BasicSpecies( "one" ) );
            final GenomeWideCombinableDomains two = BasicGenomeWideCombinableDomains
                    .createInstance( two_list, false, new BasicSpecies( "two" ) );
            final GenomeWideCombinableDomains three = BasicGenomeWideCombinableDomains
                    .createInstance( three_list, false, new BasicSpecies( "three" ) );
            final GenomeWideCombinableDomains four = BasicGenomeWideCombinableDomains
                    .createInstance( four_list, false, new BasicSpecies( "four" ) );
            final List<GenomeWideCombinableDomains> gwcd_list = new ArrayList<GenomeWideCombinableDomains>();
            gwcd_list.add( one );
            gwcd_list.add( two );
            gwcd_list.add( three );
            gwcd_list.add( four );
            final CharacterStateMatrix<BinaryStates> matrix_d = DomainParsimonyCalculator
                    .createMatrixOfDomainPresenceOrAbsence( gwcd_list );
            final CharacterStateMatrix<BinaryStates> matrix_bc = DomainParsimonyCalculator
                    .createMatrixOfBinaryDomainCombinationPresenceOrAbsence( gwcd_list );
            // 1 a b c e f g h l m
            // 2 a b c e f g i n o
            // 3 a b d e f g j p q
            // 4 a b d p r
            if ( matrix_d.getState( 0, 0 ) != X ) {
                return false;
            }
            if ( matrix_d.getState( 0, 1 ) != X ) {
                return false;
            }
            if ( matrix_d.getState( 0, 2 ) != X ) {
                return false;
            }
            if ( matrix_d.getState( 0, 3 ) != O ) {
                return false;
            }
            if ( matrix_d.getState( 0, 4 ) != X ) {
                return false;
            }
            if ( matrix_d.getState( 0, 5 ) != X ) {
                return false;
            }
            if ( matrix_d.getState( 0, 6 ) != X ) {
                return false;
            }
            if ( matrix_d.getState( 0, 7 ) != X ) {
                return false;
            }
            if ( matrix_d.getState( 0, 8 ) != O ) {
                return false;
            }
            // 1 a-a a-b a-c e-f e-g e-h f-g f-h g-h l-m
            // 2 a-b a-c e-f e-g e-i f-g f-i g-i n-o
            // 3 a-b a-d e-f e-g e-j f-g f-j g-j p-q
            // 4 a-b a-d p-r
            if ( matrix_bc.getState( 0, 0 ) != X ) {
                return false;
            }
            if ( matrix_bc.getState( 0, 1 ) != X ) {
                return false;
            }
            if ( matrix_bc.getState( 0, 2 ) != X ) {
                return false;
            }
            if ( matrix_bc.getState( 0, 3 ) != O ) {
                return false;
            }
            if ( matrix_bc.getState( 0, 4 ) != X ) {
                return false;
            }
            if ( matrix_bc.getState( 1, 0 ) != O ) {
                return false;
            }
            if ( matrix_bc.getState( 1, 1 ) != X ) {
                return false;
            }
            if ( matrix_bc.getState( 1, 2 ) != X ) {
                return false;
            }
            if ( matrix_bc.getState( 1, 3 ) != O ) {
                return false;
            }
            if ( matrix_bc.getState( 1, 4 ) != X ) {
                return false;
            }
            if ( matrix_bc.getState( 2, 0 ) != O ) {
                return false;
            }
            if ( matrix_bc.getState( 2, 1 ) != X ) {
                return false;
            }
            if ( matrix_bc.getState( 2, 2 ) != O ) {
                return false;
            }
            if ( matrix_bc.getState( 2, 3 ) != X ) {
                return false;
            }
            if ( matrix_bc.getState( 2, 4 ) != X ) {
                return false;
            }
            final PhylogenyFactory factory0 = ParserBasedPhylogenyFactory.getInstance();
            final String p0_str = "((one,two)1-2,(three,four)3-4)root";
            final Phylogeny p0 = factory0.create( p0_str, new NHXParser() )[ 0 ];
            final DomainParsimonyCalculator dp0 = DomainParsimonyCalculator.createInstance( p0, gwcd_list );
            dp0.executeDolloParsimonyOnDomainPresence();
            final CharacterStateMatrix<GainLossStates> gl_matrix_d = dp0.getGainLossMatrix();
            final CharacterStateMatrix<BinaryStates> is_matrix_d = dp0.getInternalStatesMatrix();
            dp0.executeDolloParsimonyOnBinaryDomainCombintionPresence();
            final CharacterStateMatrix<GainLossStates> gl_matrix_bc = dp0.getGainLossMatrix();
            final CharacterStateMatrix<BinaryStates> is_matrix_bc = dp0.getInternalStatesMatrix();
            if ( is_matrix_d.getState( "root", "A" ) != X ) {
                return false;
            }
            if ( is_matrix_d.getState( "root", "B" ) != X ) {
                return false;
            }
            if ( is_matrix_d.getState( "root", "C" ) != O ) {
                return false;
            }
            if ( is_matrix_d.getState( "root", "D" ) != O ) {
                return false;
            }
            if ( is_matrix_d.getState( "root", "E" ) != X ) {
                return false;
            }
            if ( is_matrix_bc.getState( "root", "A=A" ) != O ) {
                return false;
            }
            if ( is_matrix_bc.getState( "root", "A=B" ) != X ) {
                return false;
            }
            if ( is_matrix_bc.getState( "root", "A=C" ) != O ) {
                return false;
            }
            if ( is_matrix_bc.getState( "root", "A=D" ) != O ) {
                return false;
            }
            if ( is_matrix_bc.getState( "root", "G=H" ) != O ) {
                return false;
            }
            if ( is_matrix_bc.getState( "1-2", "G=H" ) != O ) {
                return false;
            }
            if ( is_matrix_bc.getState( "root", "E=F" ) != X ) {
                return false;
            }
            if ( gl_matrix_bc.getState( "root", "E=F" ) != P ) {
                return false;
            }
            if ( gl_matrix_bc.getState( "root", "A=A" ) != A ) {
                return false;
            }
            if ( gl_matrix_bc.getState( "one", "A=A" ) != G ) {
                return false;
            }
            if ( gl_matrix_bc.getState( "root", "A=B" ) != P ) {
                return false;
            }
            if ( gl_matrix_bc.getState( "3-4", "A=D" ) != G ) {
                return false;
            }
            if ( gl_matrix_bc.getState( "four", "E=F" ) != L ) {
                return false;
            }
            if ( gl_matrix_d.getState( "3-4", "P" ) != G ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testParsimonyOnSecondaryFeatures() {
        try {
            final BinaryStates X = BinaryStates.PRESENT;
            final BinaryStates O = BinaryStates.ABSENT;
            final GainLossStates G = GainLossStates.GAIN;
            final GainLossStates L = GainLossStates.LOSS;
            final GainLossStates A = GainLossStates.UNCHANGED_ABSENT;
            final GainLossStates P = GainLossStates.UNCHANGED_PRESENT;
            final Domain a = new BasicDomain( "A", 1, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain b = new BasicDomain( "B", 2, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain c = new BasicDomain( "C", 3, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain d = new BasicDomain( "D", 1, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain e = new BasicDomain( "E", 2, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain f = new BasicDomain( "F", 3, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain g = new BasicDomain( "G", 1, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain h = new BasicDomain( "H", 2, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain i = new BasicDomain( "I", 3, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain j = new BasicDomain( "J", 1, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain l = new BasicDomain( "L", 3, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain m = new BasicDomain( "M", 1, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain n = new BasicDomain( "N", 2, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain o = new BasicDomain( "O", 3, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain p = new BasicDomain( "P", 1, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain q = new BasicDomain( "Q", 2, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            final Domain r = new BasicDomain( "R", 3, 25, ( short ) 1, ( short ) 4, 0.1, -12 );
            // 1 a-a a-b a-c e-f-g-h l-m
            // 2 a-b a-c e-f-g-i n-o
            // 3 a-b a-d e-f-g-j p-q
            // 4 a-b a-d p-r
            // 1 a-a a-b a-c e-f e-g e-h f-g f-h g-h l-m
            // 2 a-b a-c e-f e-g e-i f-g f-i g-i n-o
            // 3 a-b a-d e-f e-g e-j f-g f-j g-j p-q
            // 4 a-b a-d p-r
            // 1 a b c e f g h l m
            // 2 a b c e f g i n o
            // 3 a b d e f g j p q
            // 4 a b d p r
            final Protein aa1 = new BasicProtein( "aa1", "one" );
            aa1.addProteinDomain( a );
            aa1.addProteinDomain( a );
            final Protein ab1 = new BasicProtein( "ab1", "one" );
            ab1.addProteinDomain( a );
            ab1.addProteinDomain( b );
            final Protein ac1 = new BasicProtein( "ac1", "one" );
            ac1.addProteinDomain( a );
            ac1.addProteinDomain( c );
            final Protein efgh1 = new BasicProtein( "efgh1", "one" );
            efgh1.addProteinDomain( e );
            efgh1.addProteinDomain( f );
            efgh1.addProteinDomain( g );
            efgh1.addProteinDomain( h );
            final Protein lm1 = new BasicProtein( "lm1", "one" );
            lm1.addProteinDomain( l );
            lm1.addProteinDomain( m );
            final Protein ab2 = new BasicProtein( "ab2", "two" );
            ab2.addProteinDomain( a );
            ab2.addProteinDomain( b );
            final Protein ac2 = new BasicProtein( "ac2", "two" );
            ac2.addProteinDomain( a );
            ac2.addProteinDomain( c );
            final Protein efgi2 = new BasicProtein( "efgi2", "two" );
            efgi2.addProteinDomain( e );
            efgi2.addProteinDomain( f );
            efgi2.addProteinDomain( g );
            efgi2.addProteinDomain( i );
            final Protein no2 = new BasicProtein( "no2", "two" );
            no2.addProteinDomain( n );
            no2.addProteinDomain( o );
            final Protein ab3 = new BasicProtein( "ab3", "three" );
            ab3.addProteinDomain( a );
            ab3.addProteinDomain( b );
            final Protein ad3 = new BasicProtein( "ad3", "three" );
            ad3.addProteinDomain( a );
            ad3.addProteinDomain( d );
            final Protein efgj3 = new BasicProtein( "efgj3", "three" );
            efgj3.addProteinDomain( e );
            efgj3.addProteinDomain( f );
            efgj3.addProteinDomain( g );
            efgj3.addProteinDomain( j );
            final Protein pq3 = new BasicProtein( "pq3", "three" );
            pq3.addProteinDomain( p );
            pq3.addProteinDomain( q );
            final Protein ab4 = new BasicProtein( "ab4", "four" );
            ab4.addProteinDomain( a );
            ab4.addProteinDomain( b );
            final Protein ad4 = new BasicProtein( "ad4", "four" );
            ad4.addProteinDomain( a );
            ad4.addProteinDomain( d );
            final Protein pr4 = new BasicProtein( "pr4", "four" );
            pr4.addProteinDomain( p );
            pr4.addProteinDomain( r );
            final List<Protein> one_list = new ArrayList<Protein>();
            one_list.add( aa1 );
            one_list.add( ab1 );
            one_list.add( ac1 );
            one_list.add( efgh1 );
            one_list.add( lm1 );
            final List<Protein> two_list = new ArrayList<Protein>();
            two_list.add( ab2 );
            two_list.add( ac2 );
            two_list.add( efgi2 );
            two_list.add( no2 );
            final List<Protein> three_list = new ArrayList<Protein>();
            three_list.add( ab3 );
            three_list.add( ad3 );
            three_list.add( efgj3 );
            three_list.add( pq3 );
            final List<Protein> four_list = new ArrayList<Protein>();
            four_list.add( ab4 );
            four_list.add( ad4 );
            four_list.add( pr4 );
            final GenomeWideCombinableDomains one = BasicGenomeWideCombinableDomains
                    .createInstance( one_list, false, new BasicSpecies( "one" ) );
            final GenomeWideCombinableDomains two = BasicGenomeWideCombinableDomains
                    .createInstance( two_list, false, new BasicSpecies( "two" ) );
            final GenomeWideCombinableDomains three = BasicGenomeWideCombinableDomains
                    .createInstance( three_list, false, new BasicSpecies( "three" ) );
            final GenomeWideCombinableDomains four = BasicGenomeWideCombinableDomains
                    .createInstance( four_list, false, new BasicSpecies( "four" ) );
            final List<GenomeWideCombinableDomains> gwcd_list = new ArrayList<GenomeWideCombinableDomains>();
            gwcd_list.add( one );
            gwcd_list.add( two );
            gwcd_list.add( three );
            gwcd_list.add( four );
            final Map<DomainId, Set<String>> map_same = new HashMap<DomainId, Set<String>>();
            final HashSet<String> a_s = new HashSet<String>();
            a_s.add( "AAA" );
            final HashSet<String> b_s = new HashSet<String>();
            b_s.add( "BBB" );
            final HashSet<String> c_s = new HashSet<String>();
            c_s.add( "CCC" );
            final HashSet<String> d_s = new HashSet<String>();
            d_s.add( "DDD" );
            final HashSet<String> e_s = new HashSet<String>();
            e_s.add( "EEE" );
            final HashSet<String> f_s = new HashSet<String>();
            f_s.add( "FFF" );
            final HashSet<String> g_s = new HashSet<String>();
            g_s.add( "GGG" );
            final HashSet<String> h_s = new HashSet<String>();
            h_s.add( "HHH" );
            final HashSet<String> i_s = new HashSet<String>();
            i_s.add( "III" );
            final HashSet<String> j_s = new HashSet<String>();
            j_s.add( "JJJ" );
            final HashSet<String> l_s = new HashSet<String>();
            l_s.add( "LLL" );
            final HashSet<String> m_s = new HashSet<String>();
            m_s.add( "MMM" );
            final HashSet<String> n_s = new HashSet<String>();
            n_s.add( "NNN" );
            final HashSet<String> o_s = new HashSet<String>();
            o_s.add( "OOO" );
            final HashSet<String> p_s = new HashSet<String>();
            p_s.add( "PPP" );
            final HashSet<String> q_s = new HashSet<String>();
            q_s.add( "QQQ" );
            final HashSet<String> r_s = new HashSet<String>();
            r_s.add( "RRR" );
            map_same.put( a.getDomainId(), a_s );
            map_same.put( b.getDomainId(), b_s );
            map_same.put( c.getDomainId(), c_s );
            map_same.put( d.getDomainId(), d_s );
            map_same.put( e.getDomainId(), e_s );
            map_same.put( f.getDomainId(), f_s );
            map_same.put( g.getDomainId(), g_s );
            map_same.put( h.getDomainId(), h_s );
            map_same.put( i.getDomainId(), i_s );
            map_same.put( j.getDomainId(), j_s );
            map_same.put( l.getDomainId(), l_s );
            map_same.put( m.getDomainId(), m_s );
            map_same.put( n.getDomainId(), n_s );
            map_same.put( o.getDomainId(), o_s );
            map_same.put( p.getDomainId(), p_s );
            map_same.put( q.getDomainId(), q_s );
            map_same.put( r.getDomainId(), r_s );
            final CharacterStateMatrix<BinaryStates> matrix_s = DomainParsimonyCalculator
                    .createMatrixOfSecondaryFeaturePresenceOrAbsence( gwcd_list, map_same, null );
            // 1 a b c e f g h l m
            // 2 a b c e f g i n o
            // 3 a b d e f g j p q
            // 4 a b d p r
            if ( matrix_s.getState( 0, 0 ) != X ) {
                return false;
            }
            if ( matrix_s.getState( 0, 1 ) != X ) {
                return false;
            }
            if ( matrix_s.getState( 0, 2 ) != X ) {
                return false;
            }
            if ( matrix_s.getState( 0, 3 ) != O ) {
                return false;
            }
            if ( matrix_s.getState( 0, 4 ) != X ) {
                return false;
            }
            if ( matrix_s.getState( 0, 5 ) != X ) {
                return false;
            }
            if ( matrix_s.getState( 0, 6 ) != X ) {
                return false;
            }
            if ( matrix_s.getState( 0, 7 ) != X ) {
                return false;
            }
            if ( matrix_s.getState( 0, 8 ) != O ) {
                return false;
            }
            final PhylogenyFactory factory0 = ParserBasedPhylogenyFactory.getInstance();
            final String p0_str = "((one,two)1-2,(three,four)3-4)root";
            final Phylogeny p0 = factory0.create( p0_str, new NHXParser() )[ 0 ];
            final DomainParsimonyCalculator dp0 = DomainParsimonyCalculator.createInstance( p0, gwcd_list, map_same );
            dp0.executeDolloParsimonyOnSecondaryFeatures( null );
            final CharacterStateMatrix<GainLossStates> gl_matrix_d = dp0.getGainLossMatrix();
            final CharacterStateMatrix<BinaryStates> is_matrix_d = dp0.getInternalStatesMatrix();
            if ( is_matrix_d.getState( "root", "AAA" ) != X ) {
                return false;
            }
            if ( is_matrix_d.getState( "root", "BBB" ) != X ) {
                return false;
            }
            if ( is_matrix_d.getState( "root", "CCC" ) != O ) {
                return false;
            }
            if ( is_matrix_d.getState( "root", "DDD" ) != O ) {
                return false;
            }
            if ( is_matrix_d.getState( "root", "EEE" ) != X ) {
                return false;
            }
            if ( gl_matrix_d.getState( "3-4", "PPP" ) != G ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testPaupLogParser( final File test_dir ) {
        try {
            final PaupLogParser parser = new PaupLogParser();
            parser.setSource( new File( test_dir + ForesterUtil.getFileSeparator() + "paup_log_test_1" ) );
            final CharacterStateMatrix<BinaryStates> matrix = parser.parse();
            if ( matrix.getNumberOfIdentifiers() != 8 ) {
                return false;
            }
            if ( !matrix.getIdentifier( 0 ).equals( "MOUSE" ) ) {
                return false;
            }
            if ( !matrix.getIdentifier( 1 ).equals( "NEMVE" ) ) {
                return false;
            }
            if ( !matrix.getIdentifier( 2 ).equals( "MONBE" ) ) {
                return false;
            }
            if ( !matrix.getIdentifier( 3 ).equals( "DICDI" ) ) {
                return false;
            }
            if ( !matrix.getIdentifier( 4 ).equals( "ARATH" ) ) {
                return false;
            }
            if ( !matrix.getIdentifier( 5 ).equals( "6" ) ) {
                return false;
            }
            if ( !matrix.getIdentifier( 6 ).equals( "7" ) ) {
                return false;
            }
            if ( !matrix.getIdentifier( 7 ).equals( "8" ) ) {
                return false;
            }
            if ( matrix.getNumberOfCharacters() != ( 66 + 66 + 28 ) ) {
                return false;
            }
            if ( matrix.getState( 0, 4 ) != BinaryStates.ABSENT ) {
                return false;
            }
            if ( matrix.getState( 0, 5 ) != BinaryStates.PRESENT ) {
                return false;
            }
            if ( matrix.getState( 1, 5 ) != BinaryStates.PRESENT ) {
                return false;
            }
            if ( matrix.getState( 7, 154 ) != BinaryStates.ABSENT ) {
                return false;
            }
            if ( matrix.getState( 7, 155 ) != BinaryStates.PRESENT ) {
                return false;
            }
            if ( matrix.getState( 7, 156 ) != BinaryStates.PRESENT ) {
                return false;
            }
            if ( matrix.getState( 7, 157 ) != BinaryStates.ABSENT ) {
                return false;
            }
            if ( matrix.getState( 7, 158 ) != BinaryStates.PRESENT ) {
                return false;
            }
            if ( matrix.getState( 7, 159 ) != BinaryStates.ABSENT ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testProteinId() {
        try {
            final ProteinId id1 = new ProteinId( "a" );
            final ProteinId id2 = new ProteinId( "a" );
            final ProteinId id3 = new ProteinId( "A" );
            final ProteinId id4 = new ProteinId( "b" );
            if ( !id1.equals( id1 ) ) {
                return false;
            }
            if ( id1.getId().equals( "x" ) ) {
                return false;
            }
            if ( id1.getId().equals( null ) ) {
                return false;
            }
            if ( !id1.equals( id2 ) ) {
                return false;
            }
            if ( id1.equals( id3 ) ) {
                return false;
            }
            if ( id1.hashCode() != id1.hashCode() ) {
                return false;
            }
            if ( id1.hashCode() != id2.hashCode() ) {
                return false;
            }
            if ( id1.hashCode() == id3.hashCode() ) {
                return false;
            }
            if ( id1.compareTo( id1 ) != 0 ) {
                return false;
            }
            if ( id1.compareTo( id2 ) != 0 ) {
                return false;
            }
            if ( id1.compareTo( id3 ) != 0 ) {
                return false;
            }
            if ( id1.compareTo( id4 ) >= 0 ) {
                return false;
            }
            if ( id4.compareTo( id1 ) <= 0 ) {
                return false;
            }
            if ( !id4.getId().equals( "b" ) ) {
                return false;
            }
            final ProteinId id5 = new ProteinId( " C " );
            if ( !id5.getId().equals( "C" ) ) {
                return false;
            }
            if ( id5.equals( id1 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSpecies() {
        try {
            final Species s1 = new BasicSpecies( "a" );
            final Species s2 = new BasicSpecies( "a" );
            final Species s3 = new BasicSpecies( "A" );
            final Species s4 = new BasicSpecies( "b" );
            if ( !s1.equals( s1 ) ) {
                return false;
            }
            if ( s1.getSpeciesId().equals( "x" ) ) {
                return false;
            }
            if ( s1.getSpeciesId().equals( null ) ) {
                return false;
            }
            if ( !s1.equals( s2 ) ) {
                return false;
            }
            if ( s1.equals( s3 ) ) {
                return false;
            }
            if ( s1.hashCode() != s1.hashCode() ) {
                return false;
            }
            if ( s1.hashCode() != s2.hashCode() ) {
                return false;
            }
            if ( s1.hashCode() == s3.hashCode() ) {
                return false;
            }
            if ( s1.compareTo( s1 ) != 0 ) {
                return false;
            }
            if ( s1.compareTo( s2 ) != 0 ) {
                return false;
            }
            if ( s1.compareTo( s3 ) != 0 ) {
                return false;
            }
            if ( s1.compareTo( s4 ) >= 0 ) {
                return false;
            }
            if ( s4.compareTo( s1 ) <= 0 ) {
                return false;
            }
            if ( !s4.getSpeciesId().equals( "b" ) ) {
                return false;
            }
            final Species s5 = new BasicSpecies( " C " );
            if ( !s5.getSpeciesId().equals( "C" ) ) {
                return false;
            }
            if ( s5.equals( s1 ) ) {
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
