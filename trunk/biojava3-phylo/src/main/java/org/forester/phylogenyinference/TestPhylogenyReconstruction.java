// $Id: TestPhylogenyReconstruction.java,v 1.8 2007/11/30 23:45:07 cmzmasek Exp
// $
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

package org.forester.phylogenyinference;

import java.io.File;
import java.io.StringWriter;
import java.util.Date;
import java.util.List;

import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogenyinference.CharacterStateMatrix.BinaryStates;
import org.forester.phylogenyinference.CharacterStateMatrix.GainLossStates;
import org.forester.util.ForesterUtil;

public class TestPhylogenyReconstruction {

    private final static double  ZERO_DIFF = 1.0E-9;
    private final static boolean TIME      = false;

    public static boolean isEqual( final double a, final double b ) {
        return ( ( Math.abs( a - b ) ) < ZERO_DIFF );
    }

    public static void main( final String[] args ) {
        timeNeighborJoining();
    }

    public static boolean test( final File test_dir ) {
        System.out.print( "  Basic symmetrical distance matrix: " );
        if ( !testBasicSymmetricalDistanceMatrix() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Basic character state matrix: " );
        if ( !testBasicCharacterStateMatrix() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Symmetrical distance matrix parser: " );
        if ( !testSymmetricalDistanceMatrixParser() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Neighbor Joining: " );
        if ( !testNeighborJoining() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Dollo Parsimony: " );
        if ( !testDolloParsimony() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Dollo Parsimony on non binary trees: " );
        if ( !testDolloParsimonyOnNonBinaryTree() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Fitch Parsimony: " );
        if ( !testFitchParsimony() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        return true;
    }

    private static boolean testBasicCharacterStateMatrix() {
        try {
            final CharacterStateMatrix<String> matrix_0 = new BasicCharacterStateMatrix<String>( 4, 8 );
            final CharacterStateMatrix<String> matrix_00 = new BasicCharacterStateMatrix<String>( 4, 8 );
            matrix_0.setIdentifier( 0, "A" );
            matrix_0.setIdentifier( 1, "B" );
            matrix_0.setIdentifier( 2, "C" );
            matrix_0.setIdentifier( 3, "D" );
            matrix_0.setCharacter( 0, "0" );
            matrix_0.setCharacter( 1, "1" );
            matrix_0.setCharacter( 2, "2" );
            matrix_0.setCharacter( 3, "3" );
            matrix_0.setCharacter( 4, "4" );
            matrix_0.setCharacter( 5, "5" );
            matrix_0.setCharacter( 6, "6" );
            matrix_0.setCharacter( 7, "7" );
            matrix_00.setIdentifier( 0, "A" );
            matrix_00.setIdentifier( 1, "B" );
            matrix_00.setIdentifier( 2, "C" );
            matrix_00.setIdentifier( 3, "D" );
            matrix_00.setCharacter( 3, "3" );
            matrix_00.setCharacter( 4, "4" );
            if ( !matrix_0.getCharacter( 1 ).equals( "1" ) ) {
                return false;
            }
            if ( !matrix_0.getIdentifier( 0 ).equals( "A" ) ) {
                return false;
            }
            matrix_0.setState( 0, 0, "00" );
            matrix_00.setState( 0, 0, "00" );
            if ( !matrix_0.getState( 0, 0 ).equals( "00" ) ) {
                return false;
            }
            matrix_0.setState( 0, 1, "01" );
            matrix_00.setState( 0, 1, "01" );
            if ( !matrix_0.getState( 0, 1 ).equals( "01" ) ) {
                return false;
            }
            matrix_0.setState( 1, 1, "11" );
            matrix_00.setState( 1, 1, "11" );
            if ( !matrix_0.getState( 1, 1 ).equals( "11" ) ) {
                return false;
            }
            matrix_0.setState( 1, 0, "10" );
            matrix_00.setState( 1, 0, "10" );
            if ( !matrix_0.getState( 1, 0 ).equals( "10" ) ) {
                return false;
            }
            matrix_0.setState( 1, 2, "12" );
            matrix_00.setState( 1, 2, "12" );
            if ( !matrix_0.getState( 1, 2 ).equals( "12" ) ) {
                return false;
            }
            matrix_0.setState( 3, 7, "37" );
            matrix_00.setState( 3, 7, "37" );
            if ( !matrix_0.getState( 3, 7 ).equals( "37" ) ) {
                return false;
            }
            matrix_0.setState( 2, 6, "26" );
            matrix_00.setState( 2, 6, "26" );
            if ( !matrix_0.getState( 2, 6 ).equals( "26" ) ) {
                return false;
            }
            matrix_0.setState( "D", "3", "33" );
            matrix_00.setState( "D", "3", "33" );
            if ( !matrix_0.getState( 3, 3 ).equals( "33" ) ) {
                return false;
            }
            if ( !matrix_0.getState( "D", "3" ).equals( "33" ) ) {
                return false;
            }
            matrix_0.setState( "C", "4", "24" );
            matrix_00.setState( "C", "4", "24" );
            if ( !matrix_0.getState( 2, 4 ).equals( "24" ) ) {
                return false;
            }
            if ( !matrix_0.getState( "C", "4" ).equals( "24" ) ) {
                return false;
            }
            if ( matrix_0.isEmpty() ) {
                return false;
            }
            if ( matrix_0.getNumberOfIdentifiers() != 4 ) {
                return false;
            }
            if ( matrix_0.getNumberOfCharacters() != 8 ) {
                return false;
            }
            if ( !matrix_0.equals( matrix_0 ) ) {
                return false;
            }
            if ( !matrix_0.equals( matrix_00 ) ) {
                return false;
            }
            matrix_00.setState( "C", "4", "123" );
            if ( matrix_0.equals( matrix_00 ) ) {
                return false;
            }
            final Integer[][] ints = { { 1, 2, 3, 4 }, { 5, 6, 7, 8 }, { 9, 10, 11, 12 } };
            final CharacterStateMatrix<Integer> matrix_000 = new BasicCharacterStateMatrix<Integer>( ints );
            matrix_000.toString();
            if ( matrix_000.getNumberOfCharacters() != 4 ) {
                return false;
            }
            if ( matrix_000.getNumberOfIdentifiers() != 3 ) {
                return false;
            }
            if ( matrix_000.getState( 0, 1 ) != 2 ) {
                return false;
            }
            if ( matrix_000.getState( 2, 3 ) != 12 ) {
                return false;
            }
            final Integer[][] ints0 = { { 1, 2, 3, 4 }, { 5, 6, 7, 8 }, { 9, 10, 11, 12 } };
            final CharacterStateMatrix<Integer> matrix_0000 = new BasicCharacterStateMatrix<Integer>( ints0 );
            if ( !matrix_000.equals( matrix_0000 ) ) {
                return false;
            }
            final Integer[][] ints00 = { { 1, 2, 3, -4 }, { 5, 6, 7, 8 }, { 9, 10, 11, 12 } };
            final CharacterStateMatrix<Integer> matrix_00000 = new BasicCharacterStateMatrix<Integer>( ints00 );
            if ( matrix_000.equals( matrix_00000 ) ) {
                return false;
            }
            final CharacterStateMatrix<String> clone0 = matrix_0.copy();
            final CharacterStateMatrix<String> clone00 = matrix_00.copy();
            if ( !clone0.equals( matrix_0 ) ) {
                return false;
            }
            if ( !clone00.equals( matrix_00 ) ) {
                return false;
            }
            if ( clone00.equals( clone0 ) ) {
                return false;
            }
            final CharacterStateMatrix<String> pivot0 = matrix_0.pivot();
            final CharacterStateMatrix<String> pivot00 = matrix_00.pivot();
            if ( !pivot0.getState( 1, 0 ).equals( "01" ) ) {
                return false;
            }
            if ( !pivot0.getState( 6, 2 ).equals( "26" ) ) {
                return false;
            }
            if ( !matrix_0.getState( 2, 6 ).equals( "26" ) ) {
                return false;
            }
            final CharacterStateMatrix<String> pivotpivot00 = pivot00.pivot();
            if ( !pivotpivot00.equals( matrix_00 ) ) {
                return false;
            }
            final CharacterStateMatrix<BinaryStates> nex = new BasicCharacterStateMatrix<BinaryStates>( 4, 3 );
            nex.setIdentifier( 0, "amphioxus" );
            nex.setIdentifier( 1, "sponge" );
            nex.setIdentifier( 2, "sea_anemone" );
            nex.setIdentifier( 3, "cobra" );
            nex.setCharacter( 0, "notch" );
            nex.setCharacter( 1, "homeobox" );
            nex.setCharacter( 2, "wnt" );
            nex.setState( 0, 0, BinaryStates.ABSENT );
            nex.setState( 0, 1, BinaryStates.ABSENT );
            nex.setState( 0, 2, BinaryStates.ABSENT );
            nex.setState( 1, 0, BinaryStates.PRESENT );
            nex.setState( 1, 1, BinaryStates.PRESENT );
            nex.setState( 1, 2, BinaryStates.ABSENT );
            nex.setState( 2, 0, BinaryStates.PRESENT );
            nex.setState( 2, 1, BinaryStates.PRESENT );
            nex.setState( 2, 2, BinaryStates.PRESENT );
            nex.setState( 3, 0, BinaryStates.PRESENT );
            nex.setState( 3, 1, BinaryStates.ABSENT );
            nex.setState( 3, 2, BinaryStates.ABSENT );
            StringWriter w = new StringWriter();
            nex.toWriter( w, CharacterStateMatrix.Format.NEXUS_BINARY );
            //System.out.println( w.getBuffer().toString() );
            w = new StringWriter();
            nex.pivot().toWriter( w, CharacterStateMatrix.Format.NEXUS_BINARY );
            //System.out.println( w.getBuffer().toString() );
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testBasicSymmetricalDistanceMatrix() {
        try {
            final DistanceMatrix matrix_0 = new BasicSymmetricalDistanceMatrix( 4 );
            matrix_0.setIdentifier( 0, "A" );
            matrix_0.setIdentifier( 1, "B" );
            matrix_0.setIdentifier( 2, "C" );
            matrix_0.setIdentifier( 3, "0123456789012" );
            matrix_0.setValue( 1, 0, 0.00001 );
            matrix_0.setValue( 0, 2, 0.0000009 );
            matrix_0.setValue( 3, 0, 3.0 );
            matrix_0.setValue( 1, 2, 4.0 );
            matrix_0.setValue( 3, 1, 5.0 );
            matrix_0.setValue( 2, 3, 6.0 );
            if ( !matrix_0.getIdentifier( 0 ).equals( "A" ) ) {
                return false;
            }
            if ( !matrix_0.getIdentifier( 1 ).equals( "B" ) ) {
                return false;
            }
            if ( !matrix_0.getIdentifier( 2 ).equals( "C" ) ) {
                return false;
            }
            if ( !matrix_0.getIdentifier( 3 ).equals( "0123456789012" ) ) {
                return false;
            }
            if ( matrix_0.getSize() != 4 ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 0, 0 ), 0.0 ) ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 3, 3 ), 0.0 ) ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 0, 1 ), 0.00001 ) ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 0, 2 ), 0.0000009 ) ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 0, 3 ), 3 ) ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 1, 0 ), 0.00001 ) ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 1, 2 ), 4 ) ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 1, 3 ), 5 ) ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 2, 0 ), 0.0000009 ) ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 2, 1 ), 4 ) ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 2, 3 ), 6 ) ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 3, 0 ), 3 ) ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 3, 1 ), 5 ) ) {
                return false;
            }
            if ( !isEqual( matrix_0.getValue( 3, 2 ), 6 ) ) {
                return false;
            }
            final StringBuffer matrix_0_phylip = new StringBuffer();
            matrix_0_phylip.append( "    4" );
            matrix_0_phylip.append( ForesterUtil.LINE_SEPARATOR );
            matrix_0_phylip.append( "A           0.000000  0.000010  0.000001  3.000000" );
            matrix_0_phylip.append( ForesterUtil.LINE_SEPARATOR );
            matrix_0_phylip.append( "B           0.000010  0.000000  4.000000  5.000000" );
            matrix_0_phylip.append( ForesterUtil.LINE_SEPARATOR );
            matrix_0_phylip.append( "C           0.000001  4.000000  0.000000  6.000000" );
            matrix_0_phylip.append( ForesterUtil.LINE_SEPARATOR );
            matrix_0_phylip.append( "0123456789  3.000000  5.000000  6.000000  0.000000" );
            if ( !matrix_0_phylip.toString()
                    .equals( matrix_0.toStringBuffer( DistanceMatrix.Format.PHYLIP ).toString() ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDolloParsimony() {
        try {
            final BinaryStates PRESENT = BinaryStates.PRESENT;
            final BinaryStates ABSENT = BinaryStates.ABSENT;
            final GainLossStates UNCHANGED_PRESENT = GainLossStates.UNCHANGED_PRESENT;
            final DolloParsimony dollo1 = DolloParsimony.createInstance();
            final PhylogenyFactory factory1 = ParserBasedPhylogenyFactory.getInstance();
            final String p1_str = "((((((a,b)ab,c)ac,d)ad,(e,f)ef)af,(g,h)gh)ah,i)r";
            final Phylogeny p1 = factory1.create( p1_str, new NHXParser() )[ 0 ];
            CharacterStateMatrix<CharacterStateMatrix.BinaryStates> m1 = new BasicCharacterStateMatrix<CharacterStateMatrix.BinaryStates>( 9,
                                                                                                                                           1 );
            m1.setIdentifier( 0, "a" );
            m1.setIdentifier( 1, "b" );
            m1.setIdentifier( 2, "c" );
            m1.setIdentifier( 3, "d" );
            m1.setIdentifier( 4, "e" );
            m1.setIdentifier( 5, "f" );
            m1.setIdentifier( 6, "g" );
            m1.setIdentifier( 7, "h" );
            m1.setIdentifier( 8, "i" );
            m1.setCharacter( 0, "0" );
            m1.setState( "a", "0", PRESENT );
            m1.setState( "b", "0", ABSENT );
            m1.setState( "c", "0", PRESENT );
            m1.setState( "d", "0", ABSENT );
            m1.setState( "e", "0", ABSENT );
            m1.setState( "f", "0", ABSENT );
            m1.setState( "g", "0", ABSENT );
            m1.setState( "h", "0", ABSENT );
            m1.setState( "i", "0", ABSENT );
            dollo1.execute( p1, m1 );
            if ( dollo1.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo1.getTotalLosses() != 1 ) {
                return false;
            }
            if ( dollo1.getTotalUnchanged() != 15 ) {
                return false;
            }
            m1.setState( "b", "0", PRESENT );
            dollo1.execute( p1, m1 );
            if ( dollo1.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo1.getTotalLosses() != 0 ) {
                return false;
            }
            if ( dollo1.getTotalUnchanged() != 16 ) {
                return false;
            }
            m1.setState( "b", "0", ABSENT );
            m1.setState( "e", "0", PRESENT );
            dollo1.execute( p1, m1 );
            if ( dollo1.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo1.getTotalLosses() != 3 ) {
                return false;
            }
            if ( dollo1.getTotalUnchanged() != 13 ) {
                return false;
            }
            m1.setState( "a", "0", ABSENT );
            m1.setState( "c", "0", ABSENT );
            m1.setState( "g", "0", PRESENT );
            dollo1.setReturnInternalStates( true );
            dollo1.setReturnGainLossMatrix( true );
            dollo1.execute( p1, m1 );
            if ( dollo1.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo1.getTotalLosses() != 3 ) {
                return false;
            }
            if ( dollo1.getTotalUnchanged() != 13 ) {
                return false;
            }
            final DolloParsimony dollo2 = DolloParsimony.createInstance();
            final PhylogenyFactory factory2 = ParserBasedPhylogenyFactory.getInstance();
            final String p2_str = "((((((a,b)ab,c)ac,d)ad,(e,f)ef)af,(g,h,i)gi)ai,((j,k,l)jl,(m,n,o)mo,(p,q,r)pr)jr)root";
            final Phylogeny p2 = factory2.create( p2_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<CharacterStateMatrix.BinaryStates> m2 = new BasicCharacterStateMatrix<CharacterStateMatrix.BinaryStates>( 18,
                                                                                                                                                 4 );
            m2.setIdentifier( 0, "a" );
            m2.setIdentifier( 1, "b" );
            m2.setIdentifier( 2, "c" );
            m2.setIdentifier( 3, "d" );
            m2.setIdentifier( 4, "e" );
            m2.setIdentifier( 5, "f" );
            m2.setIdentifier( 6, "g" );
            m2.setIdentifier( 7, "h" );
            m2.setIdentifier( 8, "i" );
            m2.setIdentifier( 9, "j" );
            m2.setIdentifier( 10, "k" );
            m2.setIdentifier( 11, "l" );
            m2.setIdentifier( 12, "m" );
            m2.setIdentifier( 13, "n" );
            m2.setIdentifier( 14, "o" );
            m2.setIdentifier( 15, "p" );
            m2.setIdentifier( 16, "q" );
            m2.setIdentifier( 17, "r" );
            m2.setCharacter( 0, "0" );
            m2.setCharacter( 1, "1" );
            m2.setCharacter( 2, "2" );
            m2.setCharacter( 3, "3" );
            m2.setState( "a", "0", PRESENT );
            m2.setState( "b", "0", ABSENT );
            m2.setState( "c", "0", PRESENT );
            m2.setState( "d", "0", ABSENT );
            m2.setState( "e", "0", ABSENT );
            m2.setState( "f", "0", ABSENT );
            m2.setState( "g", "0", ABSENT );
            m2.setState( "h", "0", ABSENT );
            m2.setState( "i", "0", ABSENT );
            m2.setState( "j", "0", ABSENT );
            m2.setState( "k", "0", ABSENT );
            m2.setState( "l", "0", ABSENT );
            m2.setState( "m", "0", ABSENT );
            m2.setState( "n", "0", ABSENT );
            m2.setState( "o", "0", ABSENT );
            m2.setState( "p", "0", ABSENT );
            m2.setState( "q", "0", ABSENT );
            m2.setState( "r", "0", ABSENT );
            m2.setState( "a", "1", PRESENT );
            m2.setState( "b", "1", ABSENT );
            m2.setState( "c", "1", PRESENT );
            m2.setState( "d", "1", ABSENT );
            m2.setState( "e", "1", ABSENT );
            m2.setState( "f", "1", ABSENT );
            m2.setState( "g", "1", PRESENT );
            m2.setState( "h", "1", ABSENT );
            m2.setState( "i", "1", ABSENT );
            m2.setState( "j", "1", PRESENT );
            m2.setState( "k", "1", ABSENT );
            m2.setState( "l", "1", ABSENT );
            m2.setState( "m", "1", PRESENT );
            m2.setState( "n", "1", ABSENT );
            m2.setState( "o", "1", ABSENT );
            m2.setState( "p", "1", ABSENT );
            m2.setState( "q", "1", ABSENT );
            m2.setState( "r", "1", ABSENT );
            m2.setState( "a", "2", ABSENT );
            m2.setState( "b", "2", ABSENT );
            m2.setState( "c", "2", ABSENT );
            m2.setState( "d", "2", ABSENT );
            m2.setState( "e", "2", ABSENT );
            m2.setState( "f", "2", ABSENT );
            m2.setState( "g", "2", ABSENT );
            m2.setState( "h", "2", ABSENT );
            m2.setState( "i", "2", ABSENT );
            m2.setState( "j", "2", PRESENT );
            m2.setState( "k", "2", ABSENT );
            m2.setState( "l", "2", ABSENT );
            m2.setState( "m", "2", PRESENT );
            m2.setState( "n", "2", ABSENT );
            m2.setState( "o", "2", ABSENT );
            m2.setState( "p", "2", PRESENT );
            m2.setState( "q", "2", ABSENT );
            m2.setState( "r", "2", ABSENT );
            m2.setState( "a", "3", ABSENT );
            m2.setState( "b", "3", ABSENT );
            m2.setState( "c", "3", PRESENT );
            m2.setState( "d", "3", ABSENT );
            m2.setState( "e", "3", ABSENT );
            m2.setState( "f", "3", ABSENT );
            m2.setState( "g", "3", PRESENT );
            m2.setState( "h", "3", ABSENT );
            m2.setState( "i", "3", ABSENT );
            m2.setState( "j", "3", ABSENT );
            m2.setState( "k", "3", ABSENT );
            m2.setState( "l", "3", ABSENT );
            m2.setState( "m", "3", ABSENT );
            m2.setState( "n", "3", ABSENT );
            m2.setState( "o", "3", ABSENT );
            m2.setState( "p", "3", ABSENT );
            m2.setState( "q", "3", ABSENT );
            m2.setState( "r", "3", ABSENT );
            dollo2.setReturnInternalStates( true );
            dollo2.setReturnGainLossMatrix( true );
            dollo2.execute( p2, m2 );
            final CharacterStateMatrix<BinaryStates> i_m = dollo2.getInternalStatesMatrix();
            final CharacterStateMatrix<GainLossStates> gl_m = dollo2.getGainLossMatrix();
            if ( dollo2.getTotalGains() != 3 ) {
                return false;
            }
            if ( dollo2.getTotalLosses() != 22 ) {
                return false;
            }
            if ( dollo2.getTotalUnchanged() != 95 ) {
                return false;
            }
            if ( i_m.getState( "ab", "0" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "ac", "0" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "ad", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "af", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "ef", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "ai", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "gi", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "jl", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "mo", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "pr", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "jr", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "root", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "ab", "1" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "ac", "1" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "ad", "1" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "af", "1" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "ef", "1" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "ai", "1" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "gi", "1" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "jl", "1" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "mo", "1" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "pr", "1" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "jr", "1" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "root", "1" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "ab", "2" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "ac", "2" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "ad", "2" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "af", "2" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "ef", "2" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "ai", "2" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "gi", "2" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "jl", "2" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "mo", "2" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "pr", "2" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "jr", "2" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "root", "2" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "ab", "3" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "ac", "3" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "ad", "3" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "af", "3" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "ef", "3" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "ai", "3" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "gi", "3" ) != PRESENT ) {
                return false;
            }
            if ( i_m.getState( "jl", "3" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "mo", "3" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "pr", "3" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "jr", "3" ) != ABSENT ) {
                return false;
            }
            if ( i_m.getState( "root", "3" ) != ABSENT ) {
                return false;
            }
            if ( gl_m.getState( "a", "0" ) != UNCHANGED_PRESENT ) {
                return false;
            }
            final DolloParsimony dollo9 = DolloParsimony.createInstance();
            final PhylogenyFactory factory9 = ParserBasedPhylogenyFactory.getInstance();
            final String p9_str = "((((((a,b)ab,c)ac,d)ad,(e,f)ef)af,(g,h)gh)ah,i)r";
            final Phylogeny p9 = factory9.create( p9_str, new NHXParser() )[ 0 ];
            m1 = new BasicCharacterStateMatrix<CharacterStateMatrix.BinaryStates>( 9, 3 );
            m1.setIdentifier( 0, "a" );
            m1.setIdentifier( 1, "b" );
            m1.setIdentifier( 2, "c" );
            m1.setIdentifier( 3, "d" );
            m1.setIdentifier( 4, "e" );
            m1.setIdentifier( 5, "f" );
            m1.setIdentifier( 6, "g" );
            m1.setIdentifier( 7, "h" );
            m1.setIdentifier( 8, "i" );
            m1.setState( 0, 0, PRESENT );
            m1.setState( 1, 0, ABSENT );
            m1.setState( 2, 0, PRESENT );
            m1.setState( 3, 0, ABSENT );
            m1.setState( 4, 0, ABSENT );
            m1.setState( 5, 0, ABSENT );
            m1.setState( 6, 0, ABSENT );
            m1.setState( 7, 0, ABSENT );
            m1.setState( 8, 0, ABSENT );
            m1.setState( 0, 1, PRESENT );
            m1.setState( 1, 1, PRESENT );
            m1.setState( 2, 1, PRESENT );
            m1.setState( 3, 1, PRESENT );
            m1.setState( 4, 1, ABSENT );
            m1.setState( 5, 1, ABSENT );
            m1.setState( 6, 1, ABSENT );
            m1.setState( 7, 1, ABSENT );
            m1.setState( 8, 1, ABSENT );
            m1.setState( 0, 2, PRESENT );
            m1.setState( 1, 2, ABSENT );
            m1.setState( 2, 2, ABSENT );
            m1.setState( 3, 2, ABSENT );
            m1.setState( 4, 2, ABSENT );
            m1.setState( 5, 2, ABSENT );
            m1.setState( 6, 2, ABSENT );
            m1.setState( 7, 2, PRESENT );
            m1.setState( 8, 2, ABSENT );
            dollo9.execute( p9, m1 );
            if ( dollo9.getTotalGains() != 3 ) {
                return false;
            }
            if ( dollo9.getTotalLosses() != 6 ) {
                return false;
            }
            final DolloParsimony dollo10 = DolloParsimony.createInstance();
            final PhylogenyFactory factory10 = ParserBasedPhylogenyFactory.getInstance();
            final String p10_str = "((((((a,b)ab,c)ac,d)ad,(e,f)ef)af,(g,h)gh)ah,i)r";
            final Phylogeny p10 = factory10.create( p10_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<CharacterStateMatrix.BinaryStates> m10 = new BasicCharacterStateMatrix<CharacterStateMatrix.BinaryStates>( 9,
                                                                                                                                                  1 );
            m10.setIdentifier( 0, "a" );
            m10.setIdentifier( 1, "b" );
            m10.setIdentifier( 2, "c" );
            m10.setIdentifier( 3, "d" );
            m10.setIdentifier( 4, "e" );
            m10.setIdentifier( 5, "f" );
            m10.setIdentifier( 6, "g" );
            m10.setIdentifier( 7, "h" );
            m10.setIdentifier( 8, "i" );
            m10.setState( 0, 0, PRESENT );
            m10.setState( 1, 0, ABSENT );
            m10.setState( 2, 0, PRESENT );
            m10.setState( 3, 0, ABSENT );
            m10.setState( 4, 0, ABSENT );
            m10.setState( 5, 0, ABSENT );
            m10.setState( 6, 0, ABSENT );
            m10.setState( 7, 0, ABSENT );
            m10.setState( 8, 0, ABSENT );
            dollo10.execute( p10, m10 );
            if ( dollo10.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo10.getTotalLosses() != 1 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testDolloParsimonyOnNonBinaryTree() {
        try {
            final BinaryStates PRESENT = BinaryStates.PRESENT;
            final BinaryStates ABSENT = BinaryStates.ABSENT;
            final DolloParsimony dollo1 = DolloParsimony.createInstance();
            final PhylogenyFactory factory1 = ParserBasedPhylogenyFactory.getInstance();
            final String p1_str = "((((((a,b,y)aby,c)ac,d)ad,(e,f)ef)af,(g,h)gh)ah,i)r";
            final Phylogeny p1 = factory1.create( p1_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<CharacterStateMatrix.BinaryStates> m1 = new BasicCharacterStateMatrix<CharacterStateMatrix.BinaryStates>( 10,
                                                                                                                                                 1 );
            m1.setIdentifier( 0, "a" );
            m1.setIdentifier( 1, "b" );
            m1.setIdentifier( 2, "y" );
            m1.setIdentifier( 3, "c" );
            m1.setIdentifier( 4, "d" );
            m1.setIdentifier( 5, "e" );
            m1.setIdentifier( 6, "f" );
            m1.setIdentifier( 7, "g" );
            m1.setIdentifier( 8, "h" );
            m1.setIdentifier( 9, "i" );
            m1.setCharacter( 0, "0" );
            m1.setState( "a", "0", PRESENT );
            m1.setState( "b", "0", ABSENT );
            m1.setState( "y", "0", PRESENT );
            m1.setState( "c", "0", PRESENT );
            m1.setState( "d", "0", ABSENT );
            m1.setState( "e", "0", ABSENT );
            m1.setState( "f", "0", ABSENT );
            m1.setState( "g", "0", ABSENT );
            m1.setState( "h", "0", ABSENT );
            m1.setState( "i", "0", ABSENT );
            dollo1.execute( p1, m1 );
            if ( dollo1.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo1.getTotalLosses() != 1 ) {
                return false;
            }
            if ( dollo1.getTotalUnchanged() != 16 ) {
                return false;
            }
            m1.setState( "b", "0", PRESENT );
            dollo1.execute( p1, m1 );
            if ( dollo1.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo1.getTotalLosses() != 0 ) {
                return false;
            }
            if ( dollo1.getTotalUnchanged() != 17 ) {
                return false;
            }
            m1.setState( "a", "0", ABSENT );
            m1.setState( "b", "0", ABSENT );
            dollo1.execute( p1, m1 );
            if ( dollo1.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo1.getTotalLosses() != 2 ) {
                return false;
            }
            if ( dollo1.getTotalUnchanged() != 15 ) {
                return false;
            }
            m1.setState( "y", "0", ABSENT );
            dollo1.execute( p1, m1 );
            if ( dollo1.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo1.getTotalLosses() != 0 ) {
                return false;
            }
            if ( dollo1.getTotalUnchanged() != 17 ) {
                return false;
            }
            final DolloParsimony dollo2 = DolloParsimony.createInstance();
            final PhylogenyFactory factory2 = ParserBasedPhylogenyFactory.getInstance();
            final String p2_str = "((((((a,b,y)aby,c,d)cad,e,f)af,(g,h)gh)ah,i))r";
            final Phylogeny p2 = factory2.create( p2_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<CharacterStateMatrix.BinaryStates> m2 = new BasicCharacterStateMatrix<CharacterStateMatrix.BinaryStates>( 10,
                                                                                                                                                 1 );
            m2.setIdentifier( 0, "a" );
            m2.setIdentifier( 1, "b" );
            m2.setIdentifier( 2, "y" );
            m2.setIdentifier( 3, "c" );
            m2.setIdentifier( 4, "d" );
            m2.setIdentifier( 5, "e" );
            m2.setIdentifier( 6, "f" );
            m2.setIdentifier( 7, "g" );
            m2.setIdentifier( 8, "h" );
            m2.setIdentifier( 9, "i" );
            m2.setCharacter( 0, "0" );
            m2.setState( "a", "0", PRESENT );
            m2.setState( "b", "0", ABSENT );
            m2.setState( "y", "0", PRESENT );
            m2.setState( "c", "0", PRESENT );
            m2.setState( "d", "0", ABSENT );
            m2.setState( "e", "0", ABSENT );
            m2.setState( "f", "0", ABSENT );
            m2.setState( "g", "0", ABSENT );
            m2.setState( "h", "0", ABSENT );
            m2.setState( "i", "0", ABSENT );
            dollo2.setReturnInternalStates( true );
            dollo2.execute( p2, m2 );
            CharacterStateMatrix<BinaryStates> i_m2 = dollo2.getInternalStatesMatrix();
            if ( i_m2.getState( "aby", "0" ) != PRESENT ) {
                return false;
            }
            if ( i_m2.getState( "cad", "0" ) != PRESENT ) {
                return false;
            }
            if ( i_m2.getState( "af", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m2.getState( "gh", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m2.getState( "ah", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m2.getState( "r", "0" ) != ABSENT ) {
                return false;
            }
            if ( dollo2.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo2.getTotalLosses() != 2 ) {
                return false;
            }
            if ( dollo2.getTotalUnchanged() != 14 ) {
                return false;
            }
            m2.setState( "b", "0", PRESENT );
            dollo2.execute( p2, m2 );
            if ( dollo2.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo2.getTotalLosses() != 1 ) {
                return false;
            }
            if ( dollo2.getTotalUnchanged() != 15 ) {
                return false;
            }
            m2.setState( "a", "0", ABSENT );
            m2.setState( "b", "0", ABSENT );
            dollo2.execute( p2, m2 );
            if ( dollo2.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo2.getTotalLosses() != 3 ) {
                return false;
            }
            if ( dollo2.getTotalUnchanged() != 13 ) {
                return false;
            }
            m2.setState( "y", "0", ABSENT );
            dollo2.execute( p2, m2 );
            if ( dollo2.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo2.getTotalLosses() != 0 ) {
                return false;
            }
            if ( dollo2.getTotalUnchanged() != 16 ) {
                return false;
            }
            m2.setState( "c", "0", ABSENT );
            dollo2.execute( p2, m2 );
            if ( dollo2.getTotalGains() != 0 ) {
                return false;
            }
            if ( dollo2.getTotalLosses() != 0 ) {
                return false;
            }
            if ( dollo2.getTotalUnchanged() != 17 ) {
                return false;
            }
            m2.setState( "y", "0", PRESENT );
            m2.setState( "e", "0", PRESENT );
            dollo2.execute( p2, m2 );
            if ( dollo2.getTotalGains() != 1 ) {
                return false;
            }
            if ( dollo2.getTotalLosses() != 5 ) {
                return false;
            }
            if ( dollo2.getTotalUnchanged() != 11 ) {
                return false;
            }
            i_m2 = dollo2.getInternalStatesMatrix();
            if ( i_m2.getState( "aby", "0" ) != PRESENT ) {
                return false;
            }
            if ( i_m2.getState( "cad", "0" ) != PRESENT ) {
                return false;
            }
            if ( i_m2.getState( "af", "0" ) != PRESENT ) {
                return false;
            }
            if ( i_m2.getState( "gh", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m2.getState( "ah", "0" ) != ABSENT ) {
                return false;
            }
            if ( i_m2.getState( "r", "0" ) != ABSENT ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testFitchParsimony() {
        try {
            final BinaryStates PRESENT = BinaryStates.PRESENT;
            final BinaryStates ABSENT = BinaryStates.ABSENT;
            final GainLossStates GAIN = GainLossStates.GAIN;
            final GainLossStates LOSS = GainLossStates.LOSS;
            final GainLossStates UNCHANGED_PRESENT = GainLossStates.UNCHANGED_PRESENT;
            final GainLossStates UNCHANGED_ABSENT = GainLossStates.UNCHANGED_ABSENT;
            final FitchParsimony<String> fitch1 = new FitchParsimony<String>();
            final PhylogenyFactory factory1 = ParserBasedPhylogenyFactory.getInstance();
            final String p1_str = "((((((a,b)ab,c)ac,d)ad,(e,f)ef)af,(g,h,i)gi)ai,((j,k,l)jl,(m,n,o)mo,(p,q,r)pr)jr)root";
            final Phylogeny p1 = factory1.create( p1_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<String> m1 = new BasicCharacterStateMatrix<String>( 18, 1 );
            m1.setIdentifier( 0, "a" );
            m1.setIdentifier( 1, "b" );
            m1.setIdentifier( 2, "c" );
            m1.setIdentifier( 3, "d" );
            m1.setIdentifier( 4, "e" );
            m1.setIdentifier( 5, "f" );
            m1.setIdentifier( 6, "g" );
            m1.setIdentifier( 7, "h" );
            m1.setIdentifier( 8, "i" );
            m1.setIdentifier( 9, "j" );
            m1.setIdentifier( 10, "k" );
            m1.setIdentifier( 11, "l" );
            m1.setIdentifier( 12, "m" );
            m1.setIdentifier( 13, "n" );
            m1.setIdentifier( 14, "o" );
            m1.setIdentifier( 15, "p" );
            m1.setIdentifier( 16, "q" );
            m1.setIdentifier( 17, "r" );
            m1.setCharacter( 0, "0" );
            m1.setState( "a", "0", "A" );
            m1.setState( "b", "0", "A" );
            m1.setState( "c", "0", "B" );
            m1.setState( "d", "0", "C" );
            m1.setState( "e", "0", "D" );
            m1.setState( "f", "0", "A" );
            m1.setState( "g", "0", "A" );
            m1.setState( "h", "0", "B" );
            m1.setState( "i", "0", "C" );
            m1.setState( "j", "0", "A" );
            m1.setState( "k", "0", "B" );
            m1.setState( "l", "0", "C" );
            m1.setState( "m", "0", "B" );
            m1.setState( "n", "0", "B" );
            m1.setState( "o", "0", "B" );
            m1.setState( "p", "0", "A" );
            m1.setState( "q", "0", "C" );
            m1.setState( "r", "0", "D" );
            fitch1.setReturnInternalStates( true );
            fitch1.setReturnGainLossMatrix( false );
            fitch1.setRandomize( false );
            fitch1.execute( p1, m1 );
            final CharacterStateMatrix<String> i_m = fitch1.getInternalStatesMatrix();
            final CharacterStateMatrix<List<String>> i_m_all = fitch1.getInternalStatesMatrixPriorToTraceback();
            if ( fitch1.getCost() != 10 ) {
                return false;
            }
            if ( !i_m.getState( "ab", "0" ).equals( "A" ) ) {
                return false;
            }
            if ( !i_m.getState( "ac", "0" ).equals( "A" ) ) {
                return false;
            }
            if ( !i_m.getState( "ad", "0" ).equals( "A" ) ) {
                return false;
            }
            if ( !i_m.getState( "ef", "0" ).equals( "A" ) ) {
                return false;
            }
            if ( !i_m.getState( "ai", "0" ).equals( "A" ) ) {
                return false;
            }
            if ( !i_m.getState( "gi", "0" ).equals( "A" ) ) {
                return false;
            }
            if ( !i_m.getState( "jl", "0" ).equals( "A" ) ) {
                return false;
            }
            if ( !i_m.getState( "mo", "0" ).equals( "B" ) ) {
                return false;
            }
            if ( !i_m.getState( "pr", "0" ).equals( "A" ) ) {
                return false;
            }
            if ( i_m_all.getState( "ab", "0" ).size() != 1 ) {
                return false;
            }
            if ( !i_m_all.getState( "ab", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( i_m_all.getState( "ac", "0" ).size() != 2 ) {
                return false;
            }
            if ( !i_m_all.getState( "ac", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( !i_m_all.getState( "ac", "0" ).contains( "B" ) ) {
                return false;
            }
            if ( i_m_all.getState( "ad", "0" ).size() != 3 ) {
                return false;
            }
            if ( !i_m_all.getState( "ad", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( !i_m_all.getState( "ad", "0" ).contains( "B" ) ) {
                return false;
            }
            if ( !i_m_all.getState( "ad", "0" ).contains( "C" ) ) {
                return false;
            }
            if ( i_m_all.getState( "af", "0" ).size() != 1 ) {
                return false;
            }
            if ( !i_m_all.getState( "af", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( i_m_all.getState( "ef", "0" ).size() != 2 ) {
                return false;
            }
            if ( !i_m_all.getState( "ef", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( !i_m_all.getState( "ef", "0" ).contains( "D" ) ) {
                return false;
            }
            if ( i_m_all.getState( "gi", "0" ).size() != 3 ) {
                return false;
            }
            if ( !i_m_all.getState( "gi", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( !i_m_all.getState( "gi", "0" ).contains( "B" ) ) {
                return false;
            }
            if ( !i_m_all.getState( "gi", "0" ).contains( "C" ) ) {
                return false;
            }
            if ( i_m_all.getState( "ai", "0" ).size() != 1 ) {
                return false;
            }
            if ( !i_m_all.getState( "ai", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( i_m_all.getState( "jl", "0" ).size() != 3 ) {
                return false;
            }
            if ( !i_m_all.getState( "jl", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( !i_m_all.getState( "jl", "0" ).contains( "B" ) ) {
                return false;
            }
            if ( !i_m_all.getState( "jl", "0" ).contains( "C" ) ) {
                return false;
            }
            if ( i_m_all.getState( "mo", "0" ).size() != 1 ) {
                return false;
            }
            if ( !i_m_all.getState( "mo", "0" ).contains( "B" ) ) {
                return false;
            }
            if ( i_m_all.getState( "pr", "0" ).size() != 3 ) {
                return false;
            }
            if ( !i_m_all.getState( "pr", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( !i_m_all.getState( "pr", "0" ).contains( "C" ) ) {
                return false;
            }
            if ( !i_m_all.getState( "pr", "0" ).contains( "D" ) ) {
                return false;
            }
            if ( i_m_all.getState( "jr", "0" ).size() != 4 ) {
                return false;
            }
            if ( !i_m_all.getState( "jr", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( !i_m_all.getState( "jr", "0" ).contains( "B" ) ) {
                return false;
            }
            if ( !i_m_all.getState( "jr", "0" ).contains( "C" ) ) {
                return false;
            }
            if ( !i_m_all.getState( "jr", "0" ).contains( "D" ) ) {
                return false;
            }
            final FitchParsimony<String> fitch2 = new FitchParsimony<String>();
            final PhylogenyFactory factory2 = ParserBasedPhylogenyFactory.getInstance();
            final String p2_str = "((a,b)ab,(c,(d,e)de)cde)r";
            final Phylogeny p2 = factory2.create( p2_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<String> m2 = new BasicCharacterStateMatrix<String>( 5, 1 );
            m2.setIdentifier( 0, "a" );
            m2.setIdentifier( 1, "b" );
            m2.setIdentifier( 2, "c" );
            m2.setIdentifier( 3, "d" );
            m2.setIdentifier( 4, "e" );
            m2.setCharacter( 0, "0" );
            m2.setState( "a", "0", "C" );
            m2.setState( "b", "0", "A" );
            m2.setState( "c", "0", "C" );
            m2.setState( "d", "0", "A" );
            m2.setState( "e", "0", "G" );
            fitch2.setReturnInternalStates( true );
            fitch2.setReturnGainLossMatrix( false );
            fitch2.execute( p2, m2 );
            final CharacterStateMatrix<String> i_m2 = fitch2.getInternalStatesMatrix();
            final CharacterStateMatrix<List<String>> i_m_all2 = fitch2.getInternalStatesMatrixPriorToTraceback();
            if ( fitch2.getCost() != 3 ) {
                return false;
            }
            if ( !i_m2.getState( "ab", "0" ).equals( "A" ) ) {
                return false;
            }
            if ( !i_m2.getState( "de", "0" ).equals( "A" ) ) {
                return false;
            }
            if ( !i_m2.getState( "cde", "0" ).equals( "A" ) ) {
                return false;
            }
            if ( !i_m2.getState( "r", "0" ).equals( "A" ) ) {
                return false;
            }
            if ( i_m_all2.getState( "cde", "0" ).size() != 3 ) {
                return false;
            }
            if ( !i_m_all2.getState( "cde", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( !i_m_all2.getState( "cde", "0" ).contains( "C" ) ) {
                return false;
            }
            if ( !i_m_all2.getState( "cde", "0" ).contains( "G" ) ) {
                return false;
            }
            if ( i_m_all2.getState( "ab", "0" ).size() != 2 ) {
                return false;
            }
            if ( !i_m_all2.getState( "ab", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( !i_m_all2.getState( "ab", "0" ).contains( "C" ) ) {
                return false;
            }
            fitch2.setReturnInternalStates( true );
            fitch2.setReturnGainLossMatrix( false );
            fitch2.setUseLast( true );
            fitch2.execute( p2, m2 );
            final CharacterStateMatrix<String> i_m21 = fitch2.getInternalStatesMatrix();
            final CharacterStateMatrix<List<String>> i_m_all21 = fitch2.getInternalStatesMatrixPriorToTraceback();
            if ( fitch2.getCost() != 3 ) {
                return false;
            }
            if ( !i_m21.getState( "ab", "0" ).equals( "C" ) ) {
                return false;
            }
            if ( !i_m21.getState( "de", "0" ).equals( "G" ) ) {
                return false;
            }
            if ( !i_m21.getState( "cde", "0" ).equals( "C" ) ) {
                return false;
            }
            if ( !i_m21.getState( "r", "0" ).equals( "C" ) ) {
                return false;
            }
            if ( i_m_all21.getState( "cde", "0" ).size() != 3 ) {
                return false;
            }
            if ( !i_m_all21.getState( "cde", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( !i_m_all21.getState( "cde", "0" ).contains( "C" ) ) {
                return false;
            }
            if ( !i_m_all21.getState( "cde", "0" ).contains( "G" ) ) {
                return false;
            }
            final FitchParsimony<String> fitch3 = new FitchParsimony<String>();
            final PhylogenyFactory factory3 = ParserBasedPhylogenyFactory.getInstance();
            final String p3_str = "(((a,b)ab,((c,d)cd,e)cde)abcde,f)r";
            final Phylogeny p3 = factory3.create( p3_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<String> m3 = new BasicCharacterStateMatrix<String>( 6, 1 );
            m3.setIdentifier( 0, "a" );
            m3.setIdentifier( 1, "b" );
            m3.setIdentifier( 2, "c" );
            m3.setIdentifier( 3, "d" );
            m3.setIdentifier( 4, "e" );
            m3.setIdentifier( 5, "f" );
            m3.setCharacter( 0, "0" );
            m3.setState( "a", "0", "C" );
            m3.setState( "b", "0", "U" );
            m3.setState( "c", "0", "G" );
            m3.setState( "d", "0", "U" );
            m3.setState( "e", "0", "A" );
            m3.setState( "f", "0", "A" );
            fitch3.setReturnInternalStates( true );
            fitch3.setReturnGainLossMatrix( false );
            fitch3.execute( p3, m3 );
            final CharacterStateMatrix<String> i_m3 = fitch3.getInternalStatesMatrix();
            final CharacterStateMatrix<List<String>> i_m_all3 = fitch3.getInternalStatesMatrixPriorToTraceback();
            if ( fitch3.getCost() != 4 ) {
                return false;
            }
            if ( !i_m3.getState( "ab", "0" ).equals( "U" ) ) {
                return false;
            }
            if ( !i_m3.getState( "cd", "0" ).equals( "U" ) ) {
                return false;
            }
            if ( !i_m3.getState( "cde", "0" ).equals( "U" ) ) {
                return false;
            }
            if ( !i_m3.getState( "abcde", "0" ).equals( "U" ) ) {
                return false;
            }
            if ( !i_m3.getState( "r", "0" ).equals( "A" ) ) {
                return false;
            }
            if ( i_m_all3.getState( "cde", "0" ).size() != 3 ) {
                return false;
            }
            if ( !i_m_all3.getState( "cde", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( !i_m_all3.getState( "cde", "0" ).contains( "G" ) ) {
                return false;
            }
            if ( !i_m_all3.getState( "cde", "0" ).contains( "U" ) ) {
                return false;
            }
            if ( i_m_all3.getState( "ab", "0" ).size() != 2 ) {
                return false;
            }
            if ( !i_m_all3.getState( "ab", "0" ).contains( "C" ) ) {
                return false;
            }
            if ( !i_m_all3.getState( "ab", "0" ).contains( "U" ) ) {
                return false;
            }
            if ( i_m_all3.getState( "cd", "0" ).size() != 2 ) {
                return false;
            }
            if ( !i_m_all3.getState( "cd", "0" ).contains( "G" ) ) {
                return false;
            }
            if ( !i_m_all3.getState( "cd", "0" ).contains( "U" ) ) {
                return false;
            }
            if ( i_m_all3.getState( "abcde", "0" ).size() != 1 ) {
                return false;
            }
            if ( !i_m_all3.getState( "abcde", "0" ).contains( "U" ) ) {
                return false;
            }
            if ( i_m_all3.getState( "r", "0" ).size() != 2 ) {
                return false;
            }
            if ( !i_m_all3.getState( "r", "0" ).contains( "A" ) ) {
                return false;
            }
            if ( !i_m_all3.getState( "r", "0" ).contains( "U" ) ) {
                return false;
            }
            final FitchParsimony<BinaryStates> fitch4 = new FitchParsimony<BinaryStates>();
            final PhylogenyFactory factory4 = ParserBasedPhylogenyFactory.getInstance();
            final String p4_str = "(((a,b)ab,((c,d)cd,e)cde)abcde,f)r";
            final Phylogeny p4 = factory4.create( p4_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<BinaryStates> m4 = new BasicCharacterStateMatrix<BinaryStates>( 6, 1 );
            m4.setIdentifier( 0, "a" );
            m4.setIdentifier( 1, "b" );
            m4.setIdentifier( 2, "c" );
            m4.setIdentifier( 3, "d" );
            m4.setIdentifier( 4, "e" );
            m4.setIdentifier( 5, "f" );
            m4.setCharacter( 0, "0" );
            m4.setState( "a", "0", PRESENT );
            m4.setState( "b", "0", ABSENT );
            m4.setState( "c", "0", PRESENT );
            m4.setState( "d", "0", PRESENT );
            m4.setState( "e", "0", ABSENT );
            m4.setState( "f", "0", ABSENT );
            fitch4.setReturnInternalStates( true );
            fitch4.setReturnGainLossMatrix( true );
            fitch4.execute( p4, m4 );
            final CharacterStateMatrix<GainLossStates> gl_m_4 = fitch4.getGainLossMatrix();
            if ( fitch4.getCost() != 2 ) {
                return false;
            }
            if ( fitch4.getTotalLosses() != 0 ) {
                return false;
            }
            if ( fitch4.getTotalGains() != 2 ) {
                return false;
            }
            if ( fitch4.getTotalUnchanged() != 9 ) {
                return false;
            }
            if ( gl_m_4.getState( "a", "0" ) != GAIN ) {
                return false;
            }
            if ( gl_m_4.getState( "b", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_4.getState( "ab", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_4.getState( "cd", "0" ) != GAIN ) {
                return false;
            }
            if ( gl_m_4.getState( "r", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            final FitchParsimony<BinaryStates> fitch5 = new FitchParsimony<BinaryStates>();
            final PhylogenyFactory factory5 = ParserBasedPhylogenyFactory.getInstance();
            final String p5_str = "(((a,b)ab,((c,d)cd,e)cde)abcde,f)r";
            final Phylogeny p5 = factory5.create( p5_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<BinaryStates> m5 = new BasicCharacterStateMatrix<BinaryStates>( 6, 1 );
            m5.setIdentifier( 0, "a" );
            m5.setIdentifier( 1, "b" );
            m5.setIdentifier( 2, "c" );
            m5.setIdentifier( 3, "d" );
            m5.setIdentifier( 4, "e" );
            m5.setIdentifier( 5, "f" );
            m5.setCharacter( 0, "0" );
            m5.setState( "a", "0", PRESENT );
            m5.setState( "b", "0", ABSENT );
            m5.setState( "c", "0", PRESENT );
            m5.setState( "d", "0", ABSENT );
            m5.setState( "e", "0", PRESENT );
            m5.setState( "f", "0", ABSENT );
            fitch5.setReturnInternalStates( true );
            fitch5.setReturnGainLossMatrix( true );
            fitch5.execute( p5, m5 );
            final CharacterStateMatrix<GainLossStates> gl_m_5 = fitch5.getGainLossMatrix();
            if ( fitch5.getCost() != 3 ) {
                return false;
            }
            if ( fitch5.getTotalLosses() != 2 ) {
                return false;
            }
            if ( fitch5.getTotalGains() != 1 ) {
                return false;
            }
            if ( fitch5.getTotalUnchanged() != 8 ) {
                return false;
            }
            if ( gl_m_5.getState( "abcde", "0" ) != GAIN ) {
                return false;
            }
            if ( gl_m_5.getState( "a", "0" ) != UNCHANGED_PRESENT ) {
                return false;
            }
            if ( gl_m_5.getState( "b", "0" ) != LOSS ) {
                return false;
            }
            if ( gl_m_5.getState( "d", "0" ) != LOSS ) {
                return false;
            }
            if ( gl_m_5.getState( "r", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            final FitchParsimony<BinaryStates> fitch6 = new FitchParsimony<BinaryStates>();
            final PhylogenyFactory factory6 = ParserBasedPhylogenyFactory.getInstance();
            final String p6_str = "(((a,b)ab,((c,d)cd,e)cde)abcde,f)r";
            final Phylogeny p6 = factory6.create( p6_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<BinaryStates> m6 = new BasicCharacterStateMatrix<BinaryStates>( 6, 1 );
            m6.setIdentifier( 0, "a" );
            m6.setIdentifier( 1, "b" );
            m6.setIdentifier( 2, "c" );
            m6.setIdentifier( 3, "d" );
            m6.setIdentifier( 4, "e" );
            m6.setIdentifier( 5, "f" );
            m6.setCharacter( 0, "0" );
            m6.setState( "a", "0", PRESENT );
            m6.setState( "b", "0", ABSENT );
            m6.setState( "c", "0", PRESENT );
            m6.setState( "d", "0", PRESENT );
            m6.setState( "e", "0", ABSENT );
            m6.setState( "f", "0", PRESENT );
            fitch6.setReturnInternalStates( true );
            fitch6.setReturnGainLossMatrix( true );
            fitch6.execute( p6, m6 );
            final CharacterStateMatrix<GainLossStates> gl_m_6 = fitch6.getGainLossMatrix();
            if ( fitch6.getCost() != 2 ) {
                return false;
            }
            if ( fitch6.getTotalLosses() != 2 ) {
                return false;
            }
            if ( fitch6.getTotalGains() != 0 ) {
                return false;
            }
            if ( fitch6.getTotalUnchanged() != 9 ) {
                return false;
            }
            if ( gl_m_6.getState( "abcde", "0" ) != UNCHANGED_PRESENT ) {
                return false;
            }
            if ( gl_m_6.getState( "r", "0" ) != UNCHANGED_PRESENT ) {
                return false;
            }
            if ( gl_m_6.getState( "b", "0" ) != LOSS ) {
                return false;
            }
            if ( gl_m_6.getState( "e", "0" ) != LOSS ) {
                return false;
            }
            final FitchParsimony<BinaryStates> fitch7 = new FitchParsimony<BinaryStates>();
            final PhylogenyFactory factory7 = ParserBasedPhylogenyFactory.getInstance();
            final String p7_str = "(((a,b)ab,(c,d)cd)abcd,((e,f)ef,(g,h)gh)efgh)r";
            final Phylogeny p7 = factory7.create( p7_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<BinaryStates> m7 = new BasicCharacterStateMatrix<BinaryStates>( 8, 1 );
            m7.setIdentifier( 0, "a" );
            m7.setIdentifier( 1, "b" );
            m7.setIdentifier( 2, "c" );
            m7.setIdentifier( 3, "d" );
            m7.setIdentifier( 4, "e" );
            m7.setIdentifier( 5, "f" );
            m7.setIdentifier( 6, "g" );
            m7.setIdentifier( 7, "h" );
            m7.setCharacter( 0, "0" );
            m7.setState( "a", "0", PRESENT );
            m7.setState( "b", "0", ABSENT );
            m7.setState( "c", "0", PRESENT );
            m7.setState( "d", "0", ABSENT );
            m7.setState( "e", "0", PRESENT );
            m7.setState( "f", "0", ABSENT );
            m7.setState( "g", "0", PRESENT );
            m7.setState( "h", "0", ABSENT );
            fitch7.setReturnInternalStates( true );
            fitch7.setReturnGainLossMatrix( true );
            fitch7.execute( p7, m7 );
            final CharacterStateMatrix<GainLossStates> gl_m_7 = fitch7.getGainLossMatrix();
            if ( fitch7.getCost() != 4 ) {
                return false;
            }
            if ( fitch7.getTotalLosses() != 0 ) {
                return false;
            }
            if ( fitch7.getTotalGains() != 4 ) {
                return false;
            }
            if ( fitch7.getTotalUnchanged() != 11 ) {
                return false;
            }
            if ( gl_m_7.getState( "a", "0" ) != GAIN ) {
                return false;
            }
            if ( gl_m_7.getState( "c", "0" ) != GAIN ) {
                return false;
            }
            if ( gl_m_7.getState( "e", "0" ) != GAIN ) {
                return false;
            }
            if ( gl_m_7.getState( "g", "0" ) != GAIN ) {
                return false;
            }
            if ( gl_m_7.getState( "r", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            fitch7.setReturnInternalStates( true );
            fitch7.setReturnGainLossMatrix( true );
            fitch7.setUseLast( true );
            fitch7.execute( p7, m7 );
            final CharacterStateMatrix<GainLossStates> gl_m_71 = fitch7.getGainLossMatrix();
            if ( fitch7.getCost() != 4 ) {
                return false;
            }
            if ( fitch7.getTotalLosses() != 4 ) {
                return false;
            }
            if ( fitch7.getTotalGains() != 0 ) {
                return false;
            }
            if ( fitch7.getTotalUnchanged() != 11 ) {
                return false;
            }
            if ( gl_m_71.getState( "b", "0" ) != LOSS ) {
                return false;
            }
            if ( gl_m_71.getState( "d", "0" ) != LOSS ) {
                return false;
            }
            if ( gl_m_71.getState( "f", "0" ) != LOSS ) {
                return false;
            }
            if ( gl_m_71.getState( "h", "0" ) != LOSS ) {
                return false;
            }
            if ( gl_m_71.getState( "r", "0" ) != UNCHANGED_PRESENT ) {
                return false;
            }
            final FitchParsimony<BinaryStates> fitch8 = new FitchParsimony<BinaryStates>();
            final PhylogenyFactory factory8 = ParserBasedPhylogenyFactory.getInstance();
            final String p8_str = "(((a,b)ab,(c,d)cd)abcd,((e,f)ef,(g,h)gh)efgh)r";
            final Phylogeny p8 = factory8.create( p8_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<BinaryStates> m8 = new BasicCharacterStateMatrix<BinaryStates>( 8, 1 );
            m8.setIdentifier( 0, "a" );
            m8.setIdentifier( 1, "b" );
            m8.setIdentifier( 2, "c" );
            m8.setIdentifier( 3, "d" );
            m8.setIdentifier( 4, "e" );
            m8.setIdentifier( 5, "f" );
            m8.setIdentifier( 6, "g" );
            m8.setIdentifier( 7, "h" );
            m8.setCharacter( 0, "0" );
            m8.setState( "a", "0", PRESENT );
            m8.setState( "b", "0", PRESENT );
            m8.setState( "c", "0", PRESENT );
            m8.setState( "d", "0", ABSENT );
            m8.setState( "e", "0", ABSENT );
            m8.setState( "f", "0", ABSENT );
            m8.setState( "g", "0", ABSENT );
            m8.setState( "h", "0", ABSENT );
            fitch8.setReturnInternalStates( true );
            fitch8.setReturnGainLossMatrix( true );
            fitch8.execute( p8, m8 );
            final CharacterStateMatrix<GainLossStates> gl_m_8 = fitch8.getGainLossMatrix();
            if ( fitch8.getCost() != 2 ) {
                return false;
            }
            if ( fitch8.getTotalLosses() != 1 ) {
                return false;
            }
            if ( fitch8.getTotalGains() != 1 ) {
                return false;
            }
            if ( fitch8.getTotalUnchanged() != 13 ) {
                return false;
            }
            if ( gl_m_8.getState( "d", "0" ) != LOSS ) {
                return false;
            }
            if ( gl_m_8.getState( "abcd", "0" ) != GAIN ) {
                return false;
            }
            final FitchParsimony<BinaryStates> fitch9 = new FitchParsimony<BinaryStates>();
            final PhylogenyFactory factory9 = ParserBasedPhylogenyFactory.getInstance();
            final String p9_str = "(((a,b)ab,c)abc,d)abcd";
            final Phylogeny p9 = factory9.create( p9_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<BinaryStates> m9 = new BasicCharacterStateMatrix<BinaryStates>( 4, 1 );
            m9.setIdentifier( 0, "a" );
            m9.setIdentifier( 1, "b" );
            m9.setIdentifier( 2, "c" );
            m9.setIdentifier( 3, "d" );
            m9.setCharacter( 0, "0" );
            m9.setState( "a", "0", PRESENT );
            m9.setState( "b", "0", ABSENT );
            m9.setState( "c", "0", PRESENT );
            m9.setState( "d", "0", ABSENT );
            fitch9.setReturnInternalStates( true );
            fitch9.setReturnGainLossMatrix( true );
            fitch9.setUseLast( false );
            fitch9.execute( p9, m9 );
            final CharacterStateMatrix<GainLossStates> gl_m_9a = fitch9.getGainLossMatrix();
            if ( fitch9.getCost() != 2 ) {
                return false;
            }
            if ( fitch9.getTotalLosses() != 1 ) {
                return false;
            }
            if ( fitch9.getTotalGains() != 1 ) {
                return false;
            }
            if ( fitch9.getTotalUnchanged() != 5 ) {
                return false;
            }
            if ( gl_m_9a.getState( "a", "0" ) != UNCHANGED_PRESENT ) {
                return false;
            }
            if ( gl_m_9a.getState( "b", "0" ) != LOSS ) {
                return false;
            }
            if ( gl_m_9a.getState( "c", "0" ) != UNCHANGED_PRESENT ) {
                return false;
            }
            if ( gl_m_9a.getState( "d", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_9a.getState( "ab", "0" ) != UNCHANGED_PRESENT ) {
                return false;
            }
            if ( gl_m_9a.getState( "abc", "0" ) != GAIN ) {
                return false;
            }
            if ( gl_m_9a.getState( "abcd", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            fitch9.setUseLast( true );
            fitch9.execute( p9, m9 );
            final CharacterStateMatrix<GainLossStates> gl_m_9b = fitch9.getGainLossMatrix();
            if ( fitch9.getCost() != 2 ) {
                return false;
            }
            if ( fitch9.getTotalLosses() != 2 ) {
                return false;
            }
            if ( fitch9.getTotalGains() != 0 ) {
                return false;
            }
            if ( fitch9.getTotalUnchanged() != 5 ) {
                return false;
            }
            if ( gl_m_9b.getState( "a", "0" ) != UNCHANGED_PRESENT ) {
                return false;
            }
            if ( gl_m_9b.getState( "b", "0" ) != LOSS ) {
                return false;
            }
            if ( gl_m_9b.getState( "c", "0" ) != UNCHANGED_PRESENT ) {
                return false;
            }
            if ( gl_m_9b.getState( "d", "0" ) != LOSS ) {
                return false;
            }
            if ( gl_m_9b.getState( "ab", "0" ) != UNCHANGED_PRESENT ) {
                return false;
            }
            if ( gl_m_9b.getState( "abc", "0" ) != UNCHANGED_PRESENT ) {
                return false;
            }
            if ( gl_m_9b.getState( "abcd", "0" ) != UNCHANGED_PRESENT ) {
                return false;
            }
            fitch9.setUseLast( false );
            fitch9.setRandomize( true );
            fitch9.setRandomNumberSeed( 8722445 );
            fitch9.execute( p9, m9 );
            fitch9.getGainLossMatrix();
            if ( fitch9.getCost() != 2 ) {
                return false;
            }
            if ( fitch9.getTotalLosses() != 1 ) {
                return false;
            }
            if ( fitch9.getTotalGains() != 1 ) {
                return false;
            }
            if ( fitch9.getTotalUnchanged() != 5 ) {
                return false;
            }
            final FitchParsimony<BinaryStates> fitch10 = new FitchParsimony<BinaryStates>();
            final PhylogenyFactory factory10 = ParserBasedPhylogenyFactory.getInstance();
            final String p10_str = "((((a,b)ab,c)abc,d)abcd,e)abcde";
            final Phylogeny p10 = factory10.create( p10_str, new NHXParser() )[ 0 ];
            final CharacterStateMatrix<BinaryStates> m10 = new BasicCharacterStateMatrix<BinaryStates>( 5, 1 );
            m10.setIdentifier( 0, "a" );
            m10.setIdentifier( 1, "b" );
            m10.setIdentifier( 2, "c" );
            m10.setIdentifier( 3, "d" );
            m10.setIdentifier( 4, "e" );
            m10.setCharacter( 0, "0" );
            m10.setState( "a", "0", PRESENT );
            m10.setState( "b", "0", ABSENT );
            m10.setState( "c", "0", ABSENT );
            m10.setState( "d", "0", PRESENT );
            m10.setState( "e", "0", ABSENT );
            fitch10.setReturnInternalStates( true );
            fitch10.setReturnGainLossMatrix( true );
            fitch10.setUseLast( false );
            fitch10.execute( p10, m10 );
            final CharacterStateMatrix<GainLossStates> gl_m_10a = fitch10.getGainLossMatrix();
            if ( fitch10.getCost() != 2 ) {
                return false;
            }
            if ( fitch10.getTotalLosses() != 0 ) {
                return false;
            }
            if ( fitch10.getTotalGains() != 2 ) {
                return false;
            }
            if ( fitch10.getTotalUnchanged() != 7 ) {
                return false;
            }
            if ( gl_m_10a.getState( "a", "0" ) != GAIN ) {
                return false;
            }
            if ( gl_m_10a.getState( "b", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_10a.getState( "c", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_10a.getState( "d", "0" ) != GAIN ) {
                return false;
            }
            if ( gl_m_10a.getState( "e", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_10a.getState( "ab", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_10a.getState( "abc", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_10a.getState( "abcd", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_10a.getState( "abcde", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            fitch10.setUseLast( true );
            fitch10.execute( p10, m10 );
            final CharacterStateMatrix<GainLossStates> gl_m_10b = fitch10.getGainLossMatrix();
            if ( fitch10.getCost() != 2 ) {
                return false;
            }
            if ( fitch10.getTotalLosses() != 0 ) {
                return false;
            }
            if ( fitch10.getTotalGains() != 2 ) {
                return false;
            }
            if ( fitch10.getTotalUnchanged() != 7 ) {
                return false;
            }
            if ( gl_m_10b.getState( "a", "0" ) != GAIN ) {
                return false;
            }
            if ( gl_m_10b.getState( "b", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_10b.getState( "c", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_10b.getState( "d", "0" ) != GAIN ) {
                return false;
            }
            if ( gl_m_10b.getState( "e", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_10b.getState( "ab", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_10b.getState( "abc", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_10b.getState( "abcd", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
            if ( gl_m_10b.getState( "abcde", "0" ) != UNCHANGED_ABSENT ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNeighborJoining() {
        try {
            BasicSymmetricalDistanceMatrix m = new BasicSymmetricalDistanceMatrix( 6 );
            m.setRow( "5", 1 );
            m.setRow( "4 7", 2 );
            m.setRow( "7 10 7", 3 );
            m.setRow( "6 9 6 5", 4 );
            m.setRow( "8 11 8 9 8", 5 );
            m.setIdentifier( 0, "A" );
            m.setIdentifier( 1, "B" );
            m.setIdentifier( 2, "C" );
            m.setIdentifier( 3, "D" );
            m.setIdentifier( 4, "E" );
            m.setIdentifier( 5, "F" );
            final NeighborJoining nj = NeighborJoining.createInstance();
            nj.setVerbose( false );
            nj.execute( m );
            m = new BasicSymmetricalDistanceMatrix( 7 );
            m.setIdentifier( 0, "Bovine" );
            m.setIdentifier( 1, "Mouse" );
            m.setIdentifier( 2, "Gibbon" );
            m.setIdentifier( 3, "Orang" );
            m.setIdentifier( 4, "Gorilla" );
            m.setIdentifier( 5, "Chimp" );
            m.setIdentifier( 6, "Human" );
            m.setRow( "0.00000 1.68660 1.71980 1.66060 1.52430 1.60430 1.59050", 0 );
            m.setRow( "1.68660 0.00000 1.52320 1.48410 1.44650 1.43890 1.46290", 1 );
            m.setRow( "1.71980 1.52320 0.00000 0.71150 0.59580 0.61790 0.55830", 2 );
            m.setRow( "1.66060 1.48410 0.71150 0.00000 0.46310 0.50610 0.47100", 3 );
            m.setRow( "1.52430 1.44650 0.59580 0.46310 0.00000 0.34840 0.30830", 4 );
            m.setRow( "1.60430 1.43890 0.61790 0.50610 0.34840 0.00000 0.26920", 5 );
            m.setRow( "1.59050 1.46290 0.55830 0.47100 0.30830 0.26920 0.00000", 6 );
            nj.execute( m );
            m = new BasicSymmetricalDistanceMatrix( 4 );
            m.setIdentifier( 0, "A" );
            m.setIdentifier( 1, "B" );
            m.setIdentifier( 2, "C" );
            m.setIdentifier( 3, "D" );
            m.setRow( "0.00 0.95 0.17 0.98", 0 );
            m.setRow( "0.95 0.00 1.02 1.83", 1 );
            m.setRow( "0.17 1.02 0.00 1.01", 2 );
            m.setRow( "0.98 1.83 1.01 0.00", 3 );
            final Phylogeny p3 = nj.execute( m );
            //
            // -- A 0.05
            // - |0.01
            // ----------------------- B 0.90
            //
            // --- C 0.10
            // - |0.01
            // ------------------------- D 0.91
            p3.reRoot( p3.getNode( "C" ).getParent() );
            if ( !isEqual( p3.getNode( "A" ).getDistanceToParent(), 0.05 ) ) {
                return false;
            }
            if ( !isEqual( p3.getNode( "B" ).getDistanceToParent(), 0.90 ) ) {
                return false;
            }
            if ( !isEqual( p3.getNode( "C" ).getDistanceToParent(), 0.10 ) ) {
                return false;
            }
            if ( !isEqual( p3.getNode( "D" ).getDistanceToParent(), 0.91 ) ) {
                return false;
            }
            if ( TIME ) {
                timeNeighborJoining();
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSymmetricalDistanceMatrixParser() {
        try {
            final String l = ForesterUtil.getLineSeparator();
            StringBuffer source = new StringBuffer();
            source.append( " 4" + l );
            source.append( "A 0 0 0 0" + l );
            source.append( "B 1 0 0 0" + l );
            source.append( "C 2 4 0 0" + l );
            source.append( "D 3 5 6 0" + l );
            source.append( l );
            source.append( " 4" + l );
            source.append( "A 0   11  12  13" + l );
            source.append( "B 11  0   14  15" + l );
            source.append( "C 12  14  0   16" + l );
            source.append( "D 13  15  16  0" + l );
            source.append( l );
            source.append( l );
            source.append( "     " + l );
            source.append( " 4" + l );
            source.append( " A        0     " + l );
            source.append( " B            21 0" + l );
            source.append( " C            22 24    0  " + l );
            source.append( " # 2 222 2 2 " + l );
            source.append( " D            23 25 26 0" + l );
            source.append( l );
            source.append( l );
            source.append( "     " + l );
            final SymmetricalDistanceMatrixParser p0 = SymmetricalDistanceMatrixParser.createInstance();
            final DistanceMatrix[] ma0 = p0.parse( source.toString() );
            if ( ma0.length != 3 ) {
                return false;
            }
            if ( !isEqual( ma0[ 0 ].getValue( 0, 0 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 0 ].getValue( 1, 0 ), 1 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 0 ].getValue( 2, 0 ), 2 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 0 ].getValue( 3, 0 ), 3 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 0 ].getValue( 0, 1 ), 1 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 0 ].getValue( 1, 1 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 0 ].getValue( 2, 1 ), 4 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 0 ].getValue( 3, 1 ), 5 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 1 ].getValue( 0, 0 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 1 ].getValue( 1, 0 ), 11 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 1 ].getValue( 2, 0 ), 12 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 1 ].getValue( 3, 0 ), 13 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 1 ].getValue( 0, 1 ), 11 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 1 ].getValue( 1, 1 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 1 ].getValue( 2, 1 ), 14 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 1 ].getValue( 3, 1 ), 15 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 2 ].getValue( 0, 0 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 2 ].getValue( 1, 0 ), 21 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 2 ].getValue( 2, 0 ), 22 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 2 ].getValue( 3, 0 ), 23 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 2 ].getValue( 0, 1 ), 21 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 2 ].getValue( 1, 1 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 2 ].getValue( 2, 1 ), 24 ) ) {
                return false;
            }
            if ( !isEqual( ma0[ 2 ].getValue( 3, 1 ), 25 ) ) {
                return false;
            }
            source = new StringBuffer();
            source.append( "A 0 0 0 0" + l );
            source.append( "B 1 0 0 0" + l );
            source.append( "C 2 4 0 0" + l );
            source.append( "D 3 5 6 0" + l );
            source.append( " " + l );
            source.append( "A 0   11  12  13" + l );
            source.append( "B 11  0   14  15" + l );
            source.append( "C 12  14  0   16" + l );
            source.append( "D 13  15  16  0" + l );
            source.append( l );
            source.append( " A        0     " + l );
            source.append( " B            21 0" + l );
            source.append( " C            22 24    0  " + l );
            source.append( " # 2 222 2 2 " + l );
            source.append( " D            23 25 26 0" + l );
            final DistanceMatrix[] ma1 = p0.parse( source.toString() );
            if ( ma1.length != 3 ) {
                return false;
            }
            if ( !isEqual( ma1[ 0 ].getValue( 0, 0 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 0 ].getValue( 1, 0 ), 1 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 0 ].getValue( 2, 0 ), 2 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 0 ].getValue( 3, 0 ), 3 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 0 ].getValue( 0, 1 ), 1 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 0 ].getValue( 1, 1 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 0 ].getValue( 2, 1 ), 4 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 0 ].getValue( 3, 1 ), 5 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 1 ].getValue( 0, 0 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 1 ].getValue( 1, 0 ), 11 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 1 ].getValue( 2, 0 ), 12 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 1 ].getValue( 3, 0 ), 13 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 1 ].getValue( 0, 1 ), 11 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 1 ].getValue( 1, 1 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 1 ].getValue( 2, 1 ), 14 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 1 ].getValue( 3, 1 ), 15 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 2 ].getValue( 0, 0 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 2 ].getValue( 1, 0 ), 21 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 2 ].getValue( 2, 0 ), 22 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 2 ].getValue( 3, 0 ), 23 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 2 ].getValue( 0, 1 ), 21 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 2 ].getValue( 1, 1 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 2 ].getValue( 2, 1 ), 24 ) ) {
                return false;
            }
            if ( !isEqual( ma1[ 2 ].getValue( 3, 1 ), 25 ) ) {
                return false;
            }
            source = new StringBuffer();
            source.append( "A 0" + l );
            source.append( "B 10 0" + l );
            final DistanceMatrix[] ma2 = p0.parse( source.toString() );
            if ( ma2.length != 1 ) {
                return false;
            }
            if ( !isEqual( ma2[ 0 ].getValue( 0, 1 ), 10 ) ) {
                return false;
            }
            source = new StringBuffer();
            source.append( " " + l );
            source.append( "#" + l );
            final DistanceMatrix[] ma3 = p0.parse( source.toString() );
            if ( ma3.length != 0 ) {
                return false;
            }
            source = new StringBuffer();
            source.append( " " + l );
            source.append( "A 0   11  12  13" + l );
            source.append( "B     0   14  15" + l );
            source.append( "C         0   16" + l );
            source.append( "D              0" + l );
            source.append( l );
            source.append( "A 0 21  22  23" + l );
            source.append( "B 0 24  25" + l );
            source.append( "C 0 26" + l );
            source.append( "D 0" + l );
            p0.setInputMatrixType( SymmetricalDistanceMatrixParser.InputMatrixType.UPPER_TRIANGLE );
            final DistanceMatrix[] ma4 = p0.parse( source );
            if ( ma4.length != 2 ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 0, 0 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 1, 0 ), 11 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 2, 0 ), 12 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 3, 0 ), 13 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 0, 1 ), 11 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 1, 1 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 2, 1 ), 14 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 3, 1 ), 15 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 0, 2 ), 12 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 1, 2 ), 14 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 2, 2 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 3, 2 ), 16 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 0, 3 ), 13 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 1, 3 ), 15 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 2, 3 ), 16 ) ) {
                return false;
            }
            if ( !isEqual( ma4[ 0 ].getValue( 3, 3 ), 0 ) ) {
                return false;
            }
            source = new StringBuffer();
            source.append( " 4 " + l );
            source.append( "A 0   11  12  13" + l );
            source.append( "B     0   14  15" + l );
            source.append( "C         0   16" + l );
            source.append( "D              0" + l );
            source.append( " 4" + l );
            source.append( "A 0 21  22  23" + l );
            source.append( "B 0 24  25" + l );
            source.append( "C 0 26" + l );
            source.append( "D 0" + l );
            source.append( "     " + l );
            source.append( " 4" + l );
            source.append( "A 0 21  22  23" + l );
            source.append( "B 0 24  25" + l );
            source.append( "C 0 26" + l );
            source.append( "D 0" + l );
            source.append( l );
            source.append( "A 0 21  22  23" + l );
            source.append( "B 0 24  25" + l );
            source.append( "C 0 26" + l );
            source.append( "D 0" + l );
            p0.setInputMatrixType( SymmetricalDistanceMatrixParser.InputMatrixType.UPPER_TRIANGLE );
            final DistanceMatrix[] ma5 = p0.parse( source );
            if ( ma5.length != 4 ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 0, 0 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 1, 0 ), 11 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 2, 0 ), 12 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 3, 0 ), 13 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 0, 1 ), 11 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 1, 1 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 2, 1 ), 14 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 3, 1 ), 15 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 0, 2 ), 12 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 1, 2 ), 14 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 2, 2 ), 0 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 3, 2 ), 16 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 0, 3 ), 13 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 1, 3 ), 15 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 2, 3 ), 16 ) ) {
                return false;
            }
            if ( !isEqual( ma5[ 0 ].getValue( 3, 3 ), 0 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static void timeNeighborJoining() {
        final NeighborJoining nj = NeighborJoining.createInstance();
        for( int n = 3; n <= 12; ++n ) {
            final int x = ( int ) Math.pow( 2, n );
            final BasicSymmetricalDistanceMatrix mt = new BasicSymmetricalDistanceMatrix( x );
            mt.randomize( new Date().getTime() );
            final long start_time = new Date().getTime();
            nj.execute( mt );
            System.out.println( "Size: " + x + " -> " + ( new Date().getTime() - start_time ) + "ms." );
        }
    }
}
