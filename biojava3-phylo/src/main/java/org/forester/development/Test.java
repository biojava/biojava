// $Id: Test.java,v 1.24 2009/12/11 01:13:41 cmzmasek Exp $
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

package org.forester.development;

import java.util.Date;

import org.forester.phylogeny.Phylogeny;

/*
 * *
 */
public class Test {

    public static void main( final String[] args ) {
        int failed = 0;
        int succeeded = 0;
        final long start_time = new Date().getTime();
        System.out.print( "Amino acid sequence: " );
        if ( Test.testAminoAcidSequence() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.println( "Creation of balanced phylogeny: " );
        if ( Test.testCreateBalancedPhylogeny() ) {
            System.out.println( "OK." );
            succeeded++;
        }
        else {
            System.out.println( "failed." );
            failed++;
        }
        System.out.println( "\nTime requirement:  " + ( new Date().getTime() - start_time ) + "ms." );
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
    }

    private static boolean testAminoAcidSequence() {
        try {
            final AminoAcidSequence aa0 = new AminoAcidSequence( 12 );
            if ( aa0.getLength() != 12 ) {
                return false;
            }
            final AminoAcidSequence aa1 = new AminoAcidSequence( "aa1", " aAklm-?xX*z$#" );
            final AminoAcidSequence aa2 = aa1.copy();
            aa1.setStateAt( 0, ( byte ) 1 );
            aa1.setResidueAt( 1, 'r' );
            if ( aa1.getLength() != 14 ) {
                return false;
            }
            if ( aa1.getResidueAt( 0 ) != 'R' ) {
                return false;
            }
            if ( aa1.getResidueAt( 1 ) != 'R' ) {
                return false;
            }
            if ( aa2.getResidueAt( 0 ) != '?' ) {
                return false;
            }
            if ( aa2.getResidueAt( 1 ) != 'A' ) {
                return false;
            }
            if ( aa2.getResidueAt( 6 ) != '-' ) {
                return false;
            }
            if ( aa2.getResidueAt( 7 ) != '?' ) {
                return false;
            }
            if ( aa2.getResidueAt( 8 ) != '?' ) {
                return false;
            }
            if ( aa2.getResidueAt( 9 ) != '?' ) {
                return false;
            }
            if ( aa2.getResidueAt( 10 ) != '*' ) {
                return false;
            }
            if ( !aa2.getSequenceAsString().equals( "?AAKLM-???*???" ) ) {
                return false;
            }
            if ( !aa1.getSequenceAsString().equals( "RRAKLM-???*???" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static boolean testCreateBalancedPhylogeny() {
        try {
            final Phylogeny p0 = DevelopmentTools.createBalancedPhylogeny( 6, 5 );
            if ( p0.getRoot().getNumberOfDescendants() != 5 ) {
                return false;
            }
            if ( p0.getNumberOfExternalNodes() != 15625 ) {
                return false;
            }
            final Phylogeny p1 = DevelopmentTools.createBalancedPhylogeny( 2, 10 );
            if ( p1.getRoot().getNumberOfDescendants() != 10 ) {
                return false;
            }
            if ( p1.getNumberOfExternalNodes() != 100 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return false;
        }
        return true;
    }
}
