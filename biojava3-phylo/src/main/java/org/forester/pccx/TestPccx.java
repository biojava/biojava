// $Id: TestPccx.java,v 1.9 2009/10/26 23:29:39 cmzmasek Exp $
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

package org.forester.pccx;

import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;

/*
 * @author Christian M. Zmasek
 */
public class TestPccx {

    private final static double ZERO_DIFF = 1.0E-6;

    private static boolean isEqual( final double a, final double b ) {
        return ( ( Math.abs( a - b ) ) < TestPccx.ZERO_DIFF );
    }

    public static boolean test() {
        if ( !TestPccx.testExternalNodeBasedCoverage() ) {
            return false;
        }
        return true;
    }

    private static boolean testExternalNodeBasedCoverage() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final String ps1 = "((((A:0.1,B:0.7):0.2,C:1.0):2.0,D:1.7):1.3,((E:0.3,F:0.4):1.1,(G:0.5,H:0.6):1.2):1.4,X:2.0)";
            final Phylogeny p1 = factory.create( ps1, new NHXParser() )[ 0 ];
            final List<Phylogeny> phylogenies = new ArrayList<Phylogeny>();
            final List<String> names = new ArrayList<String>();
            phylogenies.add( p1 );
            names.add( "A" );
            names.add( "A" );
            final CoverageCalculationOptions options = new ExternalNodeBasedCoverageMethodOptions( "org.forester.pccx.BranchCountingBasedScoringMethod" );
            final CoverageCalculator cc = CoverageCalculator.getInstance( new ExternalNodeBasedCoverageMethod(),
                                                                          options );
            Coverage cov = cc.calculateCoverage( phylogenies, names, false );
            if ( !TestPccx.isEqual( cov.getScore(), ( 1.0 + 1.0 / 2 + 1.0 / 3 + 1.0 / 4 + 1.0 / 7 + 1.0 / 7 + 1.0 / 7
                    + 1.0 / 7 + 1.0 / 5 ) / 9 ) ) {
                return false;
            }
            names.add( "B" );
            names.add( "B" );
            cov = cc.calculateCoverage( phylogenies, names, false );
            if ( !TestPccx.isEqual( cov.getScore(), ( 1.0 + 1.0 + 1.0 / 3 + 1.0 / 4 + 1.0 / 7 + 1.0 / 7 + 1.0 / 7 + 1.0
                    / 7 + 1.0 / 5 ) / 9 ) ) {
                return false;
            }
            names.add( "G" );
            cov = cc.calculateCoverage( phylogenies, names, false );
            if ( !TestPccx
                    .isEqual( cov.getScore(),
                              ( 1.0 + 1.0 + 1.0 / 3 + 1.0 / 4 + 1.0 / 4 + 1.0 / 4 + 1.0 + 1.0 / 2 + 1.0 / 4 ) / 9 ) ) {
                return false;
            }
            names.add( "E" );
            cov = cc.calculateCoverage( phylogenies, names, false );
            if ( !TestPccx.isEqual( cov.getScore(),
                                    ( 1.0 + 1.0 + 1.0 / 3 + 1.0 / 4 + 1.0 + 1.0 / 2 + 1.0 + 1.0 / 2 + 1.0 / 4 ) / 9 ) ) {
                return false;
            }
            names.add( "X" );
            cov = cc.calculateCoverage( phylogenies, names, false );
            if ( !TestPccx.isEqual( cov.getScore(),
                                    ( 1.0 + 1.0 + 1.0 / 3 + 1.0 / 3 + 1.0 + 1.0 / 2 + 1.0 + 1.0 / 2 + 1.0 ) / 9 ) ) {
                return false;
            }
            names.add( "C" );
            names.add( "C" );
            names.add( "C" );
            cov = cc.calculateCoverage( phylogenies, names, false );
            if ( !TestPccx.isEqual( cov.getScore(),
                                    ( 1.0 + 1.0 + 1.0 + 1.0 / 3 + 1.0 + 1.0 / 2 + 1.0 + 1.0 / 2 + 1.0 ) / 9 ) ) {
                return false;
            }
            names.add( "D" );
            cov = cc.calculateCoverage( phylogenies, names, false );
            if ( !TestPccx
                    .isEqual( cov.getScore(), ( 1.0 + 1.0 + 1.0 + 1.0 + 1.0 + 1.0 / 2 + 1.0 + 1.0 / 2 + 1.0 ) / 9 ) ) {
                return false;
            }
            names.add( "F" );
            cov = cc.calculateCoverage( phylogenies, names, false );
            if ( !TestPccx.isEqual( cov.getScore(), ( 1.0 + 1.0 + 1.0 + 1.0 + 1.0 + 1.0 + 1.0 + 1.0 / 2 + 1.0 ) / 9 ) ) {
                return false;
            }
            names.add( "H" );
            cov = cc.calculateCoverage( phylogenies, names, false );
            if ( !TestPccx.isEqual( cov.getScore(), ( 1.0 + 1.0 + 1.0 + 1.0 + 1.0 + 1.0 + 1.0 + 1.0 + 1.0 ) / 9 ) ) {
                return false;
            }
            final CoverageExtender ce = new BasicExternalNodeBasedCoverageExtender();
            List<String> l = ce
                    .find( phylogenies,
                           null,
                           0,
                           new ExternalNodeBasedCoverageMethodOptions( "org.forester.pccx.BranchCountingBasedScoringMethod" ),
                           null );
            if ( !l.get( 0 ).equals( "X" ) ) {
                return false;
            }
            if ( !l.get( 1 ).equals( "A" ) ) {
                return false;
            }
            if ( !l.get( 2 ).equals( "E" ) ) {
                return false;
            }
            if ( !l.get( 3 ).equals( "G" ) ) {
                return false;
            }
            if ( !l.get( 4 ).equals( "C" ) ) {
                return false;
            }
            if ( !l.get( 5 ).equals( "D" ) ) {
                return false;
            }
            if ( !l.get( 6 ).equals( "B" ) ) {
                return false;
            }
            if ( !l.get( 7 ).equals( "F" ) ) {
                return false;
            }
            if ( !l.get( 8 ).equals( "H" ) ) {
                return false;
            }
            final List<String> already_covered = new ArrayList<String>();
            already_covered.add( "A" );
            already_covered.add( "X" );
            already_covered.add( "H" );
            already_covered.add( "C" );
            l = ce
                    .find( phylogenies,
                           already_covered,
                           0,
                           new ExternalNodeBasedCoverageMethodOptions( "org.forester.pccx.BranchCountingBasedScoringMethod" ),
                           null );
            if ( !l.get( 0 ).equals( "E" ) ) {
                return false;
            }
            if ( !l.get( 1 ).equals( "D" ) ) {
                return false;
            }
            if ( !l.get( 2 ).equals( "B" ) ) {
                return false;
            }
            if ( !l.get( 3 ).equals( "F" ) ) {
                return false;
            }
            if ( !l.get( 4 ).equals( "G" ) ) {
                return false;
            }
            final String ps2 = "((((A:0.1,B:0.7):0.2,C:1.0):2.0,D:1.7):1.3,((E:0.3,F:0.4):1.1,(G:0.5,H:0.6):1.2):1.4,X:2.0)";
            final String ps3 = "((((A:0.1,B:0.1):0.2,C:1.0):2.0,D:1.7):1.3,((E:0.3,F:0.4):1.1,(G:0.5,H:0.6):1.2):1.4,X:2.0)";
            final String ps4 = "((((A:0.1,B:0.05):0.2,C:1.0):2.0,D:1.7):1.3,((E:0.3,F:0.4):1.1,(G:0.5,H:0.6):1.2):1.4,X:2.0)";
            final Phylogeny p2 = factory.create( ps2, new NHXParser() )[ 0 ];
            final Phylogeny p3 = factory.create( ps3, new NHXParser() )[ 0 ];
            final Phylogeny p4 = factory.create( ps4, new NHXParser() )[ 0 ];
            final List<Phylogeny> phylogenies2 = new ArrayList<Phylogeny>();
            final List<String> names2 = new ArrayList<String>();
            phylogenies2.add( p2 );
            phylogenies2.add( p3 );
            phylogenies2.add( p4 );
            names2.add( "A" );
            names2.add( "A" );
            final CoverageCalculationOptions options2 = new ExternalNodeBasedCoverageMethodOptions( "org.forester.pccx.BranchLengthBasedScoringMethod" );
            final CoverageCalculator cc2 = CoverageCalculator.getInstance( new ExternalNodeBasedCoverageMethod(),
                                                                           options2 );
            Coverage cov2 = cc2.calculateCoverage( phylogenies2, names2, false );
            final double nf = 1 / ( 1 / 0.1 + 1 / 0.7 + 1 / 1.0 + 1 / 1.7 + 1 / 0.3 + 1 / 0.4 + 1 / 0.5 + 1 / 0.6 + 1 / 2.0 );
            if ( !TestPccx.isEqual( cov2.getScore(), ( 1 / 0.1 + ( 1 / 0.8 + 1 / 0.2 + 1 / 0.15 ) / 3 + 1 / 1.3 + 1
                    / 4.0 + 1 / 6.4 + 1 / 6.5 + 1 / 6.7 + 1 / 6.8 + 1 / 5.6 )
                    * nf ) ) {
                return false;
            }
            names2.add( "C" );
            cov2 = cc2.calculateCoverage( phylogenies2, names2, false );
            if ( !TestPccx.isEqual( cov2.getScore(), ( 1 / 0.1 + ( 1 / 0.8 + 1 / 0.2 + 1 / 0.15 ) / 3 + 1 / 1.0 + 1
                    / 4.0 + 1 / 6.4 + 1 / 6.5 + 1 / 6.7 + 1 / 6.8 + 1 / 5.6 )
                    * nf ) ) {
                return false;
            }
            names2.add( "E" );
            cov2 = cc2.calculateCoverage( phylogenies2, names2, false );
            if ( !TestPccx.isEqual( cov2.getScore(), ( 1 / 0.1 + ( 1 / 0.8 + 1 / 0.2 + 1 / 0.15 ) / 3 + 1 / 1.0 + +1
                    / 4.0 + 1 / 0.3 + 1 / 0.7 + 1 / 3.1 + 1 / 3.2 + 1 / 4.8 )
                    * nf ) ) {
                return false;
            }
            final CoverageCalculationOptions options_log = new ExternalNodeBasedCoverageMethodOptions( "org.forester.pccx.LogBranchLengthBasedScoringMethod" );
            final CoverageCalculator cclog = CoverageCalculator.getInstance( new ExternalNodeBasedCoverageMethod(),
                                                                             options_log );
            final Coverage cov_log = cclog.calculateCoverage( phylogenies2, names2, false );
            if ( !TestPccx.isEqual( cov_log.getScore(), 0.8534252108361485 ) ) {
                return false;
            }
            final String ps10 = "((((A:0.1,B:0.7):0.2,C:1.0):2.0,D:1.7):1.3,((E:0.3,F:0.4):1.1,(G:0.5,H:0.6):1.2):1.4,((((I:0.1,J:0.7):0.2,K:1.0):2.0,L:1.7):1.3,((M:0.3,N:0.4,O:0.1,P:0.2):1.1,(Q:0.5,R:0.6):1.2):1.4,S:2.0):2.0)";
            final Phylogeny p10 = factory.create( ps10, new NHXParser() )[ 0 ];
            final List<Phylogeny> phylogenies10 = new ArrayList<Phylogeny>();
            final List<String> names10 = new ArrayList<String>();
            phylogenies10.add( p10 );
            names10.add( "A" );
            names10.add( "B" );
            names10.add( "N" );
            names10.add( "O" );
            final CoverageCalculationOptions options10 = new ExternalNodeBasedCoverageMethodOptions( "org.forester.pccx.BranchCountingBasedScoringMethod" );
            final CoverageCalculator cc10 = CoverageCalculator.getInstance( new ExternalNodeBasedCoverageMethod(),
                                                                            options10 );
            cc10.calculateCoverage( phylogenies10, names10, true );
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
}
