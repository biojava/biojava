// $Id: TestGSDI.java,v 1.13 2009/10/26 23:29:39 cmzmasek Exp $
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

package org.forester.sdi;

import java.io.IOException;

import org.forester.development.DevelopmentTools;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;

public final class TestGSDI {

    private final static Phylogeny createPhylogeny( final String nhx ) throws IOException {
        final Phylogeny p = ParserBasedPhylogenyFactory.getInstance().create( nhx, new NHXParser() )[ 0 ];
        p.setRooted( true );
        return p;
    }

    private final static Event getEvent( final Phylogeny p, final String n1, final String n2 ) {
        return PhylogenyMethods.getInstance().getLCA( p.getNode( n1 ), p.getNode( n2 ) ).getNodeData().getEvent();
    }

    public static boolean test() {
        if ( !TestGSDI.testGSDI_general() ) {
            return false;
        }
        if ( !TestGSDI.testGSDI_against_binary_gene_tree() ) {
            return false;
        }
        return true;
    }

    private static boolean testGSDI_against_binary_gene_tree() {
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final String multi_species_2_str = "(((((([&&NHX:S=1],[&&NHX:S=2]),"
                    + "([&&NHX:S=3],[&&NHX:S=4],[&&NHX:S=5])),"
                    + "([&&NHX:S=6],[&&NHX:S=7],[&&NHX:S=8],[&&NHX:S=9])),"
                    + "([&&NHX:S=10],[&&NHX:S=11])),"
                    + "([&&NHX:S=12],[&&NHX:S=13],[&&NHX:S=14])),"
                    + "([&&NHX:S=15],([&&NHX:S=16],[&&NHX:S=17]),([&&NHX:S=18],[&&NHX:S=19],[&&NHX:S=20]),([&&NHX:S=21],[&&NHX:S=22],[&&NHX:S=23],[&&NHX:S=24])));";
            final String gene_2_1_str = "(((((([&&NHX:S=1],[&&NHX:S=2])1_2,([&&NHX:S=3],[&&NHX:S=4])),"
                    + "([&&NHX:S=6],[&&NHX:S=7])6_7_8_9)1_9,([&&NHX:S=10],[&&NHX:S=11])),"
                    + "([&&NHX:S=12],[&&NHX:S=13])12_13_14)1_14,"
                    + "([&&NHX:S=15],([&&NHX:S=21],[&&NHX:S=24])21_22_23_24)15_24);";
            final Phylogeny multi_species_2 = factory.create( multi_species_2_str, new NHXParser() )[ 0 ];
            final Phylogeny gene_2_1 = factory.create( gene_2_1_str, new NHXParser() )[ 0 ];
            multi_species_2.setRooted( true );
            gene_2_1.setRooted( true );
            final GSDI sdi = new GSDI( gene_2_1, multi_species_2, false );
            if ( sdi.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi.getDuplicationsSum() != 0 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testGSDI_general() {
        try {
            final PhylogenyMethods pm = PhylogenyMethods.getInstance();
            final String s1_ = "((([&&NHX:S=A2],[&&NHX:S=A1]),[&&NHX:S=B],[&&NHX:S=C]),[&&NHX:S=D])";
            final Phylogeny s1 = ParserBasedPhylogenyFactory.getInstance().create( s1_, new NHXParser() )[ 0 ];
            s1.setRooted( true );
            final Phylogeny g1 = TestGSDI
                    .createPhylogeny( "((((B[&&NHX:S=B],A1[&&NHX:S=A1]),C[&&NHX:S=C]),A2[&&NHX:S=A2]),D[&&NHX:S=D])" );
            final GSDI sdi1 = new GSDI( g1, s1, false );
            if ( sdi1.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( !pm.getLCA( g1.getNode( "B" ), g1.getNode( "A1" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            if ( !pm.getLCA( g1.getNode( "C" ), g1.getNode( "A1" ) ).getNodeData().getEvent()
                    .isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !( pm.getLCA( g1.getNode( "A2" ), g1.getNode( "A1" ) ).getNodeData().getEvent().isDuplication() ) ) {
                return false;
            }
            if ( !pm.getLCA( g1.getNode( "D" ), g1.getNode( "A1" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            final Phylogeny g2 = TestGSDI
                    .createPhylogeny( "((((A2[&&NHX:S=A2],A1[&&NHX:S=A1]),B[&&NHX:S=B]),C[&&NHX:S=C]),D[&&NHX:S=D])" );
            final GSDI sdi2 = new GSDI( g2, s1, false );
            if ( sdi2.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( !pm.getLCA( g2.getNode( "A1" ), g2.getNode( "A2" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            if ( !pm.getLCA( g2.getNode( "A1" ), g2.getNode( "B" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            if ( !pm.getLCA( g2.getNode( "A1" ), g2.getNode( "C" ) ).getNodeData().getEvent()
                    .isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !pm.getLCA( g2.getNode( "A1" ), g2.getNode( "D" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            final Phylogeny g3 = TestGSDI
                    .createPhylogeny( "((((A2[&&NHX:S=A2],A1[&&NHX:S=A1]),C[&&NHX:S=C]),B[&&NHX:S=B]),D[&&NHX:S=D])" );
            final GSDI sdi3 = new GSDI( g3, s1, false );
            if ( sdi3.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( !pm.getLCA( g3.getNode( "A1" ), g3.getNode( "A2" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            if ( !pm.getLCA( g3.getNode( "A1" ), g3.getNode( "C" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            if ( !pm.getLCA( g3.getNode( "A1" ), g3.getNode( "B" ) ).getNodeData().getEvent()
                    .isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !pm.getLCA( g3.getNode( "A1" ), g3.getNode( "D" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            final Phylogeny g4 = TestGSDI
                    .createPhylogeny( "(((B[&&NHX:S=B],C1[&&NHX:S=C]),C2[&&NHX:S=C]),D[&&NHX:S=D])" );
            final GSDI sdi4 = new GSDI( g4, s1, false );
            if ( sdi4.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( !pm.getLCA( g4.getNode( "B" ), g4.getNode( "C1" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            if ( !pm.getLCA( g4.getNode( "B" ), g4.getNode( "C2" ) ).getNodeData().getEvent().isDuplication() ) {
                return false;
            }
            if ( !pm.getLCA( g4.getNode( "B" ), g4.getNode( "D" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            final Phylogeny g5 = TestGSDI
                    .createPhylogeny( "(((D1[&&NHX:S=D],A1[&&NHX:S=A1]),B[&&NHX:S=B]),((D2[&&NHX:S=D],D3[&&NHX:S=D]),C[&&NHX:S=C]))" );
            final GSDI sdi5 = new GSDI( g5, s1, false );
            if ( sdi5.getDuplicationsSum() != 3 ) {
                return false;
            }
            if ( !pm.getLCA( g5.getNode( "D1" ), g5.getNode( "A1" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            if ( !pm.getLCA( g5.getNode( "D1" ), g5.getNode( "B" ) ).getNodeData().getEvent().isDuplication() ) {
                return false;
            }
            if ( !pm.getLCA( g5.getNode( "D1" ), g5.getNode( "D2" ) ).getNodeData().getEvent().isDuplication() ) {
                return false;
            }
            if ( !pm.getLCA( g5.getNode( "D2" ), g5.getNode( "D3" ) ).getNodeData().getEvent().isDuplication() ) {
                return false;
            }
            if ( !pm.getLCA( g5.getNode( "C" ), g5.getNode( "D3" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            final Phylogeny species7 = TestGSDI.createPhylogeny( "(((((((([&&NHX:S=a1],[&&NHX:S=a2]),"
                    + "([&&NHX:S=b1],[&&NHX:S=b2])),[&&NHX:S=x]),(([&&NHX:S=m1],[&&NHX:S=m2]),"
                    + "([&&NHX:S=n1],[&&NHX:S=n2]))),(([&&NHX:S=i1],[&&NHX:S=i2]),"
                    + "([&&NHX:S=j1],[&&NHX:S=j2]))),(([&&NHX:S=e1],[&&NHX:S=e2]),"
                    + "([&&NHX:S=f1],[&&NHX:S=f2]))),[&&NHX:S=y]),[&&NHX:S=z])" );
            final Phylogeny gene7_2 = TestGSDI
                    .createPhylogeny( "(((((((((a1[&&NHX:S=a1],a2[&&NHX:S=a2]),b1[&&NHX:S=b1]),x[&&NHX:S=x]),m1[&&NHX:S=m1]),i1[&&NHX:S=i1]),j2[&&NHX:S=j2]),e1[&&NHX:S=e1]),y[&&NHX:S=y]),z[&&NHX:S=z])" );
            gene7_2.setRooted( true );
            final GSDI sdi7_2 = new GSDI( gene7_2, species7, false );
            if ( sdi7_2.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( gene7_2, "a1", "a2" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( gene7_2, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( gene7_2, "a1", "x" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( gene7_2, "a1", "m1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( gene7_2, "a1", "i1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( gene7_2, "a1", "j2" ).isDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( gene7_2, "a1", "e1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( gene7_2, "a1", "y" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( gene7_2, "a1", "z" ).isSpeciation() ) {
                return false;
            }
            final String s2_ = "((" + "([&&NHX:S=a1],[&&NHX:S=a2],[&&NHX:S=a3],[&&NHX:S=a4]),"
                    + "([&&NHX:S=b1],[&&NHX:S=b2],[&&NHX:S=b3],[&&NHX:S=b4]),"
                    + "([&&NHX:S=c1],[&&NHX:S=c2],[&&NHX:S=c3],[&&NHX:S=c4]),"
                    + "([&&NHX:S=d1],[&&NHX:S=d2],[&&NHX:S=d3],[&&NHX:S=d4])),("
                    + "([&&NHX:S=e1],[&&NHX:S=e2],[&&NHX:S=e3],[&&NHX:S=e4]),"
                    + "([&&NHX:S=f1],[&&NHX:S=f2],[&&NHX:S=f3],[&&NHX:S=f4]),"
                    + "([&&NHX:S=g1],[&&NHX:S=g2],[&&NHX:S=g3],[&&NHX:S=g4]),"
                    + "([&&NHX:S=h1],[&&NHX:S=h2],[&&NHX:S=h3],[&&NHX:S=h4])),("
                    + "([&&NHX:S=i1],[&&NHX:S=i2],[&&NHX:S=i3],[&&NHX:S=i4]),"
                    + "([&&NHX:S=j1],[&&NHX:S=j2],[&&NHX:S=j3],[&&NHX:S=j4]),"
                    + "([&&NHX:S=k1],[&&NHX:S=k2],[&&NHX:S=k3],[&&NHX:S=k4]),"
                    + "([&&NHX:S=l1],[&&NHX:S=l2],[&&NHX:S=l3],[&&NHX:S=l4])),("
                    + "([&&NHX:S=m1],[&&NHX:S=m2],[&&NHX:S=m3],[&&NHX:S=m4]),"
                    + "([&&NHX:S=n1],[&&NHX:S=n2],[&&NHX:S=n3],[&&NHX:S=n4]),"
                    + "([&&NHX:S=o1],[&&NHX:S=o2],[&&NHX:S=o3],[&&NHX:S=o4]),"
                    + "([&&NHX:S=p1],[&&NHX:S=p2],[&&NHX:S=p3],[&&NHX:S=p4])"
                    + "),[&&NHX:S=x],[&&NHX:S=y],[&&NHX:S=z])";
            final Phylogeny s2 = ParserBasedPhylogenyFactory.getInstance().create( s2_, new NHXParser() )[ 0 ];
            s2.setRooted( true );
            final Phylogeny g2_0 = TestGSDI.createPhylogeny( "(m1[&&NHX:S=m1],m3[&&NHX:S=m3])" );
            final GSDI sdi2_0 = new GSDI( g2_0, s2, false );
            if ( sdi2_0.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_0.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_0.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !pm.getLCA( g2_0.getNode( "m1" ), g2_0.getNode( "m3" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            final Phylogeny g2_1 = TestGSDI.createPhylogeny( "(e2[&&NHX:S=e2],h2[&&NHX:S=h2])" );
            final GSDI sdi2_1 = new GSDI( g2_1, s2, false );
            if ( sdi2_1.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_1.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_1.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !pm.getLCA( g2_1.getNode( "e2" ), g2_1.getNode( "h2" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            final Phylogeny g2_2 = TestGSDI.createPhylogeny( "(e2[&&NHX:S=e2],p4[&&NHX:S=p4])" );
            final GSDI sdi2_2 = new GSDI( g2_2, s2, false );
            if ( sdi2_2.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_2.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_2.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !pm.getLCA( g2_2.getNode( "e2" ), g2_2.getNode( "p4" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            final Phylogeny g2_3 = TestGSDI.createPhylogeny( "(e2a[&&NHX:S=e2],e2b[&&NHX:S=e2])" );
            final GSDI sdi2_3 = new GSDI( g2_3, s2, false );
            if ( sdi2_3.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_3.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_3.getSpeciationsSum() != 0 ) {
                return false;
            }
            if ( !pm.getLCA( g2_3.getNode( "e2a" ), g2_3.getNode( "e2b" ) ).getNodeData().getEvent().isDuplication() ) {
                return false;
            }
            final Phylogeny g2_4 = TestGSDI.createPhylogeny( "((j1[&&NHX:S=j1],j4[&&NHX:S=j4]),i3[&&NHX:S=i3])" );
            final GSDI sdi2_4 = new GSDI( g2_4, s2, false );
            if ( sdi2_4.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_4.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_4.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !pm.getLCA( g2_4.getNode( "j1" ), g2_4.getNode( "j4" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            if ( !pm.getLCA( g2_4.getNode( "j1" ), g2_4.getNode( "i3" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            final Phylogeny g2_5 = TestGSDI.createPhylogeny( "((j1[&&NHX:S=j1],j4[&&NHX:S=j4]),f3[&&NHX:S=f3])" );
            final GSDI sdi2_5 = new GSDI( g2_5, s2, false );
            if ( sdi2_5.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_5.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_5.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !pm.getLCA( g2_5.getNode( "j1" ), g2_5.getNode( "j4" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            if ( !pm.getLCA( g2_5.getNode( "j1" ), g2_5.getNode( "f3" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            final Phylogeny g2_6 = TestGSDI.createPhylogeny( "((j3[&&NHX:S=j3],i4[&&NHX:S=i4]),f3[&&NHX:S=f3])" );
            final GSDI sdi2_6 = new GSDI( g2_6, s2, false );
            if ( sdi2_6.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_6.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_6.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !pm.getLCA( g2_6.getNode( "j3" ), g2_6.getNode( "i4" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            if ( !pm.getLCA( g2_6.getNode( "j3" ), g2_6.getNode( "f3" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            final Phylogeny g2_7 = TestGSDI.createPhylogeny( "((j1[&&NHX:S=j1],k1[&&NHX:S=k1]),i1[&&NHX:S=i1])" );
            final GSDI sdi2_7 = new GSDI( g2_7, s2, false );
            if ( sdi2_7.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_7.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_7.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !pm.getLCA( g2_7.getNode( "j1" ), g2_7.getNode( "k1" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            if ( !pm.getLCA( g2_7.getNode( "j1" ), g2_7.getNode( "i1" ) ).getNodeData().getEvent()
                    .isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_8 = TestGSDI.createPhylogeny( "(j1[&&NHX:S=j1],(k1[&&NHX:S=k1],i1[&&NHX:S=i1]))" );
            final GSDI sdi2_8 = new GSDI( g2_8, s2, false );
            if ( sdi2_8.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_8.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_8.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !pm.getLCA( g2_8.getNode( "j1" ), g2_8.getNode( "k1" ) ).getNodeData().getEvent()
                    .isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !pm.getLCA( g2_8.getNode( "k1" ), g2_8.getNode( "i1" ) ).getNodeData().getEvent().isSpeciation() ) {
                return false;
            }
            final Phylogeny g2_9 = TestGSDI.createPhylogeny( "((j1[&&NHX:S=j1],k4[&&NHX:S=k4]),f2[&&NHX:S=f2])" );
            final GSDI sdi2_9 = new GSDI( g2_9, s2, false );
            if ( sdi2_9.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_9.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_9.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_9, "j1", "k4" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_9, "j1", "f2" ).isSpeciation() ) {
                return false;
            }
            final Phylogeny g2_10 = TestGSDI.createPhylogeny( "((m1[&&NHX:S=m1],k4[&&NHX:S=k4]),f2[&&NHX:S=f2])" );
            final GSDI sdi2_10 = new GSDI( g2_10, s2, false );
            if ( sdi2_10.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_10.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_10.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_10, "m1", "k4" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_10, "m1", "f2" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_11 = TestGSDI.createPhylogeny( "((m1[&&NHX:S=m1],k4[&&NHX:S=k4]),x[&&NHX:S=x])" );
            final GSDI sdi2_11 = new GSDI( g2_11, s2, false );
            if ( sdi2_11.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_11.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_11.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_11, "m1", "k4" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_11, "m1", "x" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_12 = TestGSDI.createPhylogeny( "(m1[&&NHX:S=m1],(k4[&&NHX:S=k4],x[&&NHX:S=x]))" );
            final GSDI sdi2_12 = new GSDI( g2_12, s2, false );
            if ( sdi2_12.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_12.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_12.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_12, "x", "k4" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_12, "m1", "x" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_13 = TestGSDI.createPhylogeny( "(x[&&NHX:S=x],(y[&&NHX:S=y],z[&&NHX:S=z]))" );
            final GSDI sdi2_13 = new GSDI( g2_13, s2, false );
            if ( sdi2_13.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_13.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_13.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_13, "y", "z" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_13, "x", "z" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_14 = TestGSDI.createPhylogeny( "(a1_1[&&NHX:S=a1],(b1[&&NHX:S=b1],a1[&&NHX:S=a1]))" );
            final GSDI sdi2_14 = new GSDI( g2_14, s2, false );
            if ( sdi2_14.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_14.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_14.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_14, "b1", "a1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_14, "b1", "a1_1" ).isDuplication() ) {
                return false;
            }
            final Phylogeny g2_15 = TestGSDI.createPhylogeny( "(a2[&&NHX:S=a2],(b1[&&NHX:S=b1],a1[&&NHX:S=a1]))" );
            final GSDI sdi2_15 = new GSDI( g2_15, s2, false );
            if ( sdi2_15.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_15.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_15.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_15, "b1", "a1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_15, "b1", "a2" ).isDuplication() ) {
                return false;
            }
            final Phylogeny g2_16 = TestGSDI.createPhylogeny( "(n2[&&NHX:S=n2],(j3[&&NHX:S=j3],n1[&&NHX:S=n1]))" );
            final GSDI sdi2_16 = new GSDI( g2_16, s2, false );
            if ( sdi2_16.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_16.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_16.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_16, "j3", "n1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_16, "j3", "n2" ).isDuplication() ) {
                return false;
            }
            final Phylogeny g2_17 = TestGSDI.createPhylogeny( "(p4[&&NHX:S=p4],(j3[&&NHX:S=j3],n1[&&NHX:S=n1]))" );
            final GSDI sdi2_17 = new GSDI( g2_17, s2, false );
            if ( sdi2_17.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_17.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_17.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_17, "j3", "n1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_17, "j3", "p4" ).isDuplication() ) {
                return false;
            }
            final Phylogeny g2_18 = TestGSDI
                    .createPhylogeny( "((n11[&&NHX:S=n1],n12[&&NHX:S=n1]),(n13[&&NHX:S=n1],n14[&&NHX:S=n1]))" );
            final GSDI sdi2_18 = new GSDI( g2_18, s2, false );
            if ( sdi2_18.getDuplicationsSum() != 3 ) {
                return false;
            }
            if ( sdi2_18.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_18.getSpeciationsSum() != 0 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_18, "n11", "n12" ).isDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_18, "n13", "n14" ).isDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_18, "n11", "n13" ).isDuplication() ) {
                return false;
            }
            final Phylogeny g2_19 = TestGSDI
                    .createPhylogeny( "((n11[&&NHX:S=n1],n21[&&NHX:S=n2]),(n12[&&NHX:S=n1],n22[&&NHX:S=n2]))" );
            final GSDI sdi2_19 = new GSDI( g2_19, s2, false );
            if ( sdi2_19.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_19.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_19.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_19, "n11", "n21" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_19, "n12", "n22" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_19, "n11", "n12" ).isDuplication() ) {
                return false;
            }
            final Phylogeny g2_20 = TestGSDI
                    .createPhylogeny( "((n11[&&NHX:S=n1],n2[&&NHX:S=n2]),(n12[&&NHX:S=n1],n3[&&NHX:S=n3]))" );
            final GSDI sdi2_20 = new GSDI( g2_20, s2, false );
            if ( sdi2_20.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_20.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_20.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_20, "n11", "n2" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_20, "n12", "n3" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_20, "n11", "n12" ).isDuplication() ) {
                return false;
            }
            final Phylogeny g2_21 = TestGSDI
                    .createPhylogeny( "((n1[&&NHX:S=n1],n2[&&NHX:S=n2]),(n3[&&NHX:S=n3],a1[&&NHX:S=a1]))" );
            final GSDI sdi2_21 = new GSDI( g2_21, s2, false );
            if ( sdi2_21.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_21.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_21.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_21, "n1", "n2" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_21, "n3", "a1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_21, "n2", "a1" ).isDuplication() ) {
                return false;
            }
            final Phylogeny g2_22 = TestGSDI
                    .createPhylogeny( "((n1[&&NHX:S=n1],n2[&&NHX:S=n2]),(n3[&&NHX:S=n3],n4[&&NHX:S=n4]))" );
            final GSDI sdi2_22 = new GSDI( g2_22, s2, false );
            if ( sdi2_22.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_22.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_22.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_22, "n1", "n2" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_22, "n3", "n4" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_22, "n1", "n3" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_23 = TestGSDI
                    .createPhylogeny( "((a1[&&NHX:S=a1],b1[&&NHX:S=b1]),(c1[&&NHX:S=c1],d1[&&NHX:S=d1]))" );
            final GSDI sdi2_23 = new GSDI( g2_23, s2, false );
            if ( sdi2_23.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_23.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_23.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_23, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_23, "c1", "d1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_23, "a1", "c1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_24 = TestGSDI
                    .createPhylogeny( "((a1[&&NHX:S=a1],e1[&&NHX:S=e1]),(i1[&&NHX:S=i1],m1[&&NHX:S=m1]))" );
            final GSDI sdi2_24 = new GSDI( g2_24, s2, false );
            if ( sdi2_24.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_24.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_24.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_24, "a1", "e1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_24, "i1", "m1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_24, "a1", "i1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_25 = TestGSDI
                    .createPhylogeny( "((a1[&&NHX:S=a1],a4[&&NHX:S=a4]),(b1[&&NHX:S=b1],c1[&&NHX:S=c1]))" );
            final GSDI sdi2_25 = new GSDI( g2_25, s2, false );
            if ( sdi2_25.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_25.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_25.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_25, "a1", "a4" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_25, "b1", "c1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_25, "a1", "b1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_26 = TestGSDI
                    .createPhylogeny( "(((a1[&&NHX:S=a1],a4[&&NHX:S=a4]),b1[&&NHX:S=b1]),e1[&&NHX:S=e1])" );
            final GSDI sdi2_26 = new GSDI( g2_26, s2, false );
            if ( sdi2_26.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_26.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_26.getSpeciationsSum() != 3 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_26, "a1", "a4" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_26, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_26, "a1", "e1" ).isSpeciation() ) {
                return false;
            }
            final Phylogeny g2_27 = TestGSDI
                    .createPhylogeny( "(((a1[&&NHX:S=a1],a4[&&NHX:S=a4]),b1[&&NHX:S=b1]),c1[&&NHX:S=c1])" );
            final GSDI sdi2_27 = new GSDI( g2_27, s2, false );
            if ( sdi2_27.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_27.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_27.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_27, "a1", "a4" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_27, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_27, "a1", "c1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_28 = TestGSDI
                    .createPhylogeny( "(((a1[&&NHX:S=a1],b1[&&NHX:S=b1]),c1[&&NHX:S=c1]),e1[&&NHX:S=e1])" );
            final GSDI sdi2_28 = new GSDI( g2_28, s2, false );
            if ( sdi2_28.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_28.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_28.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_28, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_28, "a1", "c1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_28, "a1", "e1" ).isSpeciation() ) {
                return false;
            }
            final Phylogeny g2_29 = TestGSDI
                    .createPhylogeny( "(((a1[&&NHX:S=a1],b1[&&NHX:S=b1]),c1[&&NHX:S=c1]),d1[&&NHX:S=d1])" );
            final GSDI sdi2_29 = new GSDI( g2_29, s2, false );
            if ( sdi2_29.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_29.getSpeciationOrDuplicationEventsSum() != 2 ) {
                return false;
            }
            if ( sdi2_29.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_29, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_29, "a1", "c1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_29, "a1", "d1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_30 = TestGSDI
                    .createPhylogeny( "(((a1[&&NHX:S=a1],b1[&&NHX:S=b1]),c1[&&NHX:S=c1]),a2[&&NHX:S=a2])" );
            final GSDI sdi2_30 = new GSDI( g2_30, s2, false );
            if ( sdi2_30.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_30.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_30.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_30, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_30, "a1", "c1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_30, "a1", "a2" ).isDuplication() ) {
                return false;
            }
            final Phylogeny g2_31 = TestGSDI
                    .createPhylogeny( "(((a1[&&NHX:S=a1],b1[&&NHX:S=b1]),c1[&&NHX:S=c1]),c2[&&NHX:S=c2])" );
            final GSDI sdi2_31 = new GSDI( g2_31, s2, false );
            if ( sdi2_31.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_31.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_31.getSpeciationsSum() != 1 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_31, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_31, "a1", "c1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_31, "a1", "c2" ).isDuplication() ) {
                return false;
            }
            final Phylogeny g2_32 = TestGSDI
                    .createPhylogeny( "((((((((((a1[&&NHX:S=a1],a2[&&NHX:S=a2]),b1[&&NHX:S=b1]),c1[&&NHX:S=c1]),d1[&&NHX:S=d1]),x[&&NHX:S=x]),p1[&&NHX:S=p1]),i1[&&NHX:S=i1]),e1[&&NHX:S=e1]),y[&&NHX:S=y]),z[&&NHX:S=z])" );
            final GSDI sdi2_32 = new GSDI( g2_32, s2, false );
            if ( sdi2_32.getDuplicationsSum() != 0 ) {
                return false;
            }
            if ( sdi2_32.getSpeciationOrDuplicationEventsSum() != 7 ) {
                return false;
            }
            if ( sdi2_32.getSpeciationsSum() != 3 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_32, "a1", "a2" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_32, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_32, "a1", "c1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_32, "a1", "d1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_32, "a1", "x" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_32, "a1", "p1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_32, "a1", "i1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_32, "a1", "e1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_32, "a1", "y" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_32, "a1", "z" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_33 = TestGSDI
                    .createPhylogeny( "(((((((((((a1[&&NHX:S=a1],a2[&&NHX:S=a2]),b1[&&NHX:S=b1]),c1[&&NHX:S=c1]),d1[&&NHX:S=d1]),x[&&NHX:S=x]),p1[&&NHX:S=p1]),i1[&&NHX:S=i1]),k2[&&NHX:S=k2]),e1[&&NHX:S=e1]),y[&&NHX:S=y]),z[&&NHX:S=z])" );
            final GSDI sdi2_33 = new GSDI( g2_33, s2, false );
            if ( sdi2_33.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_33.getSpeciationOrDuplicationEventsSum() != 7 ) {
                return false;
            }
            if ( sdi2_33.getSpeciationsSum() != 3 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_33, "a1", "a2" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_33, "a1", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_33, "a1", "c1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_33, "a1", "d1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_33, "a1", "x" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_33, "a1", "p1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_33, "a1", "i1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_33, "a1", "k2" ).isDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_33, "a1", "e1" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_33, "a1", "y" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_33, "a1", "z" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_34 = TestGSDI
                    .createPhylogeny( "(((n1_0[&&NHX:S=n1],n2_0[&&NHX:S=n2]),(n1_1[&&NHX:S=n1],n3_0[&&NHX:S=n3])),n4_0[&&NHX:S=n4])" );
            final GSDI sdi2_34 = new GSDI( g2_34, s2, false );
            if ( sdi2_34.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_34.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_34.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_34, "n1_0", "n2_0" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_34, "n1_1", "n3_0" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_34, "n1_0", "n1_1" ).isDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_34, "n1_0", "n4_0" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_35 = TestGSDI
                    .createPhylogeny( "((((n1_0[&&NHX:S=n1],n2_0[&&NHX:S=n2]),(n1_1[&&NHX:S=n1],n3_0[&&NHX:S=n3])),n4_0[&&NHX:S=n4]),a1_0[&&NHX:S=a1])" );
            final GSDI sdi2_35 = new GSDI( g2_35, s2, false );
            if ( sdi2_35.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_35.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_35.getSpeciationsSum() != 3 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_35, "n1_0", "n2_0" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_35, "n1_1", "n3_0" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_35, "n1_0", "n1_1" ).isDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_35, "n1_0", "n4_0" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_35, "n1_0", "a1_0" ).isSpeciation() ) {
                return false;
            }
            final Phylogeny g2_36 = TestGSDI
                    .createPhylogeny( "(((a1_0[&&NHX:S=a1],b1_0[&&NHX:S=b1]),(a1_1[&&NHX:S=a1],c1_0[&&NHX:S=c1])),d1_0[&&NHX:S=d1])" );
            final GSDI sdi2_36 = new GSDI( g2_36, s2, false );
            if ( sdi2_36.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_36.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_36.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_36, "a1_0", "b1_0" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_36, "a1_1", "c1_0" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_36, "a1_0", "c1_0" ).isDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_36, "a1_0", "d1_0" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_37 = TestGSDI
                    .createPhylogeny( "(((a1_0[&&NHX:S=a1],b1_0[&&NHX:S=b1]),(a2_0[&&NHX:S=a2],c1_0[&&NHX:S=c1])),d1_0[&&NHX:S=d1])" );
            final GSDI sdi2_37 = new GSDI( g2_37, s2, false );
            if ( sdi2_37.getDuplicationsSum() != 1 ) {
                return false;
            }
            if ( sdi2_37.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_37.getSpeciationsSum() != 2 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_37, "a1_0", "b1_0" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_37, "a2_0", "c1_0" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_37, "a1_0", "c1_0" ).isDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_37, "a1_0", "d1_0" ).isSpeciationOrDuplication() ) {
                return false;
            }
            final Phylogeny g2_38 = TestGSDI
                    .createPhylogeny( "(((([&&NHX:S=n1],[&&NHX:S=n1]),([&&NHX:S=n1],[&&NHX:S=n1])),[&&NHX:S=n1]),[&&NHX:S=n1])" );
            final GSDI sdi2_38 = new GSDI( g2_38, s2, false );
            if ( sdi2_38.getDuplicationsSum() != 5 ) {
                return false;
            }
            if ( sdi2_38.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_38.getSpeciationsSum() != 0 ) {
                return false;
            }
            final Phylogeny g2_100 = TestGSDI
                    .createPhylogeny( "(((e1[&&NHX:S=e1],f2[&&NHX:S=f2]),(d3[&&NHX:S=d3],g4[&&NHX:S=g4])),(((a1[&&NHX:S=a1],h2[&&NHX:S=h2]),c3[&&NHX:S=c3]),(i4[&&NHX:S=i4],b1[&&NHX:S=b1])))" );
            final GSDI sdi2_100 = new GSDI( g2_100, s2, false );
            if ( sdi2_100.getDuplicationsSum() != 4 ) {
                return false;
            }
            if ( sdi2_100.getSpeciationOrDuplicationEventsSum() != 0 ) {
                return false;
            }
            if ( sdi2_100.getSpeciationsSum() != 4 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_100, "e1", "f2" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_100, "d3", "g4" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_100, "e1", "d3" ).isDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_100, "a1", "h2" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_100, "a1", "c3" ).isDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_100, "i4", "b1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_100, "a1", "i4" ).isDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_100, "e1", "a1" ).isDuplication() ) {
                return false;
            }
            final Phylogeny g2_101 = TestGSDI
                    .createPhylogeny( "(((e1[&&NHX:S=e1],f2[&&NHX:S=f2]),(d3[&&NHX:S=d3],g4[&&NHX:S=g4])),(((a1[&&NHX:S=a1],b2[&&NHX:S=b2]),c3[&&NHX:S=c3]),(i4[&&NHX:S=i4],j1[&&NHX:S=j1])))" );
            final GSDI sdi2_101 = new GSDI( g2_101, s2, false );
            if ( sdi2_101.getDuplicationsSum() != 2 ) {
                return false;
            }
            if ( sdi2_101.getSpeciationOrDuplicationEventsSum() != 1 ) {
                return false;
            }
            if ( sdi2_101.getSpeciationsSum() != 5 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_101, "e1", "f2" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_101, "d3", "g4" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_101, "e1", "d3" ).isDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_101, "a1", "b2" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_101, "a1", "c3" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_101, "i4", "j1" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_101, "a1", "i4" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g2_101, "e1", "a1" ).isDuplication() ) {
                return false;
            }
            final Phylogeny s_7_4 = DevelopmentTools.createBalancedPhylogeny( 7, 4 );
            DevelopmentTools.numberSpeciesInOrder( s_7_4 );
            final Phylogeny g_7_4_1 = TestGSDI
                    .createPhylogeny( "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("
                            + "1[&&NHX:S=1],2[&&NHX:S=2]),3[&&NHX:S=3]),4[&&NHX:S=4]),5[&&NHX:S=5]),"
                            + "6[&&NHX:S=6]),7[&&NHX:S=7]),8[&&NHX:S=8]),9[&&NHX:S=9]),10[&&NHX:S=10]),11[&&NHX:S=11]),"
                            + "12[&&NHX:S=12]),13[&&NHX:S=13]),14[&&NHX:S=14]),15[&&NHX:S=15]),16[&&NHX:S=16]),17[&&NHX:S=17]),"
                            + "18[&&NHX:S=18]),19[&&NHX:S=19]),20[&&NHX:S=20]),21[&&NHX:S=21]),22[&&NHX:S=22]),23[&&NHX:S=23]),"
                            + "24[&&NHX:S=24]),25[&&NHX:S=25]),26[&&NHX:S=26]),27[&&NHX:S=27]),28[&&NHX:S=28]),29[&&NHX:S=29]),"
                            + "30[&&NHX:S=30]),31[&&NHX:S=31]),32[&&NHX:S=32]),33[&&NHX:S=33]),34[&&NHX:S=34]),35[&&NHX:S=35]),"
                            + "36[&&NHX:S=36]),37[&&NHX:S=37]),38[&&NHX:S=38]),39[&&NHX:S=39]),40[&&NHX:S=40]),41[&&NHX:S=41]),"
                            + "42[&&NHX:S=42]),43[&&NHX:S=43]),44[&&NHX:S=44]),45[&&NHX:S=45]),46[&&NHX:S=46]),47[&&NHX:S=47]),"
                            + "48[&&NHX:S=48]),49[&&NHX:S=49]),50[&&NHX:S=50]),51[&&NHX:S=51]),52[&&NHX:S=52]),53[&&NHX:S=53]),"
                            + "54[&&NHX:S=54]),55[&&NHX:S=55]),56[&&NHX:S=56]),57[&&NHX:S=57]),58[&&NHX:S=58]),59[&&NHX:S=59]),"
                            + "60[&&NHX:S=60]),61[&&NHX:S=61]),62[&&NHX:S=62]),63[&&NHX:S=63]),64[&&NHX:S=64]),65[&&NHX:S=65])" );
            final GSDI sdi7_4_1 = new GSDI( g_7_4_1, s_7_4, false );
            if ( sdi7_4_1.getDuplicationsSum() != 54 ) {
                return false;
            }
            if ( sdi7_4_1.getSpeciationOrDuplicationEventsSum() != 6 ) {
                return false;
            }
            if ( sdi7_4_1.getSpeciationsSum() != 4 ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g_7_4_1, "1", "2" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g_7_4_1, "1", "3" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g_7_4_1, "1", "4" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g_7_4_1, "1", "5" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g_7_4_1, "1", "6" ).isDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g_7_4_1, "1", "9" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g_7_4_1, "1", "13" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g_7_4_1, "1", "17" ).isSpeciation() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g_7_4_1, "1", "33" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g_7_4_1, "1", "49" ).isSpeciationOrDuplication() ) {
                return false;
            }
            if ( !TestGSDI.getEvent( g_7_4_1, "1", "65" ).isSpeciation() ) {
                return false;
            }
            final Phylogeny g_7_4_2 = TestGSDI
                    .createPhylogeny( "((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("
                            + "1[&&NHX:S=1],2[&&NHX:S=2]),3[&&NHX:S=3]),4[&&NHX:S=4]),5[&&NHX:S=5]),"
                            + "6[&&NHX:S=6]),7[&&NHX:S=7]),8[&&NHX:S=8]),9[&&NHX:S=9]),10[&&NHX:S=10]),11[&&NHX:S=11]),"
                            + "12[&&NHX:S=12]),13[&&NHX:S=13]),14[&&NHX:S=14]),15[&&NHX:S=15]),16[&&NHX:S=16]),17[&&NHX:S=17]),"
                            + "18[&&NHX:S=18]),19[&&NHX:S=19]),20[&&NHX:S=20]),21[&&NHX:S=21]),22[&&NHX:S=22]),23[&&NHX:S=23]),"
                            + "24[&&NHX:S=24]),25[&&NHX:S=25]),26[&&NHX:S=26]),27[&&NHX:S=27]),28[&&NHX:S=28]),29[&&NHX:S=29]),"
                            + "30[&&NHX:S=30]),31[&&NHX:S=31]),32[&&NHX:S=32]),33[&&NHX:S=33]),34[&&NHX:S=34]),35[&&NHX:S=35]),"
                            + "36[&&NHX:S=36]),37[&&NHX:S=37]),38[&&NHX:S=38]),39[&&NHX:S=39]),40[&&NHX:S=40]),41[&&NHX:S=41]),"
                            + "42[&&NHX:S=42]),43[&&NHX:S=43]),44[&&NHX:S=44]),45[&&NHX:S=45]),46[&&NHX:S=46]),47[&&NHX:S=47]),"
                            + "48[&&NHX:S=48]),49[&&NHX:S=49]),50[&&NHX:S=50]),51[&&NHX:S=51]),52[&&NHX:S=52]),53[&&NHX:S=53]),"
                            + "54[&&NHX:S=54]),55[&&NHX:S=55]),56[&&NHX:S=56]),57[&&NHX:S=57]),58[&&NHX:S=58]),59[&&NHX:S=59]),"
                            + "60[&&NHX:S=60]),61[&&NHX:S=61]),62[&&NHX:S=62]),63[&&NHX:S=63]),64[&&NHX:S=64]),65[&&NHX:S=65]),"
                            + "66[&&NHX:S=66]),257[&&NHX:S=257]),258[&&NHX:S=258]),513[&&NHX:S=513]),514[&&NHX:S=514]),769[&&NHX:S=769]),770[&&NHX:S=770])" );
            final GSDI sdi7_4_2 = new GSDI( g_7_4_2, s_7_4, false );
            if ( sdi7_4_2.getDuplicationsSum() != 58 ) {
                return false;
            }
            if ( sdi7_4_2.getSpeciationOrDuplicationEventsSum() != 8 ) {
                return false;
            }
            if ( sdi7_4_2.getSpeciationsSum() != 5 ) {
                return false;
            }
            // final String g2_0_ =
            // "(([&&NHX:S=a1],[&&NHX:S=a2]),([&&NHX:S=o2],[&&NHX:S=o4]))";
            // final Phylogeny g2_0 = factory.create( g2_0_, new NHXParser() )[
            // 0 ];
            // g2_0.setRooted( true );
            // final GSDI sdi2_0 = new GSDI( g2_0, s2, false );
            // if ( sdi2_0.getDuplicationsSum() != 0 ) {
            // return false;
            // }
            // final String g2_1_= "";
            // final Phylogeny g2_1 = factory.create( g2_1_, new NHXParser() )[
            // 0 ];
            // g2_1.setRooted( true );
            // final GSDI sdi2_1 = new GSDI( g2_1, s2, false );
            // if ( sdi2_1.getDuplicationsSum() != 0 ) {
            // return false;
            // }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }
}
