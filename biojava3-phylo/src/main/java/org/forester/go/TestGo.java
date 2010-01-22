// $Id: TestGo.java,v 1.28 2009/04/28 20:30:56 cmzmasek Exp $
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

package org.forester.go;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;

import org.forester.surfacing.DomainId;
import org.forester.util.ForesterUtil;

public class TestGo {

    private final static double ZERO_DIFF = 1.0E-9;

    public static boolean isEqual( final double a, final double b ) {
        return ( ( Math.abs( a - b ) ) < ZERO_DIFF );
    }

    public static boolean test( final File test_dir ) {
        System.out.print( "  GO ID: " );
        if ( !testGoId() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Namespace: " );
        if ( !testNamespace() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Basic GO term: " );
        if ( !testBasicGoTerm() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  OBO parser: " );
        if ( !testOBOparser( test_dir ) ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Pfam to GO mapping: " );
        if ( !testPfamToGoMapping() ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Pfam to GO parser: " );
        if ( !testPfamToGoParser( test_dir ) ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Super terms: " );
        if ( !testSuperTermGetting( test_dir ) ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        System.out.print( "  Super term counting: " );
        if ( !testSuperTermCounting( test_dir ) ) {
            System.out.println( "failed." );
            return false;
        }
        System.out.println( "OK." );
        return true;
    }

    private static boolean testBasicGoTerm() {
        try {
            final GoTerm gt1 = new BasicGoTerm( "GO:0047579",
                                                "4-hydroxymandelate oxidase activity",
                                                "molecular_function",
                                                false );
            final GoTerm gt2 = new BasicGoTerm( "GO:0047579",
                                                "4-hydroxymandelate oxidase activity",
                                                "molecular_function",
                                                false );
            final GoTerm gt3 = new BasicGoTerm( "GO:0047579", "?", "molecular_function", true );
            final GoTerm gt4 = new BasicGoTerm( "GO:0047579",
                                                "4-hydroxymandelate oxidase activity",
                                                "biological_process",
                                                false );
            final GoTerm gt5 = new BasicGoTerm( "GO:0047578",
                                                "4-hydroxymandelate oxidase activity",
                                                "molecular_function",
                                                false );
            if ( !gt1.equals( gt2 ) ) {
                return false;
            }
            if ( !gt1.equals( gt3 ) ) {
                return false;
            }
            if ( gt1.equals( gt4 ) ) {
                return false;
            }
            if ( gt1.hashCode() != gt4.hashCode() ) {
                return false;
            }
            if ( gt1.equals( gt5 ) ) {
                return false;
            }
            final GoTerm gt6 = ( GoTerm ) gt5.copy();
            if ( !gt6.equals( gt5 ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testGoId() {
        try {
            final GoId id1 = new GoId( "GO:0042617" );
            final GoId id2 = new GoId( "GO:0042630" );
            final GoId id3 = new GoId( "GO:0042630" );
            if ( id1.equals( id2 ) ) {
                return false;
            }
            if ( !id2.equals( id3 ) ) {
                return false;
            }
            if ( !id1.toString().equals( "GO:0042617" ) ) {
                return false;
            }
            if ( id2.hashCode() != id3.hashCode() ) {
                return false;
            }
            if ( id1.hashCode() == id2.hashCode() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testNamespace() {
        try {
            final GoNameSpace b = new GoNameSpace( "Biological_process" );
            final GoNameSpace c = new GoNameSpace( "Cellular_Component" );
            final GoNameSpace m = new GoNameSpace( "molecular_function" );
            final GoNameSpace m2 = new GoNameSpace( GoNameSpace.GoNamespaceType.MOLECULAR_FUNCTION );
            if ( b.equals( c ) ) {
                return false;
            }
            if ( !m.equals( m2 ) ) {
                return false;
            }
            if ( !b.toString().equals( "biological_process" ) ) {
                return false;
            }
            if ( !c.toString().equals( "cellular_component" ) ) {
                return false;
            }
            if ( !m.toString().equals( "molecular_function" ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testOBOparser( final File test_dir ) {
        try {
            final OBOparser parser = new OBOparser( new File( test_dir + ForesterUtil.getFileSeparator() + "obo_test" ),
                                                    OBOparser.ReturnType.BASIC_GO_TERM );
            final List<GoTerm> go_terms = parser.parse();
            if ( parser.getGoTermCount() != 26 ) {
                return false;
            }
            final GoTerm g0 = go_terms.get( 0 );
            final GoTerm g1 = go_terms.get( 1 );
            final GoTerm g3 = go_terms.get( 2 );
            final GoTerm g2 = go_terms.get( 25 );
            if ( !g0.getComment().equals( "" ) ) {
                return false;
            }
            if ( !g0
                    .getDefinition()
                    .equals( "\"The distribution of mitochondria, including the mitochondrial genome, into daughter cells after mitosis or meiosis, mediated by interactions between mitochondria and the cytoskeleton.\" [GOC:mcc, PMID:10873824, PMID:11389764]" ) ) {
                return false;
            }
            if ( !g0.getGoId().getId().equals( "GO:0000001" ) ) {
                return false;
            }
            if ( g0.getGoNameSpace().equals( GoNameSpace.GoNamespaceType.BIOLOGICAL_PROCESS ) ) {
                return false;
            }
            if ( g0.getGoNameSpace().getType() != GoNameSpace.GoNamespaceType.BIOLOGICAL_PROCESS ) {
                return false;
            }
            if ( g0.getGoRelationships().size() != 0 ) {
                return false;
            }
            if ( g0.getGoXRefs().size() != 0 ) {
                return false;
            }
            if ( !g0.getName().equals( "mitochondrion inheritance" ) ) {
                return false;
            }
            if ( g0.getSuperGoIds().size() != 2 ) {
                return false;
            }
            if ( !g0.isObsolete() ) {
                return false;
            }
            if ( !g1.getComment().equals( "comment" ) ) {
                return false;
            }
            if ( !g1
                    .getDefinition()
                    .equals( "\"The maintenance of the structure and integrity of the mitochondrial genome.\" [GOC:ai]" ) ) {
                return false;
            }
            if ( !g1.getGoId().getId().equals( "GO:0000002" ) ) {
                return false;
            }
            if ( g1.getGoNameSpace().equals( GoNameSpace.GoNamespaceType.BIOLOGICAL_PROCESS ) ) {
                return false;
            }
            if ( g1.getGoNameSpace().getType() != GoNameSpace.GoNamespaceType.BIOLOGICAL_PROCESS ) {
                return false;
            }
            if ( g1.getGoRelationships().size() != 1 ) {
                return false;
            }
            if ( g1.getGoXRefs().size() != 5 ) {
                return false;
            }
            if ( !g1.getName().equals( "mitochondrial genome maintenance" ) ) {
                return false;
            }
            if ( g1.getSuperGoIds().size() != 1 ) {
                return false;
            }
            if ( g1.isObsolete() ) {
                return false;
            }
            if ( !g1.getGoXRefs().get( 0 ).equals( new BasicGoXRef( "EC:2.4.1.-" ) ) ) {
                return false;
            }
            if ( !g1.getGoXRefs().get( 0 ).getXRef().equals( "2.4.1.-" ) ) {
                return false;
            }
            if ( g1.getGoXRefs().get( 0 ).getType() != GoXRef.Type.EC ) {
                return false;
            }
            if ( g1.getGoXRefs().get( 0 ).equals( new BasicGoXRef( "EC:2.4.1.1" ) ) ) {
                return false;
            }
            if ( g1.getGoXRefs().get( 0 ).equals( new BasicGoXRef( "Reactome:2.4.1.-" ) ) ) {
                return false;
            }
            if ( !g1.getGoXRefs().get( 1 ).equals( new BasicGoXRef( "Reactome:7672" ) ) ) {
                return false;
            }
            if ( !g1.getGoXRefs().get( 2 ).equals( new BasicGoXRef( "MetaCyc:SIROHEME-FERROCHELAT-RXN" ) ) ) {
                return false;
            }
            if ( !g1.getGoXRefs().get( 3 ).equals( new BasicGoXRef( "RESID:AA02376" ) ) ) {
                return false;
            }
            if ( !g1.getGoXRefs().get( 4 ).equals( new BasicGoXRef( "UM-BBD_enzymeID:e0271" ) ) ) {
                return false;
            }
            if ( !g1.getGoRelationships().get( 0 ).equals( new BasicGoRelationship( "part_of GO:0007052" ) ) ) {
                return false;
            }
            if ( !g1.getGoRelationships().get( 0 ).getGoId().equals( new GoId( "GO:0007052" ) ) ) {
                return false;
            }
            if ( !g1.getGoRelationships().get( 0 ).getGoId().getId().equals( "GO:0007052" ) ) {
                return false;
            }
            if ( g1.getGoRelationships().get( 0 ).getType() != GoRelationship.Type.PART_OF ) {
                return false;
            }
            if ( g1.getGoRelationships().get( 0 ).equals( new BasicGoRelationship( "part_of GO:1007052" ) ) ) {
                return false;
            }
            if ( !g1.getSuperGoIds().get( 0 ).equals( new GoId( "GO:0007005" ) ) ) {
                return false;
            }
            if ( g1.getSuperGoIds().get( 0 ).equals( new GoId( "GO:1007005" ) ) ) {
                return false;
            }
            if ( !g2.getGoId().getId().equals( "GO:0000030" ) ) {
                return false;
            }
            if ( !g2.getGoId().equals( new GoId( "GO:0000030" ) ) ) {
                return false;
            }
            if ( g2.getGoId().getId().equals( "GO:0000031" ) ) {
                return false;
            }
            if ( g2.getGoId().equals( new GoId( "GO:0000031" ) ) ) {
                return false;
            }
            if ( g3.getGoSubsets().size() != 3 ) {
                return false;
            }
            if ( !g3.getGoSubsets().contains( new BasicGoSubset( "goslim_generic" ) ) ) {
                return false;
            }
            if ( !g3.getGoSubsets().contains( new BasicGoSubset( "goslim_plant" ) ) ) {
                return false;
            }
            if ( !g3.getGoSubsets().contains( new BasicGoSubset( "gosubset_prok" ) ) ) {
                return false;
            }
            if ( g3.getGoSubsets().contains( new BasicGoSubset( "goslim_candida" ) ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testPfamToGoMapping() {
        try {
            final PfamToGoMapping pg0 = new PfamToGoMapping( new DomainId( "A" ), new GoId( "GO:0000001" ) );
            final PfamToGoMapping pg1 = new PfamToGoMapping( new DomainId( "A" ), new GoId( "GO:0000001" ) );
            final PfamToGoMapping pg2 = new PfamToGoMapping( new DomainId( "B" ), new GoId( "GO:0000001" ) );
            final PfamToGoMapping pg3 = new PfamToGoMapping( new DomainId( "A" ), new GoId( "GO:0000002" ) );
            final PfamToGoMapping pg4 = new PfamToGoMapping( new DomainId( "B" ), new GoId( "GO:0000002" ) );
            if ( !pg0.equals( pg0 ) ) {
                return false;
            }
            if ( !pg0.equals( pg1 ) ) {
                return false;
            }
            if ( pg0.equals( pg2 ) ) {
                return false;
            }
            if ( pg0.equals( pg3 ) ) {
                return false;
            }
            if ( pg0.equals( pg4 ) ) {
                return false;
            }
            if ( pg0.compareTo( pg3 ) != 0 ) {
                return false;
            }
            if ( pg0.compareTo( pg2 ) >= 0 ) {
                return false;
            }
            if ( pg2.compareTo( pg0 ) <= 0 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testPfamToGoParser( final File test_dir ) {
        try {
            final PfamToGoParser parser = new PfamToGoParser( new File( test_dir + ForesterUtil.getFileSeparator()
                    + "pfam_to_go_test" ) );
            final List<PfamToGoMapping> mappings = parser.parse();
            if ( parser.getMappingCount() != 426 ) {
                return false;
            }
            if ( mappings.size() != 426 ) {
                return false;
            }
            final PfamToGoMapping m0 = mappings.get( 0 );
            final PfamToGoMapping m1 = mappings.get( 1 );
            final PfamToGoMapping m2 = mappings.get( 2 );
            final PfamToGoMapping m3 = mappings.get( 3 );
            final PfamToGoMapping m4 = mappings.get( 4 );
            final PfamToGoMapping m5 = mappings.get( 5 );
            final PfamToGoMapping m424 = mappings.get( 424 );
            final PfamToGoMapping m425 = mappings.get( 425 );
            if ( !m0.getKey().equals( new DomainId( "7tm_1" ) ) ) {
                return false;
            }
            if ( !m0.getValue().equals( new GoId( "GO:0001584" ) ) ) {
                return false;
            }
            if ( m0.getKey().equals( new DomainId( "7tm_x" ) ) ) {
                return false;
            }
            if ( m0.getValue().equals( new GoId( "GO:0001585" ) ) ) {
                return false;
            }
            if ( !m1.getKey().equals( new DomainId( "7tm_1" ) ) ) {
                return false;
            }
            if ( !m1.getValue().equals( new GoId( "GO:0007186" ) ) ) {
                return false;
            }
            if ( !m2.getKey().equals( new DomainId( "7tm_1" ) ) ) {
                return false;
            }
            if ( !m2.getValue().equals( new GoId( "GO:0016021" ) ) ) {
                return false;
            }
            if ( !m3.getKey().equals( new DomainId( "7tm_2" ) ) ) {
                return false;
            }
            if ( !m3.getValue().equals( new GoId( "GO:0004930" ) ) ) {
                return false;
            }
            if ( !m4.getKey().equals( new DomainId( "7tm_2" ) ) ) {
                return false;
            }
            if ( !m4.getValue().equals( new GoId( "GO:0016020" ) ) ) {
                return false;
            }
            if ( !m5.getKey().equals( new DomainId( "7tm_3" ) ) ) {
                return false;
            }
            if ( !m5.getValue().equals( new GoId( "GO:0008067" ) ) ) {
                return false;
            }
            if ( !m424.getKey().equals( new DomainId( "OMPdecase" ) ) ) {
                return false;
            }
            if ( !m424.getValue().equals( new GoId( "GO:0006207" ) ) ) {
                return false;
            }
            if ( !m425.getKey().equals( new DomainId( "Bac_DNA_binding" ) ) ) {
                return false;
            }
            if ( !m425.getValue().equals( new GoId( "GO:0003677" ) ) ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSuperTermCounting( final File test_dir ) {
        try {
            final OBOparser parser = new OBOparser( new File( test_dir + ForesterUtil.getFileSeparator()
                    + "gene_ontology_edit.obo" ), OBOparser.ReturnType.BASIC_GO_TERM );
            final List<GoTerm> all_go_terms = parser.parse();
            if ( parser.getGoTermCount() != 27748 ) {
                return false;
            }
            final Map<GoId, GoTerm> goid_to_term_map = GoUtils.createGoIdToGoTermMap( all_go_terms );
            final List<GoTerm> categories = new ArrayList<GoTerm>();
            final List<GoTerm> experiment_set = new ArrayList<GoTerm>();
            experiment_set.add( new BasicGoTerm( new GoId( "GO:0005690" ), "snRNP U4atac", GoNameSpace
                    .createUnassigned(), false ) );
            experiment_set.add( new BasicGoTerm( new GoId( "GO:0009698" ),
                                                 "phenylpropanoid metabolic process",
                                                 GoNameSpace.createUnassigned(),
                                                 false ) );
            experiment_set.add( new BasicGoTerm( new GoId( "GO:0008150" ), "biological_process", GoNameSpace
                    .createUnassigned(), false ) );
            experiment_set.add( new BasicGoTerm( new GoId( "GO:0006915" ),
                                                 "apoptosis",
                                                 GoNameSpace.createUnassigned(),
                                                 false ) );
            experiment_set.add( new BasicGoTerm( new GoId( "GO:0001783" ), "B cell apoptosis", GoNameSpace
                    .createUnassigned(), false ) );
            experiment_set.add( new BasicGoTerm( new GoId( "GO:0010657" ), "muscle cell apoptosis", GoNameSpace
                    .createUnassigned(), false ) );
            experiment_set.add( new BasicGoTerm( new GoId( "GO:0010657" ), "muscle cell apoptosis", GoNameSpace
                    .createUnassigned(), false ) );
            experiment_set.add( new BasicGoTerm( new GoId( "GO:0010658" ),
                                                 "striated muscle cell apoptosis",
                                                 GoNameSpace.createUnassigned(),
                                                 false ) );
            experiment_set.add( new BasicGoTerm( new GoId( "GO:0043065" ),
                                                 "positive regulation of apoptosis",
                                                 GoNameSpace.createUnassigned(),
                                                 false ) );
            categories
                    .add( new BasicGoTerm( new GoId( "GO:0016265" ), "death", GoNameSpace.createUnassigned(), false ) );
            categories.add( new BasicGoTerm( new GoId( "GO:0006915" ),
                                             "apoptosis",
                                             GoNameSpace.createUnassigned(),
                                             false ) );
            categories.add( new BasicGoTerm( new GoId( "GO:0008150" ), "biological_process", GoNameSpace
                    .createUnassigned(), false ) );
            categories.add( new BasicGoTerm( new GoId( "GO:0010657" ), "muscle cell apoptosis", GoNameSpace
                    .createUnassigned(), false ) );
            categories.add( new BasicGoTerm( new GoId( "GO:0010658" ), "striated muscle cell apoptosis", GoNameSpace
                    .createUnassigned(), false ) );
            categories.add( new BasicGoTerm( new GoId( "GO:0046242" ), "o-xylene biosynthetic process", GoNameSpace
                    .createUnassigned(), false ) );
            categories.add( new BasicGoTerm( new GoId( "GO:0016326" ), "kinesin motor activity", GoNameSpace
                    .createUnassigned(), false ) );
            categories.add( new BasicGoTerm( new GoId( "GO:0005575" ), "cellular_component", GoNameSpace
                    .createUnassigned(), false ) );
            categories.add( new BasicGoTerm( new GoId( "GO:0032502" ), "developmental process", GoNameSpace
                    .createUnassigned(), false ) );
            categories.add( new BasicGoTerm( new GoId( "GO:0051094" ),
                                             "positive regulation of developmental process",
                                             GoNameSpace.createUnassigned(),
                                             false ) );
            categories.add( new BasicGoTerm( new GoId( "GO:0048522" ),
                                             "positive regulation of cellular process",
                                             GoNameSpace.createUnassigned(),
                                             false ) );
            final Map<GoId, Integer> counts = GoUtils.countCategories( categories, experiment_set, goid_to_term_map );
            // death
            if ( counts.get( new GoId( "GO:0016265" ) ) != 5 ) {
                return false;
            }
            // apoptosis
            if ( counts.get( new GoId( "GO:0006915" ) ) != 5 ) {
                return false;
            }
            // biological_process
            if ( counts.get( new GoId( "GO:0008150" ) ) != 8 ) {
                return false;
            }
            // muscle cell apoptosis
            if ( counts.get( new GoId( "GO:0010657" ) ) != 3 ) {
                return false;
            }
            // striated muscle cell apoptosis
            if ( counts.get( new GoId( "GO:0010658" ) ) != 1 ) {
                return false;
            }
            // o-xylene biosynthetic process
            if ( counts.get( new GoId( "GO:0046242" ) ) != 0 ) {
                return false;
            }
            // kinesin motor activity
            if ( counts.get( new GoId( "GO:0016326" ) ) != 0 ) {
                return false;
            }
            // cellular_component
            if ( counts.get( new GoId( "GO:0005575" ) ) != 1 ) {
                return false;
            }
            // developmental process
            if ( counts.get( new GoId( "GO:0032502" ) ) != 5 ) {
                return false;
            }
            // positive regulation of developmental process
            if ( counts.get( new GoId( "GO:0051094" ) ) != 1 ) {
                return false;
            }
            // positive regulation of cellular process
            if ( counts.get( new GoId( "GO:0048522" ) ) != 1 ) {
                return false;
            }
            final List<GoId> categories_id = new ArrayList<GoId>();
            final List<GoId> experiment_set_id = new ArrayList<GoId>();
            experiment_set_id.add( new GoId( "GO:0005690" ) );
            experiment_set_id.add( new GoId( "GO:0009698" ) );
            experiment_set_id.add( new GoId( "GO:0008150" ) );
            experiment_set_id.add( new GoId( "GO:0006915" ) );
            experiment_set_id.add( new GoId( "GO:0001783" ) );
            experiment_set_id.add( new GoId( "GO:0010657" ) );
            experiment_set_id.add( new GoId( "GO:0010657" ) );
            experiment_set_id.add( new GoId( "GO:0010658" ) );
            categories_id.add( new GoId( "GO:0016265" ) );
            categories_id.add( new GoId( "GO:0006915" ) );
            categories_id.add( new GoId( "GO:0008150" ) );
            categories_id.add( new GoId( "GO:0010657" ) );
            categories_id.add( new GoId( "GO:0010658" ) );
            categories_id.add( new GoId( "GO:0046242" ) );
            categories_id.add( new GoId( "GO:0016326" ) );
            categories_id.add( new GoId( "GO:0005575" ) );
            final Map<GoId, Integer> counts_id = GoUtils.countCategoriesId( categories_id,
                                                                            experiment_set_id,
                                                                            goid_to_term_map );
            // death
            if ( counts_id.get( new GoId( "GO:0016265" ) ) != 5 ) {
                return false;
            }
            // apoptosis
            if ( counts_id.get( new GoId( "GO:0006915" ) ) != 5 ) {
                return false;
            }
            // biological_process
            if ( counts_id.get( new GoId( "GO:0008150" ) ) != 7 ) {
                return false;
            }
            // muscle cell apoptosis
            if ( counts_id.get( new GoId( "GO:0010657" ) ) != 3 ) {
                return false;
            }
            // striated muscle cell apoptosis
            if ( counts_id.get( new GoId( "GO:0010658" ) ) != 1 ) {
                return false;
            }
            // o-xylene biosynthetic process
            if ( counts_id.get( new GoId( "GO:0046242" ) ) != 0 ) {
                return false;
            }
            // kinesin motor activity
            if ( counts_id.get( new GoId( "GO:0016326" ) ) != 0 ) {
                return false;
            }
            // cellular_componen
            if ( counts_id.get( new GoId( "GO:0005575" ) ) != 1 ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
        return true;
    }

    private static boolean testSuperTermGetting( final File test_dir ) {
        try {
            final OBOparser parser = new OBOparser( new File( test_dir + ForesterUtil.getFileSeparator()
                    + "gene_ontology_edit.obo" ), OBOparser.ReturnType.BASIC_GO_TERM );
            final List<GoTerm> go_terms = parser.parse();
            if ( parser.getGoTermCount() != 27748 ) {
                return false;
            }
            final Map<GoId, GoTerm> goid_to_term_map = GoUtils.createGoIdToGoTermMap( go_terms );
            final SortedSet<GoTerm> b_cell_selection = GoUtils.getAllSuperGoTerms( new GoId( "GO:0002339" ),
                                                                                   goid_to_term_map );
            if ( b_cell_selection.size() != 2 ) {
                return false;
            }
            if ( !b_cell_selection.contains( new BasicGoTerm( new GoId( "GO:0002376" ),
                                                              "immune system process",
                                                              GoNameSpace.createBiologicalProcess(),
                                                              false ) ) ) {
                return false;
            }
            if ( !b_cell_selection.contains( new BasicGoTerm( new GoId( "GO:0008150" ),
                                                              "biological process",
                                                              GoNameSpace.createBiologicalProcess(),
                                                              false ) ) ) {
                return false;
            }
            final SortedSet<GoTerm> b_cell_differentation = GoUtils.getAllSuperGoTerms( new GoId( "GO:0030183" ),
                                                                                        goid_to_term_map );
            if ( b_cell_differentation.size() != 12 ) {
                return false;
            }
            final SortedSet<GoTerm> biological_process = GoUtils.getAllSuperGoTerms( new GoId( "GO:0008150" ),
                                                                                     goid_to_term_map );
            if ( biological_process.size() != 0 ) {
                return false;
            }
            final SortedSet<GoTerm> protein_aa_phosphorylation = GoUtils.getAllSuperGoTerms( new GoId( "GO:0006468" ),
                                                                                             goid_to_term_map );
            if ( protein_aa_phosphorylation.size() != 16 ) {
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
