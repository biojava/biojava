// $Id: Shin.java,v 1.3 2009/10/26 23:29:39 cmzmasek Exp $
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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class Shin {

    public Shin() {
    }

    private void analyze( final Phylogeny gene_tree,
                          final String gene_tree_file_name,
                          final Phylogeny[] species_trees,
                          final File out_dir ) throws IOException {
        final boolean minimize_cost = true;
        final boolean minimize_sum_of_dup = true;
        final boolean minimize_height = true;
        final int trees_to_return = 1;
        System.out.println( gene_tree_file_name + ": " + gene_tree.getName() );
        final Set<Taxonomy> species_tree_species = getAllExternalSpecies( species_trees[ 0 ] );
        final PhylogenyWriter w = new PhylogenyWriter();
        for( final Phylogeny species_tree : species_trees ) {
            PhylogenyMethods.deleteExternalNodesPositiveSelection( species_tree_species, gene_tree );
            if ( gene_tree.isEmpty() ) {
                System.out.println( " >> empty: " + gene_tree_file_name + ": " + gene_tree.getName() );
                continue;
            }
            final File outfile = new File( out_dir + ForesterUtil.FILE_SEPARATOR + gene_tree_file_name );
            if ( outfile.exists() ) {
                System.out
                        .println( " >> already exists, skipping: " + gene_tree_file_name + ": " + gene_tree.getName() );
            }
            final SDIR sdir = new SDIR();
            final Phylogeny[] analyzed_gene_trees = sdir.infer( gene_tree,
                                                                species_tree,
                                                                minimize_cost,
                                                                minimize_sum_of_dup,
                                                                minimize_height,
                                                                true,
                                                                trees_to_return );
            final int duplications = sdir.getMinimalDuplications();
            final int mapping_cost = sdir.getMinimalMappingCost();
            final List<Phylogeny> phys = new ArrayList<Phylogeny>();
            for( final Phylogeny phy : analyzed_gene_trees ) {
                phys.add( phy );
            }
            w.toPhyloXML( outfile, phys, 0, ForesterUtil.LINE_SEPARATOR );
        }
    }

    private void checkSpeciesTreesForEqualNumberOfExtNodes( final Phylogeny[] species_trees ) {
        int ext_nodes = -1;
        for( final Phylogeny phylogeny : species_trees ) {
            if ( ext_nodes < 0 ) {
                ext_nodes = phylogeny.getNumberOfExternalNodes();
            }
            else if ( ext_nodes != phylogeny.getNumberOfExternalNodes() ) {
                throw new IllegalArgumentException( "species trees must have all the same number of external nodes" );
            }
        }
    }

    public void method1( final List<File> gene_tree_files, final Phylogeny[] species_trees, final File out_dir )
            throws IOException {
        checkSpeciesTreesForEqualNumberOfExtNodes( species_trees );
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        for( final File gene_tree_file : gene_tree_files ) {
            if ( ForesterUtil.isReadableFile( gene_tree_file ) != null ) {
                throw new IOException( "[" + gene_tree_file + "] is not readable" );
            }
            Phylogeny[] gene_trees = null;
            gene_trees = factory.create( gene_tree_file, new PhyloXmlParser() );
            if ( gene_trees.length != 1 ) {
                throw new IOException( "[" + gene_tree_file + "] contains " + gene_trees.length
                        + " gene trees, expecting precisely one" );
            }
            analyze( gene_trees[ 0 ], gene_tree_file.getName(), species_trees, out_dir );
        }
    }

    private static Set<Taxonomy> getAllExternalSpecies( final Phylogeny phy ) {
        final Set<Taxonomy> specs = new HashSet<Taxonomy>();
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy() ) {
                specs.add( n.getNodeData().getTaxonomy() );
            }
            else {
                throw new IllegalArgumentException( "node " + n.getNodeId() + " has no taxonomic data" );
            }
        }
        return specs;
    }
}
