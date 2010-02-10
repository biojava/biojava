// $Id: support_transfer.java,v 1.14 2009/11/20 22:22:09 cmzmasek Exp $
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

package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public final class support_transfer {

    /**
     * Transfers branch length values from one Phylogeny to another. It is
     * mainly a "main method" for method "copyBranchLengthValuesFrom( Phylogeny )"
     * of org.forester.phylogeny.Phylogeny, to be used in other (Perl) programs.
     * 
     * @param args[0]
     *            Filename (String) for Phylogeny which has correct branch
     *            length values
     * @param args[1]
     *            String Filename (String) for Phylogeny to which the branch
     *            lengths of the first Phylogeny are to be copied, both Trees
     *            must only differ in their branch length values, i.e. topology
     *            and sequence names, etc. must be the same
     * @param args[2]
     *            String Filename (String) for outputfile
     * @param args[3]
     *            String [number of tree with correct bl to use in case treefile contains more than one, default 0]            
     
     */
    public static void main( final String args[] ) {
        Phylogeny phylogeny_w_bl = null; // Has correct branch lengths
        Phylogeny phylogeny_w_support_vals = null; // Has bootsrap in the b.l.
        // field (will be
        // transferred
        // to the bootstrap field by the Phylogeny constructor) or
        // has regular boostraps (NHX, :B=...).
        File infile_bl = null;
        File infile_support_vals = null;
        File outfile = null;
        int index_of_tree_w_bl = 0;
        if ( ( args.length != 3 ) && ( args.length != 4 ) ) {
            System.err.println( "SupportTransfer: Wrong number" + " of arguments. Usage: \"java transfersBranchLenghts"
                    + " <treefile with correct b.l.> <treefile with bootstraps>" + "<outputfile> "
                    + "[number of tree with correct bl to use in case treefile contains more than one, default 0]\"" );
            System.exit( -1 );
        }
        if ( args.length == 4 ) {
            index_of_tree_w_bl = ( new Integer( args[ 3 ] ) ).intValue();
        }
        try {
            infile_bl = new File( args[ 0 ] );
            infile_support_vals = new File( args[ 1 ] );
            outfile = new File( args[ 2 ] );
            if ( outfile.exists() ) {
                System.out.println( "transfersBranchLenghts: " + outfile.getAbsolutePath() + " does already exist." );
                System.exit( -1 );
            }
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp_bl = ForesterUtil.createParserDependingOnFileType( infile_bl, true );
            final PhylogenyParser pp_s = ForesterUtil.createParserDependingOnFileType( infile_support_vals, true );
            if ( pp_bl instanceof NHXParser ) {
                ( ( NHXParser ) pp_bl ).setTaxonomyExtraction( TAXONOMY_EXTRACTION.YES );
            }
            phylogeny_w_bl = factory.create( infile_bl, pp_bl )[ index_of_tree_w_bl ];
            phylogeny_w_support_vals = factory.create( infile_support_vals, pp_s )[ 0 ];
        }
        catch ( final IOException e ) {
            System.out.println( "SupportTransfer: Could not read tree(s): " + e );
            System.exit( -1 );
        }
        try {
            final double max_bs = PhylogenyMethods.getMaximumConfidenceValue( phylogeny_w_support_vals );
            PhylogenyMethods.normalizeBootstrapValues( phylogeny_w_support_vals, max_bs, 100 );
            support_transfer.transferSupportValues( phylogeny_w_support_vals, phylogeny_w_bl );
        }
        catch ( final IllegalArgumentException e ) {
            System.out.println( e.getMessage() );
            System.exit( -1 );
        }
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( outfile, phylogeny_w_bl, 0 );
        }
        catch ( final IOException e ) {
            System.out.println( "Failure to write phylogeny \'" + outfile + "\" [" + e.getMessage() + "]" );
            System.exit( -1 );
        }
    }

    /**
     * Moves the values in the branch length field to the bootstrap field, for
     * each PhylogenyNode of this Phylogeny. Converts a Phylogeny originating
     * from a phylip treefile after bootstrapping and which therefore has its
     * bootstrap values where the branch lenghts would be.
     */
    public final static void moveBranchLengthsToBootstrap( final Phylogeny p ) {
        for( final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isInternal() && ( node.getDistanceToParent() > 0 ) ) {
                PhylogenyMethods.setBootstrapConfidence( node, node.getDistanceToParent() );
            }
            else {
                PhylogenyMethods.setBootstrapConfidence( node, Confidence.CONFIDENCE_DEFAULT_VALUE );
            }
            node.setDistanceToParent( PhylogenyNode.DISTANCE_DEFAULT );
        }
    } // moveBranchLengthsToBootstrap()

    /**
     * Modifies Phylogeny to with the support values from Phylogeny from.
     * Important (but obvious): The topology of both trees needs to be the same.
     * The method is not robust, and might produce wrong results if the internal
     * topology differs or if the external node names are not unique.
     * 
     * @param from
     *            the Phylogeny to copy the support values from
     * @param to
     *            the Phylogeny to copy the support values to
     */
    public final static void transferSupportValues( final Phylogeny from, final Phylogeny to ) {
        to: for( final PhylogenyNodeIterator it_to = to.iteratorPostorder(); it_to.hasNext(); ) {
            final PhylogenyNode node_to = it_to.next();
            if ( !node_to.isExternal() ) {
                final List<String> ext_children_to = node_to.getAllExternalDescendantsNames();
                for( final PhylogenyNodeIterator it_from = from.iteratorPostorder(); it_from.hasNext(); ) {
                    final PhylogenyNode node_from = it_from.next();
                    final List<String> ext_children_from = node_from.getAllExternalDescendantsNames();
                    if ( ( ext_children_from.size() == ext_children_to.size() )
                            && ext_children_from.containsAll( ext_children_to ) ) {
                        PhylogenyMethods.setBootstrapConfidence( node_to, PhylogenyMethods
                                .getConfidenceValue( node_from ) );
                        continue to;
                    }
                }
                final String message = "Attempt to transfer support values from nonidentical topologies";
                throw new IllegalArgumentException( message );
            }
        }
    }
}