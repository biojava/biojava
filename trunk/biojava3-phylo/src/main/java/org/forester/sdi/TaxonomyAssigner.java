// $Id: TaxonomyAssigner.java,v 1.4 2009/10/17 02:42:51 cmzmasek Exp $
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

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

public class TaxonomyAssigner extends SDI {

    public TaxonomyAssigner( final Phylogeny gene_tree, final Phylogeny species_tree ) {
        super( gene_tree, species_tree );
        getSpeciesTree().preOrderReId( 0 );
        linkNodesOfG();
        geneTreePostOrderTraversal( getGeneTree().getRoot() );
    }

    void geneTreePostOrderTraversal( final PhylogenyNode g ) {
        if ( !g.isExternal() ) {
            for( final PhylogenyNodeIterator iter = g.iterateChildNodesForward(); iter.hasNext(); ) {
                geneTreePostOrderTraversal( iter.next() );
            }
            final PhylogenyNode[] linked_nodes = new PhylogenyNode[ g.getNumberOfDescendants() ];
            for( int i = 0; i < linked_nodes.length; ++i ) {
                linked_nodes[ i ] = g.getChildNode( i ).getLink();
            }
            final int[] min_max = GSDI.obtainMinMaxIdIndices( linked_nodes );
            int min_i = min_max[ 0 ];
            int max_i = min_max[ 1 ];
            while ( linked_nodes[ min_i ] != linked_nodes[ max_i ] ) {
                linked_nodes[ max_i ] = linked_nodes[ max_i ].getParent();
                final int[] min_max_ = GSDI.obtainMinMaxIdIndices( linked_nodes );
                min_i = min_max_[ 0 ];
                max_i = min_max_[ 1 ];
            }
            final PhylogenyNode s = linked_nodes[ max_i ];
            g.setLink( s );
            if ( s.getNodeData().isHasTaxonomy() ) {
                g.getNodeData().setTaxonomy( ( Taxonomy ) s.getNodeData().getTaxonomy().copy() );
            }
        }
    }

    public static void execute( final Phylogeny gene_tree, final Phylogeny species_tree ) {
        new TaxonomyAssigner( gene_tree, species_tree );
    }
}
