// $Id: RIOn.java,v 1.6 2009/10/26 23:29:39 cmzmasek Exp $
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

import java.util.ArrayList;
import java.util.List;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.GeneralTable;

public class RIOn {

    private final static boolean  ROOT_BY_MINIMIZING_MAPPING_COST = false;
    private final static boolean  ROOT_BY_MINIMIZING_SUM_OF_DUPS  = true;
    private final static boolean  ROOT_BY_MINIMIZING_TREE_HEIGHT  = true;
    GeneralTable<String, Integer> _orthologs                      = null;
    GeneralTable<String, Integer> _paralogs                       = null;
    GeneralTable<String, Integer> _super_orthologs                = null;
    GeneralTable<String, Integer> _ultra_paralogs                 = null;

    private void doInferOrthologs( final Phylogeny gene_tree, final Phylogeny species_tree ) {
        final SDIR sdiunrooted = new SDIR();
        final Phylogeny assigned_tree = sdiunrooted.infer( gene_tree,
                                                           species_tree,
                                                           ROOT_BY_MINIMIZING_MAPPING_COST,
                                                           ROOT_BY_MINIMIZING_SUM_OF_DUPS,
                                                           ROOT_BY_MINIMIZING_TREE_HEIGHT,
                                                           true,
                                                           1 )[ 0 ];
        final List<PhylogenyNode> external_nodes = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator iterator = assigned_tree.iteratorExternalForward(); iterator.hasNext(); ) {
            external_nodes.add( iterator.next() );
        }
        final PhylogenyMethods methods = PhylogenyMethods.getInstance();
        for( int i = 0; i < external_nodes.size(); ++i ) {
            for( int j = 0; j < external_nodes.size(); ++j ) {
                if ( i != j ) {
                    final PhylogenyNode node_i = external_nodes.get( i );
                    final PhylogenyNode node_j = external_nodes.get( j );
                    final PhylogenyNode lca = methods.getLCA( node_i, node_j );
                    final Event event = lca.getNodeData().getEvent();
                    final String node_i_name = node_i.getNodeData().getSequence().getName();
                    final String node_j_name = node_j.getNodeData().getSequence().getName();
                    if ( event.isDuplication() ) {
                        increaseCounter( getOrthologs(), node_i_name, node_j_name );
                    }
                    else {
                        increaseCounter( getParalogs(), node_i_name, node_j_name );
                    }
                }
            }
        }
    }

    public GeneralTable<String, Integer> getOrthologs() {
        return _orthologs;
    }

    public GeneralTable<String, Integer> getParalogs() {
        return _paralogs;
    }

    public GeneralTable<String, Integer> getSuperOrthologs() {
        return _super_orthologs;
    }

    public GeneralTable<String, Integer> getUltraParalogs() {
        return _ultra_paralogs;
    }

    private void increaseCounter( final GeneralTable<String, Integer> table,
                                  final String node_i_name,
                                  final String node_j_name ) {
        final Integer value = table.getValue( node_i_name, node_j_name );
        if ( value == null ) {
            table.setValue( node_i_name, node_j_name, 1 );
        }
        else {
            table.setValue( node_i_name, node_j_name, value.intValue() + 1 );
        }
    }

    private void init() {
        _orthologs = new GeneralTable<String, Integer>();
        _paralogs = new GeneralTable<String, Integer>();
        _super_orthologs = new GeneralTable<String, Integer>();
        _ultra_paralogs = new GeneralTable<String, Integer>();
    }

    private void setOrthologs( final GeneralTable<String, Integer> orthologs ) {
        _orthologs = orthologs;
    }

    private void setParalogs( final GeneralTable<String, Integer> paralogs ) {
        _paralogs = paralogs;
    }

    private void setSuperOrthologs( final GeneralTable<String, Integer> super_orthologs ) {
        _super_orthologs = super_orthologs;
    }

    private void setUltraParalogs( final GeneralTable<String, Integer> ultra_paralogs ) {
        _ultra_paralogs = ultra_paralogs;
    }
}
