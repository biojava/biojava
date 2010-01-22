// $Id: SDI.java,v 1.14 2009/10/26 23:29:39 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
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
// WWW: www.phylosoft.org

package org.forester.sdi;

import java.util.HashMap;
import java.util.Map;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public abstract class SDI {

    final Phylogeny _gene_tree;
    final Phylogeny _species_tree;
    int             _duplications_sum; // Sum of duplications.
    int             _mapping_cost;    // Mapping cost "L".

    /**
     * Constructor which sets the gene tree and the species tree to be compared.
     * species_tree is the species tree to which the gene tree gene_tree will be
     * compared to.
     * Infers for each PhylogenyNode of gene_tree whether
     * it represents a speciation or duplication event by calculating and
     * interpreting the mapping function M. The most parsimonious sequence of
     * speciation and duplication events is assumed.
     * The mapping cost L can be
     * calculated with method "computeMappingCost()".
     * <p>
     * Conditions:
     * </p>
     * <ul>
     * <li>Both Trees must be rooted
     * <li>Both Trees must have species names in the species name fields of all
     * their external nodes
     * </ul>
     * 
     * @param gene_tree
     *            reference to a rooted binary gene Phylogeny to which assign
     *            duplication vs speciation, must have species names in the
     *            species name fields for all external nodes
     * @param species_tree
     *            reference to a rooted binary species Phylogeny which might get
     *            stripped in the process, must have species names in the
     *            species name fields for all external nodes
     */
    public SDI( final Phylogeny gene_tree, final Phylogeny species_tree ) {
        if ( species_tree.isEmpty() || gene_tree.isEmpty() ) {
            throw new IllegalArgumentException( "attempt to infer duplications using empty tree(s)" );
        }
        if ( !gene_tree.isRooted() ) {
            throw new IllegalArgumentException( "attempt to infer duplications on unrooted gene tree" );
        }
        if ( !species_tree.isRooted() ) {
            throw new IllegalArgumentException( "attempt to infer duplications on unrooted species tree" );
        }
        _gene_tree = gene_tree;
        _species_tree = species_tree;
        _duplications_sum = 0;
        _mapping_cost = -1;
    }

    // Helper method for "computeMappingCost()".
    private void computeMappingCostHelper( final PhylogenyNode g ) {
        if ( !g.isExternal() ) {
            computeMappingCostHelper( g.getChildNode1() );
            computeMappingCostHelper( g.getChildNode2() );
            if ( ( g.getLink() != g.getChildNode1().getLink() ) && ( g.getLink() != g.getChildNode2().getLink() ) ) {
                _mapping_cost += ( g.getChildNode1().getLink().getNodeId() + g.getChildNode2().getLink().getNodeId()
                        - ( 2 * g.getLink().getNodeId() ) - 2 );
            }
            else if ( ( g.getLink() != g.getChildNode1().getLink() ) && ( g.getLink() == g.getChildNode2().getLink() ) ) {
                _mapping_cost += ( g.getChildNode1().getLink().getNodeId() - g.getLink().getNodeId() + 1 );
            }
            else if ( ( g.getLink() == g.getChildNode1().getLink() ) && ( g.getLink() != g.getChildNode2().getLink() ) ) {
                _mapping_cost += ( g.getChildNode2().getLink().getNodeId() - g.getLink().getNodeId() + 1 );
            }
            else {
                _mapping_cost++;
            }
        }
    }

    /**
     * Computes the cost of mapping the gene tree gene_tree onto the species
     * tree species_tree. Before this method can be called, the mapping has to
     * be calculated with method "infer(boolean)".
     * <p>
     * Reference. Zhang, L. (1997) On a Mirkin-Muchnik-Smith Conjecture for
     * Comparing Molecular Phylogenies. Journal of Computational Biology 4
     * 177-187.
     * 
     * @return the mapping cost "L"
     */
    public int computeMappingCostL() {
        _species_tree.levelOrderReID( 0 );
        _mapping_cost = 0;
        computeMappingCostHelper( _gene_tree.getRoot() );
        return _mapping_cost;
    }

    private TaxonomyComparisonBase determineTaxonomyComparisonBase() {
        TaxonomyComparisonBase base = null;
        boolean all_have_id = true;
        boolean all_have_code = true;
        boolean all_have_sn = true;
        boolean all_have_cn = true;
        for( final PhylogenyNodeIterator iter = _species_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( n.getNodeData().isHasTaxonomy() ) {
                final Taxonomy tax = n.getNodeData().getTaxonomy();
                if ( ( tax.getIdentifier() == null ) || ForesterUtil.isEmpty( tax.getIdentifier().getValue() ) ) {
                    all_have_id = false;
                }
                if ( ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                    all_have_code = false;
                }
                if ( ForesterUtil.isEmpty( tax.getScientificName() ) ) {
                    all_have_sn = false;
                }
                if ( ForesterUtil.isEmpty( tax.getCommonName() ) ) {
                    all_have_cn = false;
                }
            }
            else {
                throw new IllegalArgumentException( "species tree node [" + n + "] has no taxonomic data" );
            }
        }
        for( final PhylogenyNodeIterator iter = _gene_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( n.getNodeData().isHasTaxonomy() ) {
                final Taxonomy tax = n.getNodeData().getTaxonomy();
                if ( ( tax.getIdentifier() == null ) || ForesterUtil.isEmpty( tax.getIdentifier().getValue() ) ) {
                    all_have_id = false;
                }
                if ( ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                    all_have_code = false;
                }
                if ( ForesterUtil.isEmpty( tax.getScientificName() ) ) {
                    all_have_sn = false;
                }
                if ( ForesterUtil.isEmpty( tax.getCommonName() ) ) {
                    all_have_cn = false;
                }
            }
            else {
                throw new IllegalArgumentException( "gene tree node [" + n + "] has no taxonomic data" );
            }
        }
        if ( all_have_id ) {
            base = TaxonomyComparisonBase.ID;
        }
        else if ( all_have_code ) {
            base = TaxonomyComparisonBase.CODE;
        }
        else if ( all_have_sn ) {
            base = TaxonomyComparisonBase.SCIENTIFIC_NAME;
        }
        else if ( all_have_cn ) {
            base = TaxonomyComparisonBase.COMMON_NAME;
        }
        else {
            throw new IllegalArgumentException( "gene tree and species tree have incomparable taxonomies" );
        }
        return base;
    }

    /**
     * Returns the number of duplications.
     * 
     * @return number of duplications
     */
    public int getDuplicationsSum() {
        return _duplications_sum;
    }

    /**
     * Returns the gene tree.
     * 
     * @return gene tree
     */
    public Phylogeny getGeneTree() {
        return _gene_tree;
    }

    /**
     * Returns the species tree.
     * 
     * @return species tree
     */
    public Phylogeny getSpeciesTree() {
        return _species_tree;
    }

    /**
     * Calculates the mapping function for the external nodes of the gene tree:
     * links (sets the field "link" of PhylogenyNode) each external
     * PhylogenyNode of gene_tree to the external PhylogenyNode of species_tree
     * which has the same species name.
     */
    void linkNodesOfG() {
        final Map<String, PhylogenyNode> speciestree_ext_nodes = new HashMap<String, PhylogenyNode>();
        final TaxonomyComparisonBase tax_comp_base = determineTaxonomyComparisonBase();
        // Put references to all external nodes of the species tree into a map.
        // Stringyfied taxonomy is the key, node is the value.
        for( final PhylogenyNodeIterator iter = _species_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode s = iter.next();
            final String tax_str = taxonomyToString( s, tax_comp_base );
            if ( speciestree_ext_nodes.containsKey( tax_str ) ) {
                throw new IllegalArgumentException( "taxonomy [" + s.getNodeData().getTaxonomy()
                        + "] is not unique in species phylogeny" );
            }
            speciestree_ext_nodes.put( tax_str, s );
        }
        // Retrieve the reference to the node with a matching stringyfied taxonomy.
        for( final PhylogenyNodeIterator iter = _gene_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode g = iter.next();
            final String tax_str = taxonomyToString( g, tax_comp_base );
            final PhylogenyNode s = speciestree_ext_nodes.get( tax_str );
            if ( s == null ) {
                throw new IllegalArgumentException( "taxonomy [" + g.getNodeData().getTaxonomy()
                        + "] not present in species tree" );
            }
            g.setLink( s );
        }
    }

    /**
     * Calculates the mapping function for the external nodes of the gene tree:
     * links (sets the field "link" of PhylogenyNode) each external by taxonomy
     * identifier
     * PhylogenyNode of gene_tree to the external PhylogenyNode of species_tree
     * which has the same species name.
     * Olivier CHABROL : olivier.chabrol@univ-provence.fr
     */
    void linkNodesOfGByTaxonomyIdentifier() {
        final HashMap<String, PhylogenyNode> speciestree_ext_nodes = new HashMap<String, PhylogenyNode>();
        if ( _species_tree.getFirstExternalNode().isRoot() ) {
            speciestree_ext_nodes.put( _species_tree.getFirstExternalNode().getNodeData().getTaxonomy().getIdentifier()
                    .getValue(), _species_tree.getFirstExternalNode() );
        }
        else {
            for( final PhylogenyNodeIterator iter = _species_tree.iteratorExternalForward(); iter.hasNext(); ) {
                final PhylogenyNode s = iter.next();
                speciestree_ext_nodes.put( s.getNodeData().getTaxonomy().getIdentifier().getValue(), s );
            }
        }
        for( final PhylogenyNodeIterator iter = _gene_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode g = iter.next();
            final PhylogenyNode s = speciestree_ext_nodes
                    .get( g.getNodeData().getTaxonomy().getIdentifier().getValue() );
            if ( s == null ) {
                String message = "species [" + g.getNodeData().getTaxonomy().getIdentifier().getValue();
                message += "] not present in species tree";
                throw new IllegalArgumentException( message );
            }
            g.setLink( s );
        }
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append( getClass() );
        sb.append( ForesterUtil.LINE_SEPARATOR );
        sb.append( "Duplications sum                   : " + getDuplicationsSum() );
        sb.append( ForesterUtil.LINE_SEPARATOR );
        sb.append( "mapping cost L                     : " + computeMappingCostL() );
        return sb.toString();
    }

    private static String taxonomyToString( final PhylogenyNode n, final TaxonomyComparisonBase base ) {
        final Taxonomy tax = n.getNodeData().getTaxonomy();
        switch ( base ) {
            case ID:
                return tax.getIdentifier().getValue();
            case CODE:
                return tax.getTaxonomyCode();
            case SCIENTIFIC_NAME:
                return tax.getScientificName();
            case COMMON_NAME:
                return tax.getCommonName();
            default:
                throw new IllegalArgumentException( "unknown comparison base for taxonomies: " + base );
        }
    }

    enum TaxonomyComparisonBase {
        ID, CODE, SCIENTIFIC_NAME, COMMON_NAME;
    }
}
