// $Id: SDIse.java,v 1.12 2009/10/26 23:29:39 cmzmasek Exp $
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
// WWW: www.phylosoft.org/forester

package org.forester.sdi;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Event;

/*
 * Implements our algorithm for speciation - duplication inference (SDI). <p>
 * Reference: </p> <ul> <li>Zmasek, C.M. and Eddy, S.R. (2001) "A simple
 * algorithm to infer gene duplication and speciation events on a gene tree".
 * Bioinformatics, in press. </ul> <p> The initialization is accomplished by:
 * </p> <ul> <li>method "linkExtNodesOfG()" of class SDI: setting the links for
 * the external nodes of the gene tree <li>"preorderReID(int)" from class
 * Phylogeny: numbering of nodes of the species tree in preorder <li>the
 * optional stripping of the species tree is accomplished by method
 * "stripTree(Phylogeny,Phylogeny)" of class Phylogeny </ul> <p> The recursion
 * part is accomplished by this class' method
 * "geneTreePostOrderTraversal(PhylogenyNode)". <p> Requires JDK 1.2 or greater.
 * 
 * @see SDI#linkNodesOfG()
 * 
 * @see Phylogeny#preorderReID(int)
 * 
 * @see
 * PhylogenyMethods#taxonomyBasedDeletionOfExternalNodes(Phylogeny,Phylogeny)
 * 
 * @see #geneTreePostOrderTraversal(PhylogenyNode)
 * 
 * @author Christian M. Zmasek
 * 
 * @version 1.102 -- last modified: 10/02/01
 */
public class SDIse extends SDI {

    /**
     * Constructor which sets the gene tree and the species tree to be compared.
     * species_tree is the species tree to which the gene tree gene_tree will be
     * compared to - with method "infer(boolean)". Both Trees must be completely
     * binary and rooted. The actual inference is accomplished with method
     * "infer(boolean)". The mapping cost L can then be calculated with method
     * "computeMappingCost()".
     * <p>
     * (Last modified: 01/11/01)
     * 
     * @see #infer(boolean)
     * @see SDI#computeMappingCostL()
     * @param gene_tree
     *            reference to a rooted binary gene Phylogeny to which assign
     *            duplication vs speciation, must have species names in the
     *            species name fields for all external nodes
     * @param species_tree
     *            reference to a rooted binary species Phylogeny which might get
     *            stripped in the process, must have species names in the
     *            species name fields for all external nodes
     */
    public SDIse( final Phylogeny gene_tree, final Phylogeny species_tree ) {
        super( gene_tree, species_tree );
        _duplications_sum = 0;
        getSpeciesTree().preOrderReId( 0 );
        linkNodesOfG();
        geneTreePostOrderTraversal( getGeneTree().getRoot() );
    }

    // Helper method for updateM( boolean, PhylogenyNode, PhylogenyNode )
    // Calculates M for PhylogenyNode n, given that M for the two children
    // of n has been calculated.
    // (Last modified: 10/02/01)
    private void calculateMforNode( final PhylogenyNode n ) {
        if ( !n.isExternal() ) {
            final boolean was_duplication = n.isDuplication();
            PhylogenyNode a = n.getChildNode1().getLink(), b = n.getChildNode2().getLink();
            while ( a != b ) {
                if ( a.getNodeId() > b.getNodeId() ) {
                    a = a.getParent();
                }
                else {
                    b = b.getParent();
                }
            }
            n.setLink( a );
            Event event = null;
            if ( ( a == n.getChildNode1().getLink() ) || ( a == n.getChildNode2().getLink() ) ) {
                event = Event.createSingleDuplicationEvent();
                if ( !was_duplication ) {
                    ++_duplications_sum;
                }
            }
            else {
                event = Event.createSingleSpeciationEvent();
                if ( was_duplication ) {
                    --_duplications_sum;
                }
            }
            n.getNodeData().setEvent( event );
        }
    } // calculateMforNode( PhylogenyNode )

    /**
     * Traverses the subtree of PhylogenyNode g in postorder, calculating the
     * mapping function M, and determines which nodes represent speciation
     * events and which ones duplication events.
     * <p>
     * Preconditions: Mapping M for external nodes must have been calculated and
     * the species tree must be labelled in preorder.
     * <p>
     * (Last modified: 01/11/01)
     * 
     * @param g
     *            starting node of a gene tree - normally the root
     */
    void geneTreePostOrderTraversal( final PhylogenyNode g ) {
        PhylogenyNode a, b;
        if ( !g.isExternal() ) {
            geneTreePostOrderTraversal( g.getChildNode( 0 ) );
            geneTreePostOrderTraversal( g.getChildNode( 1 ) );
            a = g.getChildNode( 0 ).getLink();
            b = g.getChildNode( 1 ).getLink();
            while ( a != b ) {
                if ( a.getNodeId() > b.getNodeId() ) {
                    a = a.getParent();
                }
                else {
                    b = b.getParent();
                }
            }
            g.setLink( a );
            // Determines whether dup. or spec.
            Event event = null;
            if ( ( a == g.getChildNode( 0 ).getLink() ) || ( a == g.getChildNode( 1 ).getLink() ) ) {
                event = Event.createSingleDuplicationEvent();
                ++_duplications_sum;
            }
            else {
                event = Event.createSingleSpeciationEvent();
            }
            g.getNodeData().setEvent( event );
        }
    } // geneTreePostOrderTraversal( PhylogenyNode )

    /**
     * Updates the mapping function M after the root of the gene tree has been
     * moved by one branch. It calculates M for the root of the gene tree and
     * one of its two children.
     * <p>
     * To be used ONLY by method "SDIunrooted.fastInfer(Phylogeny,Phylogeny)".
     * <p>
     * (Last modfied: 10/02/01)
     * 
     * @param prev_root_was_dup
     *            true if the previous root was a duplication, false otherwise
     * @param prev_root_c1
     *            child 1 of the previous root
     * @param prev_root_c2
     *            child 2 of the previous root
     * @return number of duplications which have been assigned in gene tree
     */
    int updateM( final boolean prev_root_was_dup, final PhylogenyNode prev_root_c1, final PhylogenyNode prev_root_c2 ) {
        final PhylogenyNode root = getGeneTree().getRoot();
        if ( ( root.getChildNode1() == prev_root_c1 ) || ( root.getChildNode2() == prev_root_c1 ) ) {
            calculateMforNode( prev_root_c1 );
        }
        else {
            calculateMforNode( prev_root_c2 );
        }
        Event event = null;
        if ( prev_root_was_dup ) {
            event = Event.createSingleDuplicationEvent();
        }
        else {
            event = Event.createSingleSpeciationEvent();
        }
        root.getNodeData().setEvent( event );
        calculateMforNode( root );
        return getDuplicationsSum();
    } // updateM( boolean, PhylogenyNode, PhylogenyNode )
} // End of class SDIse.
