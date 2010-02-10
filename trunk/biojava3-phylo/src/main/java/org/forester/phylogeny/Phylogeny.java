// $Id: Phylogeny.java,v 1.92 2009/10/26 23:29:40 cmzmasek Exp $
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

package org.forester.phylogeny;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.data.BranchData;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.iterators.ExternalForwardIterator;
import org.forester.phylogeny.iterators.LevelOrderTreeIterator;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PostorderTreeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.ForesterUtil;

public class Phylogeny {

    public final static boolean             ALLOW_MULTIPLE_PARENTS_DEFAULT = false;
    private PhylogenyNode                   _root;
    private boolean                         _rooted;
    private boolean                         _allow_multiple_parents;
    private String                          _name;
    private String                          _type;
    private String                          _description;
    private String                          _distance_unit;
    private Confidence                      _confidence;
    private Identifier                      _identifier;
    private boolean                         _rerootable;
    private HashMap<Integer, PhylogenyNode> _idhash;
    private Set<PhylogenyNode>              _external_nodes_set;

    /**
     * Default Phylogeny constructor. Constructs an empty Phylogeny.
     */
    public Phylogeny() {
        init();
    }

    /**
     * Adds this Phylogeny to the list of child nodes of PhylogenyNode parent
     * and sets the parent of this to parent.
     * 
     * @param n
     *            the PhylogenyNode to add
     */
    public void addAsChild( final PhylogenyNode parent ) {
        if ( isEmpty() ) {
            throw new IllegalArgumentException( "Attempt to add an empty tree." );
        }
        if ( !isRooted() ) {
            throw new IllegalArgumentException( "Attempt to add an unrooted tree." );
        }
        parent.addAsChild( getRoot() );
        externalNodesHaveChanged();
    }

    public void addAsSibling( final PhylogenyNode sibling ) {
        if ( isEmpty() ) {
            throw new IllegalArgumentException( "Attempt to add an empty tree." );
        }
        if ( !isRooted() ) {
            throw new IllegalArgumentException( "Attempt to add an unrooted tree." );
        }
        final int sibling_index = sibling.getChildNodeIndex();
        final PhylogenyNode new_node = new PhylogenyNode();
        final PhylogenyNode sibling_parent = sibling.getParent();
        new_node.setChild1( sibling );
        new_node.setChild2( getRoot() );
        new_node.setParent( sibling_parent );
        sibling.setParent( new_node );
        sibling_parent.setChildNode( sibling_index, new_node );
        final double new_dist = sibling.getDistanceToParent() == PhylogenyNode.DISTANCE_DEFAULT ? PhylogenyNode.DISTANCE_DEFAULT
                : sibling.getDistanceToParent() / 2;
        new_node.setDistanceToParent( new_dist );
        sibling.setDistanceToParent( new_dist );
        externalNodesHaveChanged();
    }

    /**
     * This calculates the height of the subtree emanating at n for rooted,
     * tree-shaped phylogenies
     * 
     * @param n
     *            the root-node of a subtree
     * @return the height of the subtree emanating at n
     */
    public double calculateSubtreeHeight( final PhylogenyNode n ) {
        if ( n.isExternal() || n.isCollapse() ) {
            return ForesterUtil.isLargerOrEqualToZero( n.getDistanceToParent() );
        }
        else {
            double max = -Double.MAX_VALUE;
            for( int i = 0; i < n.getNumberOfDescendants(); ++i ) {
                final double l = calculateSubtreeHeight( n.getChildNode( i ) );
                if ( l > max ) {
                    max = l;
                }
            }
            return max + ForesterUtil.isLargerOrEqualToZero( n.getDistanceToParent() );
        }
    }

    /**
     * Returns a deep copy of this Phylogeny.
     * <p>
     * (The resulting Phylogeny has its references in the external nodes
     * corrected, if they are lacking/obsolete in this.)
     */
    public Phylogeny copy() {
        final Phylogeny tree = new Phylogeny();
        if ( isEmpty() ) {
            tree.init();
            return tree;
        }
        tree._rooted = _rooted;
        tree._name = new String( _name );
        tree._description = new String( _description );
        tree._type = new String( _type );
        tree._rerootable = _rerootable;
        tree._distance_unit = new String( _distance_unit );
        if ( _confidence != null ) {
            tree._confidence = ( Confidence ) _confidence.copy();
        }
        if ( _identifier != null ) {
            tree._identifier = ( Identifier ) _identifier.copy();
        }
        tree.setAllowMultipleParents( isAllowMultipleParents() );
        tree._root = PhylogenyMethods.copySubTree( _root );
        return tree;
    }

    /**
     * Need the delete and/or rehash _idhash (not done automatically
     * to allow client multiple deletions in linear time).
     * Need to call 'recalculateNumberOfExternalDescendants(boolean)' after this 
     * if tree is to be displayed.
     * 
     * @param remove_us the parent node of the subtree to be deleted
     */
    public void deleteSubtree( final PhylogenyNode remove_us, final boolean collapse_resulting_node_with_one_desc ) {
        if ( isEmpty() ) {
            return;
        }
        if ( remove_us.isRoot() ) {
            init();
            return;
        }
        if ( !collapse_resulting_node_with_one_desc ) {
            remove_us.getParent().removeChildNode( remove_us );
        }
        else {
            final PhylogenyNode removed_node = remove_us;
            final PhylogenyNode p = remove_us.getParent();
            if ( p.isRoot() ) {
                if ( p.getNumberOfDescendants() == 2 ) {
                    if ( removed_node.isFirstChildNode() ) {
                        setRoot( getRoot().getChildNode( 1 ) );
                        getRoot().setParent( null );
                    }
                    else {
                        setRoot( getRoot().getChildNode( 0 ) );
                        getRoot().setParent( null );
                    }
                }
                else {
                    p.removeChildNode( removed_node.getChildNodeIndex() );
                }
            }
            else {
                final PhylogenyNode pp = removed_node.getParent().getParent();
                if ( p.getNumberOfDescendants() == 2 ) {
                    final int pi = p.getChildNodeIndex();
                    if ( removed_node.isFirstChildNode() ) {
                        p.getChildNode( 1 ).setDistanceToParent( PhylogenyMethods.addPhylogenyDistances( p
                                .getDistanceToParent(), p.getChildNode( 1 ).getDistanceToParent() ) );
                        pp.setChildNode( pi, p.getChildNode( 1 ) );
                    }
                    else {
                        p.getChildNode( 0 ).setDistanceToParent( PhylogenyMethods.addPhylogenyDistances( p
                                .getDistanceToParent(), p.getChildNode( 0 ).getDistanceToParent() ) );
                        pp.setChildNode( pi, p.getChildNode( 0 ) );
                    }
                }
                else {
                    p.removeChildNode( removed_node.getChildNodeIndex() );
                }
            }
        }
        remove_us.setParent( null );
        setIdHash( null );
        externalNodesHaveChanged();
    }

    public void externalNodesHaveChanged() {
        _external_nodes_set = null;
    }

    public String[] getAllExternalNodeNames() {
        int i = 0;
        if ( isEmpty() ) {
            return null;
        }
        final String[] names = new String[ getNumberOfExternalNodes() ];
        for( final PhylogenyNodeIterator iter = iteratorExternalForward(); iter.hasNext(); ) {
            names[ i++ ] = new String( iter.next().getNodeName() );
        }
        return names;
    }

    public Confidence getConfidence() {
        return _confidence;
    }

    public String getDescription() {
        return _description;
    }

    public String getDistanceUnit() {
        return _distance_unit;
    }

    /**
     * 
     * Warning. The order of the returned nodes is random
     * -- and hence cannot be relied on.
     * 
     * @return Unordered set of PhylogenyNode
     */
    public Set<PhylogenyNode> getExternalNodes() {
        if ( _external_nodes_set == null ) {
            _external_nodes_set = new HashSet<PhylogenyNode>();
            for( final PhylogenyNodeIterator it = iteratorPostorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( n.isExternal() ) {
                    _external_nodes_set.add( n );
                }
            }
        }
        return _external_nodes_set;
    }

    /**
     * Returns the number of duplications of this Phylogeny (int). A return
     * value of -1 indicates that the number of duplications is unknown.
     */
    // public int getNumberOfDuplications() {
    // return _number_of_duplications;
    // } // getNumberOfDuplications()
    /**
     * Sets the number of duplications of this Phylogeny (int). A value of -1
     * indicates that the number of duplications is unknown.
     * 
     * @param clean_nh
     *            set to true for clean NH format
     */
    // public void setNumberOfDuplications( int i ) {
    // if ( i < 0 ) {
    // _number_of_duplications = -1;
    // }
    // else {
    // _number_of_duplications = i;
    // }
    // } // setNumberOfDuplications( int )
    /**
     * Returns the first external PhylogenyNode.
     */
    public PhylogenyNode getFirstExternalNode() {
        if ( isEmpty() ) {
            throw new IllegalStateException( "Attempt to obtain first external node of empty Phylogeney." );
        }
        PhylogenyNode node = getRoot();
        while ( node.isInternal() ) {
            node = node.getFirstChildNode();
        }
        return node;
    }

    /**
     * This calculates the height for rooted, tree-shaped phylogenies. The
     * height is the longest distance from the root to an external node. Please
     * note. Child nodes of collapsed nodes are ignored -- which is useful for
     * display purposes but might be misleading for other applications.
     * 
     * @return the height for rooted, tree-shaped phylogenies
     */
    public double getHeight() {
        if ( isEmpty() ) {
            return 0.0;
        }
        return calculateSubtreeHeight( getRoot() );
    }

    public Identifier getIdentifier() {
        return _identifier;
    }

    // ---------------------------------------------------------
    // Modification of Phylogeny topology and Phylogeny appearance
    // ---------------------------------------------------------
    private HashMap<Integer, PhylogenyNode> getIdHash() {
        return _idhash;
    }

    /**
     * Returns the name of this Phylogeny.
     */
    public String getName() {
        return _name;
    }

    /**
     * Finds the PhylogenyNode of this Phylogeny which has a matching ID number.
     * Takes O(n) time. After method hashIDs() has been called it runs in
     * constant time.
     * 
     * @param id
     *            ID number (int) of the PhylogenyNode to find
     * @return PhylogenyNode with matching ID, null if not found
     */
    public PhylogenyNode getNode( final int id ) throws NoSuchElementException {
        if ( isEmpty() ) {
            throw new NoSuchElementException( "attempt to get node in an empty phylogeny" );
        }
        if ( _idhash != null ) {
            return _idhash.get( id );
        }
        else {
            for( final PhylogenyNodeIterator iter = iteratorPreorder(); iter.hasNext(); ) {
                final PhylogenyNode node = iter.next();
                if ( node.getNodeId() == id ) {
                    return node;
                }
            }
        }
        return null;
    }

    /**
     * Returns a PhylogenyNode of this Phylogeny which has a matching name.
     * Throws an Exception if seqname is not present in this or not unique.
     * 
     * @param name
     *            name (String) of PhylogenyNode to find
     * @return PhylogenyNode with matchin name
     */
    public PhylogenyNode getNode( final String name ) {
        if ( isEmpty() ) {
            return null;
        }
        final List<PhylogenyNode> nodes = getNodes( name );
        if ( ( nodes == null ) || ( nodes.size() < 1 ) ) {
            throw new IllegalArgumentException( "node named [" + name + "] not found" );
        }
        if ( nodes.size() > 1 ) {
            throw new IllegalArgumentException( "node named [" + name + "] not unique" );
        }
        return nodes.get( 0 );
    }

    /**
     * Return Node by TaxonomyId Olivier CHABROL :
     * olivier.chabrol@univ-provence.fr
     * 
     * @param taxonomyID
     *            search taxonomy identifier
     * @param nodes
     *            sublist node to search
     * @return List node with the same taxonomy identifier
     */
    private List<PhylogenyNode> getNodeByTaxonomyID( final String taxonomyID, final List<PhylogenyNode> nodes ) {
        final List<PhylogenyNode> retour = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNode node : nodes ) {
            if ( taxonomyID.equals( PhylogenyMethods.getTaxonomyIdentifier( node ) ) ) {
                retour.add( node );
            }
        }
        return retour;
    }

    /**
     * Returns a List with references to all Nodes of this Phylogeny which have
     * a matching name.
     * 
     * @param name
     *            name (String) of Nodes to find
     * @return Vector of references to Nodes of this Phylogeny with matching
     *         names
     * @see #getNodesWithMatchingSpecies(String)
     */
    public List<PhylogenyNode> getNodes( final String name ) {
        if ( isEmpty() ) {
            return null;
        }
        final List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator iter = iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( n.getNodeName().equals( name ) ) {
                nodes.add( n );
            }
        }
        return nodes;
    }

    public List<PhylogenyNode> getNodesViaSequenceName( final String seq_name ) {
        if ( isEmpty() ) {
            return null;
        }
        final List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator iter = iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( n.getNodeData().isHasSequence() && n.getNodeData().getSequence().getName().equals( seq_name ) ) {
                nodes.add( n );
            }
        }
        return nodes;
    }

    public List<PhylogenyNode> getNodesViaTaxonomyCode( final String taxonomy_code ) {
        if ( isEmpty() ) {
            return null;
        }
        final List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator iter = iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( n.getNodeData().isHasTaxonomy()
                    && n.getNodeData().getTaxonomy().getTaxonomyCode().equals( taxonomy_code ) ) {
                nodes.add( n );
            }
        }
        return nodes;
    }

    /**
     * Returns a Vector with references to all Nodes of this Phylogeny which
     * have a matching species name.
     * 
     * @param specname
     *            species name (String) of Nodes to find
     * @return Vector of references to Nodes of this Phylogeny with matching
     *         species names.
     * @see #getNodes(String)
     */
    public List<PhylogenyNode> getNodesWithMatchingSpecies( final String specname ) {
        if ( isEmpty() ) {
            return null;
        }
        final List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator iter = iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( PhylogenyMethods.getSpecies( n ).equals( specname ) ) {
                nodes.add( n );
            }
        }
        return nodes;
    }

    public PhylogenyNode getNodeViaSequenceName( final String seq_name ) {
        if ( isEmpty() ) {
            return null;
        }
        final List<PhylogenyNode> nodes = getNodesViaSequenceName( seq_name );
        if ( ( nodes == null ) || ( nodes.size() < 1 ) ) {
            throw new IllegalArgumentException( "node with sequence named [" + seq_name + "] not found" );
        }
        if ( nodes.size() > 1 ) {
            throw new IllegalArgumentException( "node with sequence named [" + seq_name + "] not unique" );
        }
        return nodes.get( 0 );
    }

    public PhylogenyNode getNodeViaTaxonomyCode( final String taxonomy_code ) {
        if ( isEmpty() ) {
            return null;
        }
        final List<PhylogenyNode> nodes = getNodesViaTaxonomyCode( taxonomy_code );
        if ( ( nodes == null ) || ( nodes.size() < 1 ) ) {
            throw new IllegalArgumentException( "node with taxonomy code \"" + taxonomy_code + "\" not found" );
        }
        if ( nodes.size() > 1 ) {
            throw new IllegalArgumentException( "node with taxonomy code \"" + taxonomy_code + "\" not unique" );
        }
        return nodes.get( 0 );
    }

    public int getNumberOfBranches() {
        if ( isEmpty() ) {
            return 0;
        }
        int c = 0;
        for( final PhylogenyNodeIterator iter = iteratorPreorder(); iter.hasNext(); iter.next() ) {
            ++c;
        }
        if ( !isRooted() ) {
            --c;
        }
        return c;
    }

    /**
     * Returns the sum of external Nodes of this Phylogeny (int).
     */
    public int getNumberOfExternalNodes() {
        if ( isEmpty() ) {
            return 0;
        }
        return getExternalNodes().size();
    }

    /**
     * Returns all paralogs of the external PhylogenyNode n of this Phylogeny.
     * paralog are returned as List of node references.
     * <p>
     * PRECONDITION: This tree must be binary and rooted, and speciation -
     * duplication need to be assigned for each of its internal Nodes.
     * <p>
     * Returns null if this Phylogeny is empty or if n is internal.
     * <p>
     * (Last modified: 11/22/00) Olivier CHABROL :
     * olivier.chabrol@univ-provence.fr
     * 
     * @param n
     *            external PhylogenyNode whose orthologs are to be returned
     * @return Vector of references to all orthologous Nodes of PhylogenyNode n
     *         of this Phylogeny, null if this Phylogeny is empty or if n is
     *         internal
     */
    public List<PhylogenyNode> getParalogousNodes( final PhylogenyNode n, final String[] taxonomyCodeRange ) {
        PhylogenyNode node = n;
        PhylogenyNode prev = null;
        final List<PhylogenyNode> v = new ArrayList<PhylogenyNode>();
        final Map<PhylogenyNode, List<String>> map = new HashMap<PhylogenyNode, List<String>>();
        getTaxonomyMap( getRoot(), map );
        if ( !node.isExternal() || isEmpty() ) {
            return null;
        }
        final String searchNodeSpeciesId = PhylogenyMethods.getTaxonomyIdentifier( n );
        if ( !node.isExternal() || isEmpty() ) {
            return null;
        }
        List<String> taxIdList = null;
        final List<String> taxonomyCodeRangeList = Arrays.asList( taxonomyCodeRange );
        while ( !node.isRoot() ) {
            prev = node;
            node = node.getParent();
            taxIdList = map.get( node );
            if ( node.isDuplication() && isContains( taxIdList, taxonomyCodeRangeList ) ) {
                if ( node.getChildNode1() == prev ) {
                    v.addAll( getNodeByTaxonomyID( searchNodeSpeciesId, node.getChildNode2()
                            .getAllExternalDescendants() ) );
                }
                else {
                    v.addAll( getNodeByTaxonomyID( searchNodeSpeciesId, node.getChildNode1()
                            .getAllExternalDescendants() ) );
                }
            }
        }
        return v;
    }

    /**
     * Returns the root PhylogenyNode of this Phylogeny.
     */
    public PhylogenyNode getRoot() {
        return _root;
    }

    /**
     * List all species contains in all leaf under a node Olivier CHABROL :
     * olivier.chabrol@univ-provence.fr
     * 
     * @param node
     *            PhylogenyNode whose sub node species are returned
     * @return species contains in all leaf under the param node
     */
    private List<String> getSubNodeTaxonomy( final PhylogenyNode node ) {
        final List<String> taxonomyList = new ArrayList<String>();
        final List<PhylogenyNode> childs = node.getAllExternalDescendants();
        String speciesId = null;
        for( final PhylogenyNode phylogenyNode : childs ) {
            // taxId = new Long(phylogenyNode.getTaxonomyID());
            speciesId = PhylogenyMethods.getTaxonomyIdentifier( phylogenyNode );
            if ( !taxonomyList.contains( speciesId ) ) {
                taxonomyList.add( speciesId );
            }
        }
        return taxonomyList;
    }

    /**
     * Create a map [<PhylogenyNode, List<String>], the list contains the
     * species contains in all leaf under phylogeny node Olivier CHABROL :
     * olivier.chabrol@univ-provence.fr
     * 
     * @param node
     *            the tree root node
     * @param map
     *            map to fill
     */
    private void getTaxonomyMap( final PhylogenyNode node, final Map<PhylogenyNode, List<String>> map ) {
        // node is leaf
        if ( node.isExternal() ) {
            return;
        }
        map.put( node, getSubNodeTaxonomy( node ) );
        getTaxonomyMap( node.getChildNode1(), map );
        getTaxonomyMap( node.getChildNode2(), map );
    }

    public String getType() {
        return _type;
    }

    /**
     * Hashes the ID number of each PhylogenyNode of this Phylogeny to its
     * corresponding PhylogenyNode, in order to make method getNode( id ) run in
     * constant time. Important: The user is responsible for calling this method
     * (again) after this Phylogeny has been changed/created/renumbered.
     */
    public void hashIDs() {
        if ( isEmpty() ) {
            return;
        }
        setIdHash( new HashMap<Integer, PhylogenyNode>() );
        for( final PhylogenyNodeIterator iter = iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            getIdHash().put( node.getNodeId(), node );
        }
    }

    /**
     * Deletes this Phylogeny.
     */
    public void init() {
        _root = null;
        _rooted = false;
        _name = "";
        _description = "";
        _type = "";
        _distance_unit = "";
        _idhash = null;
        _confidence = null;
        _identifier = null;
        _rerootable = true;
        setAllowMultipleParents( Phylogeny.ALLOW_MULTIPLE_PARENTS_DEFAULT );
    }

    private boolean isAllowMultipleParents() {
        return _allow_multiple_parents;
    }

    /**
     * Returns whether this is a completely binary tree (i.e. all internal nodes
     * are bifurcations).
     * 
     */
    public boolean isCompletelyBinary() {
        if ( isEmpty() ) {
            return false;
        }
        for( final PhylogenyNodeIterator iter = iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isInternal() && ( node.getNumberOfDescendants() != 2 ) ) {
                return false;
            }
        }
        return true;
    }

    /**
     * Util method to check if all element of a list is contains in the
     * rangeList. Olivier CHABROL : olivier.chabrol@univ-provence.fr
     * 
     * @param list
     *            list to be check
     * @param rangeList
     *            the range list to compare
     * @return <code>true</code> if all param list element are contains in param
     *         rangeList, <code>false</code> otherwise.
     */
    private boolean isContains( final List<String> list, final List<String> rangeList ) {
        if ( list.size() > rangeList.size() ) {
            return false;
        }
        String l = null;
        for( final Iterator<String> iterator = list.iterator(); iterator.hasNext(); ) {
            l = iterator.next();
            if ( !rangeList.contains( l ) ) {
                return false;
            }
        }
        return true;
    }

    /**
     * Checks whether a Phylogeny object is deleted (or empty).
     * 
     * @return true if the tree is deleted (or empty), false otherwise
     */
    public boolean isEmpty() {
        return ( getRoot() == null );
    }

    public boolean isRerootable() {
        return _rerootable;
    }

    /**
     * Returns true is this Phylogeny is rooted.
     */
    public boolean isRooted() {
        return _rooted;
    } // isRooted()

    public boolean isTree() {
        return true;
    }

    public PhylogenyNodeIterator iteratorExternalForward() {
        return new ExternalForwardIterator( this );
    }

    public PhylogenyNodeIterator iteratorLevelOrder() {
        return new LevelOrderTreeIterator( this );
    }

    public PhylogenyNodeIterator iteratorPostorder() {
        return new PostorderTreeIterator( this );
    }

    public PhylogenyNodeIterator iteratorPreorder() {
        return new PreorderTreeIterator( this );
    }

    /**
     * Resets the ID numbers of the nodes of this Phylogeny in level order,
     * starting with start_label (for the root). <br>
     * WARNING. After this method has been called, node IDs are no longer
     * unique. <br>
     * WARNING. Afterwards static method getNodeCount of PhylogenyNode is might
     * not be correct anymore (depending on what start_label was used). Returns
     * the maximum used ID number which is used for the child nodes furthest
     * away from the root.
     * 
     * @param start_label
     *            the starting value for the ID numbers
     * @return the maximum used ID number
     */
    public int levelOrderReID( final int start_label ) {
        if ( isEmpty() ) {
            return start_label;
        }
        _idhash = null;
        int max = Integer.MIN_VALUE;
        for( final PhylogenyNodeIterator it = iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( node.isRoot() ) {
                node.setNodeId( start_label );
            }
            else {
                node.setNodeId( node.getParent().getNodeId() + 1 );
                if ( node.getNodeId() > max ) {
                    max = node.getNodeId();
                }
            }
        }
        return max;
    }

    /**
     * Arranges the order of childern for each node of this Phylogeny in such a
     * way that either the branch with more children is on top (right) or on
     * bottom (left), dependent on the value of boolean order.
     * 
     * @param order
     *            decides in which direction to order
     */
    public void orderAppearance( final boolean order ) throws IllegalStateException {
        if ( !isTree() ) {
            throw new IllegalStateException( "Attempt to order appearance on phylogeny which is not tree-like." );
        }
        if ( isEmpty() ) {
            return;
        }
        orderAppearanceHelper( getRoot(), order );
    }

    // Helper method for "orderAppearance(boolean)".
    // Traverses this Phylogeny recusively.
    private void orderAppearanceHelper( final PhylogenyNode n, final boolean order ) {
        if ( n.isExternal() ) {
            return;
        }
        else {
            PhylogenyNode temp = null;
            // FIXME
            if ( ( n.getNumberOfDescendants() == 2 )
                    && ( n.getChildNode1().getNumberOfExternalNodes() != n.getChildNode2().getNumberOfExternalNodes() )
                    && ( ( n.getChildNode1().getNumberOfExternalNodes() < n.getChildNode2().getNumberOfExternalNodes() ) == order ) ) {
                temp = n.getChildNode1();
                n.setChild1( n.getChildNode2() );
                n.setChild2( temp );
            }
            for( int i = 0; i < n.getNumberOfDescendants(); ++i ) {
                orderAppearanceHelper( n.getChildNode( i ), order );
            }
        }
    }

    public void preOrderReId() {
        if ( isEmpty() ) {
            return;
        }
        setIdHash( null );
        int i = PhylogenyNode.getNodeCount() + 1;
        for( final PhylogenyNodeIterator it = iteratorPreorder(); it.hasNext(); ) {
            it.next().setNodeId( i++ );
        }
        PhylogenyNode.setNodeCount( i );
    }

    /**
     * Resets the ID numbers of the Nodes of this Phylogeny in preorder,
     * starting with i.<br>
     * WARNING. Afterwards static method getNodeCount of PhylogenyNode is might
     * not be correct anymore (depending on what start_label was used).
     * 
     * @param start_label
     *            the starting value (int)
     * @return start_label plus the total number of Nodes of this Phylogeny
     *         (int)
     */
    public int preOrderReId( int start_label ) {
        if ( isEmpty() ) {
            return start_label;
        }
        setIdHash( null );
        for( final PhylogenyNodeIterator it = iteratorPreorder(); it.hasNext(); ) {
            it.next().setNodeId( start_label++ );
        }
        return start_label;
    }

    /**
     * Prints descriptions of all external Nodes of this Phylogeny to
     * System.out.
     */
    public void printExtNodes() {
        if ( isEmpty() ) {
            return;
        }
        for( final PhylogenyNodeIterator iter = iteratorExternalForward(); iter.hasNext(); ) {
            System.out.println( iter.next() + "\n" );
        }
    }

    /**
     * (Re)counts the number of children for each PhylogenyNode of this
     * Phylogeny. As an example, this method needs to be called after a
     * Phylogeny has been reRooted and it is to be displayed.
     * 
     * @param consider_collapsed_nodes
     *            set to true to take into account collapsed nodes (collapsed
     *            nodes have 1 child).
     */
    public void recalculateNumberOfExternalDescendants( final boolean consider_collapsed_nodes ) {
        if ( isEmpty() ) {
            return;
        }
        for( final PhylogenyNodeIterator iter = iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() || ( consider_collapsed_nodes && node.isCollapse() ) ) {
                node.setSumExtNodes( 1 );
            }
            else {
                int sum = 0;
                for( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
                    sum += node.getChildNode( i ).getNumberOfExternalNodes();
                }
                node.setSumExtNodes( sum );
            }
        }
    }

    /**
     * Places the root of this Phylogeny on the parent branch of the
     * PhylogenyNode with a corresponding ID. The new root is always placed on
     * the middle of the branch. If the resulting reRooted Phylogeny is to be
     * used any further, in most cases the following methods have to be called
     * on the resulting Phylogeny:
     * <p>
     * <li>recalculateNumberOfExternalDescendants(boolean)
     * <li>recalculateAndReset()
     * 
     * @param id
     *            ID (int) of PhylogenyNode of this Phylogeny
     */
    public void reRoot( final int id ) {
        reRoot( getNode( id ) );
    }

    /**
     * Places the root of this Phylogeny on Branch b. The new root is always
     * placed on the middle of the branch b.
     * 
     */
    public void reRoot( final PhylogenyBranch b ) {
        final PhylogenyNode n1 = b.getFirstNode();
        final PhylogenyNode n2 = b.getSecondNode();
        if ( n1.isExternal() ) {
            reRoot( n1 );
        }
        else if ( n2.isExternal() ) {
            reRoot( n2 );
        }
        else if ( ( n2 == n1.getChildNode1() ) || ( n2 == n1.getChildNode2() ) ) {
            reRoot( n2 );
        }
        else if ( ( n1 == n2.getChildNode1() ) || ( n1 == n2.getChildNode2() ) ) {
            reRoot( n1 );
        }
        else if ( ( n1.getParent() != null ) && n1.getParent().isRoot()
                && ( ( n1.getParent().getChildNode1() == n2 ) || ( n1.getParent().getChildNode2() == n2 ) ) ) {
            reRoot( n1 );
        }
        else {
            throw new IllegalArgumentException( "reRoot( Branch b ): b is not a branch." );
        }
    }

    /**
     * Places the root of this Phylogeny on the parent branch PhylogenyNode n.
     * The new root is always placed on the middle of the branch.
     * <p>
     * If the resulting reRooted Phylogeny is to be used any further, in most
     * cases the following three methods have to be called on the resulting
     * Phylogeny:
     * <ul>
     * <li>recalculateNumberOfExternalDescendants(boolean) <li>recalculateAndReset()
     * </ul>
     * <p>
     * (Last modified: 10/01/01)
     * 
     * @param n
     *            PhylogenyNode of this Phylogeny\
     */
    public void reRoot( final PhylogenyNode n ) {
        reRoot( n, -1 );
    }

    public void reRoot( final PhylogenyNode n, final double distance_n_to_parent ) {
        if ( isEmpty() || ( getNumberOfExternalNodes() < 2 ) ) {
            return;
        }
        setRooted( true );
        if ( n.isRoot() ) {
            return;
        }
        else if ( n.getParent().isRoot() ) {
            if ( ( n.getParent().getNumberOfDescendants() == 2 ) && ( distance_n_to_parent >= 0 ) ) {
                final double d = n.getParent().getChildNode1().getDistanceToParent()
                        + n.getParent().getChildNode2().getDistanceToParent();
                PhylogenyNode other;
                if ( n.getChildNodeIndex() == 0 ) {
                    other = n.getParent().getChildNode2();
                }
                else {
                    other = n.getParent().getChildNode1();
                }
                n.setDistanceToParent( distance_n_to_parent );
                final double dm = d - distance_n_to_parent;
                if ( dm >= 0 ) {
                    other.setDistanceToParent( dm );
                }
                else {
                    other.setDistanceToParent( 0 );
                }
            }
            if ( n.getParent().getNumberOfDescendants() > 2 ) {
                final int index = n.getChildNodeIndex();
                final double dn = n.getDistanceToParent();
                final PhylogenyNode prev_root = getRoot();
                prev_root.getDescendants().remove( index );
                final PhylogenyNode new_root = new PhylogenyNode();
                new_root.setChildNode( 0, n );
                new_root.setChildNode( 1, prev_root );
                if ( n.getBranchDataDirectly() != null ) {
                    prev_root.setBranchData( ( BranchData ) n.getBranchDataDirectly().copy() );
                }
                setRoot( new_root );
                if ( distance_n_to_parent >= 0 ) {
                    n.setDistanceToParent( distance_n_to_parent );
                    final double d = dn - distance_n_to_parent;
                    if ( d >= 0 ) {
                        prev_root.setDistanceToParent( d );
                    }
                    else {
                        prev_root.setDistanceToParent( 0 );
                    }
                }
                else {
                    if ( dn >= 0 ) {
                        final double d = dn / 2.0;
                        n.setDistanceToParent( d );
                        prev_root.setDistanceToParent( d );
                    }
                }
            }
        }
        else {
            PhylogenyNode a = n;
            PhylogenyNode b = null;
            PhylogenyNode c = null;
            final PhylogenyNode new_root = new PhylogenyNode();
            double distance1 = 0.0;
            double distance2 = 0.0;
            BranchData branch_data_1 = null;
            BranchData branch_data_2 = null;
            b = a.getParent();
            c = b.getParent();
            new_root.setChildNode( 0, a );
            new_root.setChildNode( 1, b );
            distance1 = c.getDistanceToParent();
            if ( c.getBranchDataDirectly() != null ) {
                branch_data_1 = ( BranchData ) c.getBranchDataDirectly().copy();
            }
            c.setDistanceToParent( b.getDistanceToParent() );
            if ( b.getBranchDataDirectly() != null ) {
                c.setBranchData( ( BranchData ) b.getBranchDataDirectly().copy() );
            }
            if ( a.getBranchDataDirectly() != null ) {
                b.setBranchData( ( BranchData ) a.getBranchDataDirectly().copy() );
            }
            // New root is always placed in the middle of the branch:
            if ( a.getDistanceToParent() == PhylogenyNode.DISTANCE_DEFAULT ) {
                b.setDistanceToParent( PhylogenyNode.DISTANCE_DEFAULT );
            }
            else {
                if ( distance_n_to_parent >= 0.0 ) {
                    final double diff = a.getDistanceToParent() - distance_n_to_parent;
                    a.setDistanceToParent( distance_n_to_parent );
                    b.setDistanceToParent( diff >= 0.0 ? diff : 0.0 );
                }
                else {
                    final double d = a.getDistanceToParent() / 2.0;
                    a.setDistanceToParent( d );
                    b.setDistanceToParent( d );
                }
            }
            b.setChildNodeOnly( a.getChildNodeIndex( b ), c );
            // moving to the old root, swapping references:
            while ( !c.isRoot() ) {
                a = b;
                b = c;
                c = c.getParent();
                b.setChildNodeOnly( a.getChildNodeIndex( b ), c );
                b.setParent( a );
                distance2 = c.getDistanceToParent();
                branch_data_2 = c.getBranchDataDirectly();
                c.setDistanceToParent( distance1 );
                c.setBranchData( branch_data_1 );
                distance1 = distance2;
                branch_data_1 = branch_data_2;
            }
            // removing the old root:
            if ( c.getNumberOfDescendants() == 2 ) {
                final PhylogenyNode node = c.getChildNode( 1 - b.getChildNodeIndex( c ) );
                node.setParent( b );
                if ( ( c.getDistanceToParent() == PhylogenyNode.DISTANCE_DEFAULT )
                        && ( node.getDistanceToParent() == PhylogenyNode.DISTANCE_DEFAULT ) ) {
                    node.setDistanceToParent( PhylogenyNode.DISTANCE_DEFAULT );
                }
                else {
                    node.setDistanceToParent( ( c.getDistanceToParent() >= 0.0 ? c.getDistanceToParent() : 0.0 )
                            + ( node.getDistanceToParent() >= 0.0 ? node.getDistanceToParent() : 0.0 ) );
                }
                if ( c.getBranchDataDirectly() != null ) {
                    node.setBranchData( ( BranchData ) c.getBranchDataDirectly().copy() );
                }
                for( int i = 0; i < b.getNumberOfDescendants(); ++i ) {
                    if ( b.getChildNode( i ) == c ) {
                        b.setChildNodeOnly( i, node );
                        break;
                    }
                }
            }
            else {
                c.setParent( b );
                c.removeChildNode( b.getChildNodeIndex( c ) );
            }
            setRoot( new_root );
        }
    }

    /**
     * Sets all Nodes of this Phylogeny to not-collapsed.
     * <p>
     * In most cases methods adjustNodeCount(false) and recalculateAndReset()
     * need to be called after this method has been called.
     */
    public void setAllNodesToNotCollapse() {
        if ( isEmpty() ) {
            return;
        }
        for( final PhylogenyNodeIterator iter = iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            node.setCollapse( false );
        }
    }

    private void setAllowMultipleParents( final boolean allow_multiple_parents ) {
        _allow_multiple_parents = allow_multiple_parents;
    }

    public void setConfidence( final Confidence confidence ) {
        _confidence = confidence;
    }

    public void setDescription( final String description ) {
        _description = description;
    }

    public void setDistanceUnit( final String _distance_unit ) {
        this._distance_unit = _distance_unit;
    }

    public void setIdentifier( final Identifier identifier ) {
        _identifier = identifier;
    }

    void setIdHash( final HashMap<Integer, PhylogenyNode> idhash ) {
        _idhash = idhash;
    }

    /**
     * Sets the indicators of all Nodes of this Phylogeny to 0.
     */
    public void setIndicatorsToZero() {
        if ( isEmpty() ) {
            return;
        }
        for( final PhylogenyNodeIterator iter = iteratorPreorder(); iter.hasNext(); ) {
            iter.next().setIndicator( ( byte ) 0 );
        }
    } // setIndicatorsToZero()

    /**
     * Sets the name of this Phylogeny to s.
     */
    public void setName( final String s ) {
        _name = s;
    }

    public void setRerootable( final boolean rerootable ) {
        _rerootable = rerootable;
    }

    public void setRoot( final PhylogenyNode n ) {
        _root = n;
    } // setRoot( PhylogenyNode )

    /**
     * Sets whether this Phylogeny is rooted or not.
     */
    public void setRooted( final boolean b ) {
        _rooted = b;
    } // setRooted( boolean )

    public void setType( final String type ) {
        _type = type;
    }

    public Phylogeny subTree( final PhylogenyNode node ) throws IllegalStateException {
        if ( !isTree() ) {
            throw new IllegalStateException( "attempt to get sub tree on phylogeny which is not tree-like." );
        }
        if ( ( node == null ) || isEmpty() ) {
            return null;
        }
        Phylogeny sub_tree = null;
        sub_tree = copy();
        final PhylogenyNode new_root = sub_tree.getNode( node.getNodeId() );
        new_root.setParent( null );
        //node.setDistanceToParent( PhylogenyNode.DISTANCE_DEFAULT );
        sub_tree.setName( "" );
        sub_tree.setDescription( "" );
        sub_tree.setIdHash( null );
        sub_tree.setConfidence( null );
        sub_tree.setIdentifier( null );
        sub_tree.setRooted( true );
        sub_tree.setRoot( new_root );
        sub_tree.recalculateNumberOfExternalDescendants( true );
        return sub_tree;
    }

    /**
     * Swaps the the two childern of a PhylogenyNode node of this Phylogeny.
     * <p>
     * (Last modified: 06/13/01)
     * 
     * @param node
     *            a PhylogenyNode of this Phylogeny
     */
    public void swapChildren( final PhylogenyNode node ) throws IllegalStateException {
        if ( !isTree() ) {
            throw new IllegalStateException( "Attempt to swap children on phylogeny which is not tree-like." );
        }
        if ( isEmpty() || node.isExternal() || ( node.getNumberOfDescendants() < 2 ) ) {
            return;
        }
        final PhylogenyNode first = node.getFirstChildNode();
        for( int i = 1; i < node.getNumberOfDescendants(); ++i ) {
            node.setChildNode( i - 1, node.getChildNode( i ) );
        }
        node.setChildNode( node.getNumberOfDescendants() - 1, first );
    } // swapChildren( PhylogenyNode )

    public String toNewHampshire( final boolean simple_nh ) {
        try {
            return new PhylogenyWriter().toNewHampshire( this, simple_nh, true ).toString();
        }
        catch ( final IOException e ) {
            throw new Error( "this should not have happend: " + e.getMessage() );
        }
    }

    public String toNewHampshireX() {
        try {
            return new PhylogenyWriter().toNewHampshireX( this ).toString();
        }
        catch ( final IOException e ) {
            throw new Error( "this should not have happend: " + e.getMessage() );
        }
    }

    public String toNexus() {
        try {
            return new PhylogenyWriter().toNexus( this ).toString();
        }
        catch ( final IOException e ) {
            throw new Error( "this should not have happend: " + e.getMessage() );
        }
    }

    public String toPhyloXML( final int phyloxml_level ) {
        try {
            return new PhylogenyWriter().toPhyloXML( this, phyloxml_level ).toString();
        }
        catch ( final IOException e ) {
            throw new Error( "this should not have happend: " + e.getMessage() );
        }
    }

    // ---------------------------------------------------------
    // Writing of Phylogeny to Strings
    // ---------------------------------------------------------
    /**
     * Converts this Phylogeny to a New Hampshire X (String) representation.
     * 
     * @return New Hampshire X (String) representation of this
     * @see #toNewHampshireX()
     */
    @Override
    public String toString() {
        return toNewHampshireX();
    }

    /**
     * Removes the root PhylogenyNode this Phylogeny.
     */
    public void unRoot() throws IllegalStateException {
        if ( !isTree() ) {
            throw new IllegalStateException( "Attempt to unroot a phylogeny which is not tree-like." );
        }
        if ( isEmpty() ) {
            return;
        }
        setIndicatorsToZero();
        if ( !isRooted() || ( getNumberOfExternalNodes() <= 1 ) ) {
            return;
        }
        setRooted( false );
        return;
    } // unRoot()
}
