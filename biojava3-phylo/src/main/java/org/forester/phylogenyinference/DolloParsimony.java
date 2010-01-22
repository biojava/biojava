// $Id: DolloParsimony.java,v 1.15 2009/10/26 23:29:39 cmzmasek Exp $
//
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

package org.forester.phylogenyinference;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogenyinference.CharacterStateMatrix.BinaryStates;
import org.forester.phylogenyinference.CharacterStateMatrix.GainLossStates;
import org.forester.util.ForesterUtil;

public class DolloParsimony {

    final static private BinaryStates            PRESENT                         = BinaryStates.PRESENT;
    final static private BinaryStates            ABSENT                          = BinaryStates.ABSENT;
    final static private BinaryStates            UNKNOWN                         = BinaryStates.UNKNOWN;
    final static private GainLossStates          LOSS                            = GainLossStates.LOSS;
    final static private GainLossStates          GAIN                            = GainLossStates.GAIN;
    final static private GainLossStates          UNCHANGED_PRESENT               = GainLossStates.UNCHANGED_PRESENT;
    final static private GainLossStates          UNCHANGED_ABSENT                = GainLossStates.UNCHANGED_ABSENT;
    private static final boolean                 RETURN_INTERNAL_STATES_DEFAULT  = false;
    private static final boolean                 RETURN_GAIN_LOSS_MATRIX_DEFAULT = false;
    private boolean                              _return_internal_states         = false;
    private boolean                              _return_gain_loss               = false;
    private int                                  _total_gains;
    private int                                  _total_losses;
    private int                                  _total_unchanged;
    private CharacterStateMatrix<BinaryStates>   _internal_states_matrix;
    private CharacterStateMatrix<GainLossStates> _gain_loss_matrix;

    private DolloParsimony() {
        init();
    }

    public void execute( final Phylogeny p, final CharacterStateMatrix<BinaryStates> external_node_states_matrix ) {
        if ( !p.isRooted() ) {
            throw new IllegalArgumentException( "attempt to execute Dollo parsimony on unroored phylogeny" );
        }
        if ( external_node_states_matrix.isEmpty() ) {
            throw new IllegalArgumentException( "character matrix is empty" );
        }
        if ( external_node_states_matrix.getNumberOfIdentifiers() != p.getNumberOfExternalNodes() ) {
            throw new IllegalArgumentException( "number of external nodes in phylogeny ["
                    + p.getNumberOfExternalNodes() + "] and number of indentifiers ["
                    + external_node_states_matrix.getNumberOfIdentifiers() + "] in matrix are not equal" );
        }
        reset();
        if ( isReturnInternalStates() ) {
            initializeInternalStates( p, external_node_states_matrix );
        }
        if ( isReturnGainLossMatrix() ) {
            initializeGainLossMatrix( p, external_node_states_matrix );
        }
        for( int character_index = 0; character_index < external_node_states_matrix.getNumberOfCharacters(); ++character_index ) {
            executeForOneCharacter( p,
                                    getStatesForCharacter( p, external_node_states_matrix, character_index ),
                                    character_index );
        }
        if ( ( external_node_states_matrix.getNumberOfCharacters() * p.getNumberOfBranches() ) != ( getTotalGains()
                + getTotalLosses() + getTotalUnchanged() ) ) {
            throw new AssertionError( "this should not have happened: something is deeply wrong with Dollo parsimony implementation" );
        }
    }

    private void executeForOneCharacter( final Phylogeny p,
                                         final Map<PhylogenyNode, BinaryStates> states,
                                         final int character_state_column ) {
        postOrderTraversal( p, states );
        preOrderTraversal( p, states, character_state_column );
    }

    /* (non-Javadoc)
     * @see org.forester.phylogenyinference.Parsimony#getCost()
     */
    public int getCost() {
        return getTotalGains() + getTotalLosses();
    }

    /* (non-Javadoc)
     * @see org.forester.phylogenyinference.Parsimony#getGainLossMatrix()
     */
    public CharacterStateMatrix<CharacterStateMatrix.GainLossStates> getGainLossMatrix() {
        if ( !isReturnGainLossMatrix() ) {
            throw new IllegalStateException( "creation of gain-loss matrix has not been enabled" );
        }
        return _gain_loss_matrix;
    }

    /* (non-Javadoc)
     * @see org.forester.phylogenyinference.Parsimony#getInternalStatesMatrix()
     */
    public CharacterStateMatrix<BinaryStates> getInternalStatesMatrix() {
        if ( !isReturnInternalStates() ) {
            throw new IllegalStateException( "creation of internal state matrix has not been enabled" );
        }
        return _internal_states_matrix;
    }

    private Map<PhylogenyNode, BinaryStates> getStatesForCharacter( final Phylogeny p,
                                                                    final CharacterStateMatrix<BinaryStates> matrix,
                                                                    final int character_index ) {
        final Map<PhylogenyNode, BinaryStates> states = new HashMap<PhylogenyNode, BinaryStates>( matrix
                .getNumberOfIdentifiers() );
        for( int indentifier_index = 0; indentifier_index < matrix.getNumberOfIdentifiers(); ++indentifier_index ) {
            final BinaryStates state = matrix.getState( indentifier_index, character_index );
            if ( state == null ) {
                throw new IllegalArgumentException( "value at [" + indentifier_index + ", " + character_index
                        + "] is null" );
            }
            states.put( p.getNode( matrix.getIdentifier( indentifier_index ) ), state );
        }
        return states;
    }

    /* (non-Javadoc)
     * @see org.forester.phylogenyinference.Parsimony#getTotalGains()
     */
    public int getTotalGains() {
        return _total_gains;
    }

    /* (non-Javadoc)
     * @see org.forester.phylogenyinference.Parsimony#getTotalLosses()
     */
    public int getTotalLosses() {
        return _total_losses;
    }

    /* (non-Javadoc)
     * @see org.forester.phylogenyinference.Parsimony#getTotalUnchanged()
     */
    public int getTotalUnchanged() {
        return _total_unchanged;
    }

    private void init() {
        setReturnInternalStates( RETURN_INTERNAL_STATES_DEFAULT );
        setReturnGainLossMatrix( RETURN_GAIN_LOSS_MATRIX_DEFAULT );
        reset();
    }

    private void initializeGainLossMatrix( final Phylogeny p,
                                           final CharacterStateMatrix<BinaryStates> external_node_states_matrix ) {
        final List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator postorder = p.iteratorPostorder(); postorder.hasNext(); ) {
            nodes.add( postorder.next() );
        }
        setGainLossMatrix( new BasicCharacterStateMatrix<CharacterStateMatrix.GainLossStates>( nodes.size(),
                                                                                               external_node_states_matrix
                                                                                                       .getNumberOfCharacters() ) );
        int identifier_index = 0;
        for( final PhylogenyNode node : nodes ) {
            getGainLossMatrix().setIdentifier( identifier_index++,
                                               ForesterUtil.isEmpty( node.getNodeName() ) ? node.getNodeId() + ""
                                                       : node.getNodeName() );
        }
        for( int character_index = 0; character_index < external_node_states_matrix.getNumberOfCharacters(); ++character_index ) {
            getGainLossMatrix().setCharacter( character_index,
                                              external_node_states_matrix.getCharacter( character_index ) );
        }
    }

    private void initializeInternalStates( final Phylogeny p,
                                           final CharacterStateMatrix<BinaryStates> external_node_states_matrix ) {
        final List<PhylogenyNode> internal_nodes = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator postorder = p.iteratorPostorder(); postorder.hasNext(); ) {
            final PhylogenyNode node = postorder.next();
            if ( node.isInternal() ) {
                internal_nodes.add( node );
            }
        }
        setInternalStatesMatrix( new BasicCharacterStateMatrix<BinaryStates>( internal_nodes.size(),
                                                                              external_node_states_matrix
                                                                                      .getNumberOfCharacters() ) );
        int identifier_index = 0;
        for( final PhylogenyNode node : internal_nodes ) {
            getInternalStatesMatrix().setIdentifier( identifier_index++,
                                                     ForesterUtil.isEmpty( node.getNodeName() ) ? node.getNodeId() + ""
                                                             : node.getNodeName() );
        }
        for( int character_index = 0; character_index < external_node_states_matrix.getNumberOfCharacters(); ++character_index ) {
            getInternalStatesMatrix().setCharacter( character_index,
                                                    external_node_states_matrix.getCharacter( character_index ) );
        }
    }

    private boolean isReturnGainLossMatrix() {
        return _return_gain_loss;
    }

    private boolean isReturnInternalStates() {
        return _return_internal_states;
    }

    private void postOrderTraversal( final Phylogeny p, final Map<PhylogenyNode, BinaryStates> states )
            throws AssertionError {
        for( final PhylogenyNodeIterator postorder = p.iteratorPostorder(); postorder.hasNext(); ) {
            final PhylogenyNode node = postorder.next();
            if ( !node.isExternal() ) {
                final int present_unknown = getNumberOfChildNodesWithPresentOrUnknownState( states, node );
                if ( present_unknown < 1 ) {
                    states.put( node, ABSENT );
                }
                else if ( present_unknown == 1 ) {
                    states.put( node, UNKNOWN );
                }
                else {
                    states.put( node, PRESENT );
                }
            }
        }
    }

    private void preOrderTraversal( final Phylogeny p,
                                    final Map<PhylogenyNode, BinaryStates> states,
                                    final int character_state_column ) throws AssertionError {
        boolean gain = false;
        for( final PhylogenyNodeIterator preorder = p.iteratorPreorder(); preorder.hasNext(); ) {
            final PhylogenyNode node = preorder.next();
            BinaryStates parent_state = null;
            if ( !node.isRoot() ) {
                parent_state = states.get( node.getParent() );
            }
            if ( !node.isExternal() ) {
                if ( states.get( node ) == UNKNOWN ) {
                    if ( parent_state == PRESENT ) {
                        states.put( node, PRESENT );
                    }
                    else {
                        if ( isCharacterPresentOrUnknownInAtLeastTwoChildNodes( states, node ) ) {
                            states.put( node, PRESENT );
                        }
                        else {
                            states.put( node, ABSENT );
                        }
                    }
                }
                if ( isReturnInternalStates() ) {
                    setInternalNodeState( states, character_state_column, node );
                }
            }
            final BinaryStates current_state = states.get( node );
            if ( ( parent_state == PRESENT ) && ( current_state == ABSENT ) ) {
                ++_total_losses;
                if ( isReturnGainLossMatrix() ) {
                    setGainLossState( character_state_column, node, LOSS );
                }
            }
            else if ( ( ( parent_state == ABSENT ) ) && ( current_state == PRESENT ) ) {
                if ( gain ) {
                    throw new RuntimeException( "this should not have happened: dollo parsimony cannot have more than one gain" );
                }
                gain = true;
                ++_total_gains;
                if ( isReturnGainLossMatrix() ) {
                    setGainLossState( character_state_column, node, GAIN );
                }
            }
            else {
                ++_total_unchanged;
                if ( isReturnGainLossMatrix() ) {
                    if ( current_state == PRESENT ) {
                        setGainLossState( character_state_column, node, UNCHANGED_PRESENT );
                    }
                    else if ( current_state == ABSENT ) {
                        setGainLossState( character_state_column, node, UNCHANGED_ABSENT );
                    }
                }
            }
        }
    }

    private void reset() {
        setTotalLosses( 0 );
        setTotalGains( 0 );
        setTotalUnchanged( 0 );
    }

    private void setGainLossMatrix( final CharacterStateMatrix<CharacterStateMatrix.GainLossStates> gain_loss_matrix ) {
        _gain_loss_matrix = gain_loss_matrix;
    }

    private void setGainLossState( final int character_state_column,
                                   final PhylogenyNode node,
                                   final GainLossStates state ) {
        getGainLossMatrix().setState( ForesterUtil.isEmpty( node.getNodeName() ) ? node.getNodeId() + "" : node
                                              .getNodeName(),
                                      character_state_column,
                                      state );
    }

    private void setInternalNodeState( final Map<PhylogenyNode, BinaryStates> states,
                                       final int character_state_column,
                                       final PhylogenyNode node ) {
        getInternalStatesMatrix().setState( ForesterUtil.isEmpty( node.getNodeName() ) ? node.getNodeId() + "" : node
                                                    .getNodeName(),
                                            character_state_column,
                                            states.get( node ) );
    }

    private void setInternalStatesMatrix( final CharacterStateMatrix<BinaryStates> internal_states_matrix ) {
        _internal_states_matrix = internal_states_matrix;
    }

    /* (non-Javadoc)
     * @see org.forester.phylogenyinference.Parsimony#setReturnGainLossMatrix(boolean)
     */
    public void setReturnGainLossMatrix( final boolean return_gain_loss ) {
        _return_gain_loss = return_gain_loss;
    }

    /* (non-Javadoc)
     * @see org.forester.phylogenyinference.Parsimony#setReturnInternalStates(boolean)
     */
    public void setReturnInternalStates( final boolean return_internal_states ) {
        _return_internal_states = return_internal_states;
    }

    private void setTotalGains( final int total_gains ) {
        _total_gains = total_gains;
    }

    private void setTotalLosses( final int total_losses ) {
        _total_losses = total_losses;
    }

    private void setTotalUnchanged( final int total_unchanged ) {
        _total_unchanged = total_unchanged;
    }

    public static DolloParsimony createInstance() {
        return new DolloParsimony();
    }

    private static int getNumberOfChildNodesWithPresentOrUnknownState( final Map<PhylogenyNode, BinaryStates> states,
                                                                       final PhylogenyNode node ) {
        int presents = 0;
        for( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
            final PhylogenyNode node_child = node.getChildNode( i );
            if ( !states.containsKey( node_child ) ) {
                throw new RuntimeException( "this should not have happened: node [" + node_child.getNodeName()
                        + "] not found in node state map" );
            }
            if ( ( states.get( node_child ) == BinaryStates.PRESENT )
                    || ( states.get( node_child ) == BinaryStates.UNKNOWN ) ) {
                ++presents;
            }
        }
        return presents;
    }

    private static boolean isCharacterPresentOrUnknownInAtLeastTwoChildNodes( final Map<PhylogenyNode, BinaryStates> states,
                                                                              final PhylogenyNode node ) {
        int count = 0;
        for( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
            final PhylogenyNode node_child = node.getChildNode( i );
            if ( ( states.get( node_child ) == PRESENT ) || ( states.get( node_child ) == UNKNOWN ) ) {
                ++count;
                if ( count > 1 ) {
                    return true;
                }
            }
        }
        return false;
    }
}
