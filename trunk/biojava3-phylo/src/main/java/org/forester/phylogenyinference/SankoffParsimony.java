// $Id: SankoffParsimony.java,v 1.5 2009/10/26 23:29:39 cmzmasek Exp $
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
// WWW: www.phylosoft.org

package org.forester.phylogenyinference;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogenyinference.CharacterStateMatrix.BinaryStates;
import org.forester.phylogenyinference.CharacterStateMatrix.GainLossStates;
import org.forester.util.ForesterUtil;

public class SankoffParsimony<STATE_TYPE> {

    final static private BinaryStates              PRESENT                         = BinaryStates.PRESENT;
    final static private BinaryStates              ABSENT                          = BinaryStates.ABSENT;
    final static private GainLossStates            LOSS                            = GainLossStates.LOSS;
    final static private GainLossStates            GAIN                            = GainLossStates.GAIN;
    final static private GainLossStates            UNCHANGED_PRESENT               = GainLossStates.UNCHANGED_PRESENT;
    final static private GainLossStates            UNCHANGED_ABSENT                = GainLossStates.UNCHANGED_ABSENT;
    private static final boolean                   RETURN_INTERNAL_STATES_DEFAULT  = false;
    private static final boolean                   RETURN_GAIN_LOSS_MATRIX_DEFAULT = false;
    private static final boolean                   RANDOMIZE_DEFAULT               = false;
    private static final long                      RANDOM_NUMBER_SEED_DEFAULT      = 21;
    private static final boolean                   USE_LAST_DEFAULT                = false;
    private boolean                                _return_internal_states         = false;
    private boolean                                _return_gain_loss               = false;
    private int                                    _total_gains;
    private int                                    _total_losses;
    private int                                    _total_unchanged;
    private CharacterStateMatrix<List<STATE_TYPE>> _internal_states_matrix_prior_to_traceback;
    private CharacterStateMatrix<STATE_TYPE>       _internal_states_matrix_after_traceback;
    private CharacterStateMatrix<GainLossStates>   _gain_loss_matrix;
    private boolean                                _randomize;
    private boolean                                _use_last;
    private int                                    _cost;
    private long                                   _random_number_seed;
    private Random                                 _random_generator;

    public SankoffParsimony() {
        init();
    }

    private int determineIndex( final SortedSet<STATE_TYPE> current_node_states, int i ) {
        if ( isRandomize() ) {
            i = getRandomGenerator().nextInt( current_node_states.size() );
        }
        else if ( isUseLast() ) {
            i = current_node_states.size() - 1;
        }
        return i;
    }

    public void execute( final Phylogeny p, final CharacterStateMatrix<STATE_TYPE> external_node_states_matrix ) {
        if ( !p.isRooted() ) {
            throw new IllegalArgumentException( "attempt to execute Fitch parsimony on unroored phylogeny" );
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
                                    getStatesForCharacterForTraceback( p, external_node_states_matrix, character_index ),
                                    character_index );
        }
        if ( external_node_states_matrix.getState( 0, 0 ) instanceof BinaryStates ) {
            if ( ( external_node_states_matrix.getNumberOfCharacters() * p.getNumberOfBranches() ) != ( getTotalGains()
                    + getTotalLosses() + getTotalUnchanged() ) ) {
                throw new IllegalStateException( "this should not have happened: something is deeply wrong with Fitch parsimony implementation" );
            }
        }
    }

    private void executeForOneCharacter( final Phylogeny p,
                                         final Map<PhylogenyNode, SortedSet<STATE_TYPE>> states,
                                         final Map<PhylogenyNode, STATE_TYPE> traceback_states,
                                         final int character_state_column ) {
        postOrderTraversal( p, states );
        preOrderTraversal( p, states, traceback_states, character_state_column );
    }

    public int getCost() {
        return _cost;
    }

    public CharacterStateMatrix<CharacterStateMatrix.GainLossStates> getGainLossMatrix() {
        if ( !isReturnGainLossMatrix() ) {
            throw new IllegalStateException( "creation of gain-loss matrix has not been enabled" );
        }
        return _gain_loss_matrix;
    }

    public CharacterStateMatrix<STATE_TYPE> getInternalStatesMatrix() {
        if ( !isReturnInternalStates() ) {
            throw new IllegalStateException( "creation of internal state matrix has not been enabled" );
        }
        return _internal_states_matrix_after_traceback;
    }

    /**
     * Returns a view of the internal states prior to trace-back.
     * 
     * @return
     */
    public CharacterStateMatrix<List<STATE_TYPE>> getInternalStatesMatrixPriorToTraceback() {
        if ( !isReturnInternalStates() ) {
            throw new IllegalStateException( "creation of internal state matrix has not been enabled" );
        }
        return _internal_states_matrix_prior_to_traceback;
    }

    private SortedSet<STATE_TYPE> getIntersectionOfStatesOfChildNodes( final Map<PhylogenyNode, SortedSet<STATE_TYPE>> states,
                                                                       final PhylogenyNode node ) throws AssertionError {
        final SortedSet<STATE_TYPE> states_in_child_nodes = new TreeSet<STATE_TYPE>();
        for( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
            final PhylogenyNode node_child = node.getChildNode( i );
            if ( !states.containsKey( node_child ) ) {
                throw new AssertionError( "this should not have happened: node [" + node_child.getNodeName()
                        + "] not found in node state map" );
            }
            if ( i == 0 ) {
                states_in_child_nodes.addAll( states.get( node_child ) );
            }
            else {
                states_in_child_nodes.retainAll( states.get( node_child ) );
            }
        }
        return states_in_child_nodes;
    }

    private Random getRandomGenerator() {
        return _random_generator;
    }

    private long getRandomNumberSeed() {
        return _random_number_seed;
    }

    private STATE_TYPE getStateAt( final int i, final SortedSet<STATE_TYPE> states ) {
        final Iterator<STATE_TYPE> it = states.iterator();
        for( int j = 0; j < i; ++j ) {
            it.next();
        }
        return it.next();
    }

    private Map<PhylogenyNode, SortedSet<STATE_TYPE>> getStatesForCharacter( final Phylogeny p,
                                                                             final CharacterStateMatrix<STATE_TYPE> matrix,
                                                                             final int character_index ) {
        final Map<PhylogenyNode, SortedSet<STATE_TYPE>> states = new HashMap<PhylogenyNode, SortedSet<STATE_TYPE>>( matrix
                .getNumberOfIdentifiers() );
        for( int indentifier_index = 0; indentifier_index < matrix.getNumberOfIdentifiers(); ++indentifier_index ) {
            final STATE_TYPE state = matrix.getState( indentifier_index, character_index );
            if ( state == null ) {
                throw new IllegalArgumentException( "value at [" + indentifier_index + ", " + character_index
                        + "] is null" );
            }
            final SortedSet<STATE_TYPE> l = new TreeSet<STATE_TYPE>();
            l.add( state );
            states.put( p.getNode( matrix.getIdentifier( indentifier_index ) ), l );
        }
        return states;
    }

    private Map<PhylogenyNode, STATE_TYPE> getStatesForCharacterForTraceback( final Phylogeny p,
                                                                              final CharacterStateMatrix<STATE_TYPE> matrix,
                                                                              final int character_index ) {
        final Map<PhylogenyNode, STATE_TYPE> states = new HashMap<PhylogenyNode, STATE_TYPE>( matrix
                .getNumberOfIdentifiers() );
        for( int indentifier_index = 0; indentifier_index < matrix.getNumberOfIdentifiers(); ++indentifier_index ) {
            final STATE_TYPE state = matrix.getState( indentifier_index, character_index );
            if ( state == null ) {
                throw new IllegalArgumentException( "value at [" + indentifier_index + ", " + character_index
                        + "] is null" );
            }
            states.put( p.getNode( matrix.getIdentifier( indentifier_index ) ), state );
        }
        return states;
    }

    public int getTotalGains() {
        return _total_gains;
    }

    public int getTotalLosses() {
        return _total_losses;
    }

    public int getTotalUnchanged() {
        return _total_unchanged;
    }

    private SortedSet<STATE_TYPE> getUnionOfStatesOfChildNodes( final Map<PhylogenyNode, SortedSet<STATE_TYPE>> states,
                                                                final PhylogenyNode node ) throws AssertionError {
        final SortedSet<STATE_TYPE> states_in_child_nodes = new TreeSet<STATE_TYPE>();
        for( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
            final PhylogenyNode node_child = node.getChildNode( i );
            if ( !states.containsKey( node_child ) ) {
                throw new AssertionError( "this should not have happened: node [" + node_child.getNodeName()
                        + "] not found in node state map" );
            }
            states_in_child_nodes.addAll( states.get( node_child ) );
        }
        return states_in_child_nodes;
    }

    private void increaseCost() {
        ++_cost;
    }

    private void init() {
        setReturnInternalStates( RETURN_INTERNAL_STATES_DEFAULT );
        setReturnGainLossMatrix( RETURN_GAIN_LOSS_MATRIX_DEFAULT );
        setRandomize( RANDOMIZE_DEFAULT );
        setUseLast( USE_LAST_DEFAULT );
        _random_number_seed = RANDOM_NUMBER_SEED_DEFAULT;
        reset();
    }

    private void initializeGainLossMatrix( final Phylogeny p,
                                           final CharacterStateMatrix<STATE_TYPE> external_node_states_matrix ) {
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
                                           final CharacterStateMatrix<STATE_TYPE> external_node_states_matrix ) {
        final List<PhylogenyNode> internal_nodes = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator postorder = p.iteratorPostorder(); postorder.hasNext(); ) {
            final PhylogenyNode node = postorder.next();
            if ( node.isInternal() ) {
                internal_nodes.add( node );
            }
        }
        setInternalStatesMatrixPriorToTraceback( new BasicCharacterStateMatrix<List<STATE_TYPE>>( internal_nodes.size(),
                                                                                                  external_node_states_matrix
                                                                                                          .getNumberOfCharacters() ) );
        setInternalStatesMatrixTraceback( new BasicCharacterStateMatrix<STATE_TYPE>( internal_nodes.size(),
                                                                                     external_node_states_matrix
                                                                                             .getNumberOfCharacters() ) );
        int identifier_index = 0;
        for( final PhylogenyNode node : internal_nodes ) {
            getInternalStatesMatrix().setIdentifier( identifier_index,
                                                     ForesterUtil.isEmpty( node.getNodeName() ) ? node.getNodeId() + ""
                                                             : node.getNodeName() );
            getInternalStatesMatrixPriorToTraceback().setIdentifier( identifier_index,
                                                                     ForesterUtil.isEmpty( node.getNodeName() ) ? node
                                                                             .getNodeId()
                                                                             + "" : node.getNodeName() );
            ++identifier_index;
        }
        for( int character_index = 0; character_index < external_node_states_matrix.getNumberOfCharacters(); ++character_index ) {
            getInternalStatesMatrix().setCharacter( character_index,
                                                    external_node_states_matrix.getCharacter( character_index ) );
            getInternalStatesMatrixPriorToTraceback().setCharacter( character_index,
                                                                    external_node_states_matrix
                                                                            .getCharacter( character_index ) );
        }
    }

    private boolean isRandomize() {
        return _randomize;
    }

    private boolean isReturnGainLossMatrix() {
        return _return_gain_loss;
    }

    private boolean isReturnInternalStates() {
        return _return_internal_states;
    }

    private boolean isUseLast() {
        return _use_last;
    }

    private void postOrderTraversal( final Phylogeny p, final Map<PhylogenyNode, SortedSet<STATE_TYPE>> states )
            throws AssertionError {
        for( final PhylogenyNodeIterator postorder = p.iteratorPostorder(); postorder.hasNext(); ) {
            final PhylogenyNode node = postorder.next();
            if ( !node.isExternal() ) {
                SortedSet<STATE_TYPE> states_in_children = getIntersectionOfStatesOfChildNodes( states, node );
                if ( states_in_children.isEmpty() ) {
                    states_in_children = getUnionOfStatesOfChildNodes( states, node );
                }
                states.put( node, states_in_children );
            }
        }
    }

    private void preOrderTraversal( final Phylogeny p,
                                    final Map<PhylogenyNode, SortedSet<STATE_TYPE>> states,
                                    final Map<PhylogenyNode, STATE_TYPE> traceback_states,
                                    final int character_state_column ) throws AssertionError {
        for( final PhylogenyNodeIterator preorder = p.iteratorPreorder(); preorder.hasNext(); ) {
            final PhylogenyNode current_node = preorder.next();
            final SortedSet<STATE_TYPE> current_node_states = states.get( current_node );
            STATE_TYPE parent_state = null;
            if ( current_node.isRoot() ) {
                int i = 0;
                i = determineIndex( current_node_states, i );
                traceback_states.put( current_node, getStateAt( i, current_node_states ) );
            }
            else {
                parent_state = traceback_states.get( current_node.getParent() );
                if ( current_node_states.contains( parent_state ) ) {
                    traceback_states.put( current_node, parent_state );
                }
                else {
                    increaseCost();
                    int i = 0;
                    i = determineIndex( current_node_states, i );
                    traceback_states.put( current_node, getStateAt( i, current_node_states ) );
                }
            }
            if ( isReturnInternalStates() ) {
                if ( !current_node.isExternal() ) {
                    setInternalNodeStatePriorToTraceback( states, character_state_column, current_node );
                    setInternalNodeState( traceback_states, character_state_column, current_node );
                }
            }
            if ( isReturnGainLossMatrix() && !current_node.isRoot() ) {
                if ( !( parent_state instanceof BinaryStates ) ) {
                    throw new IllegalStateException( "attempt to create gain loss matrix for not binary states" );
                }
                final BinaryStates parent_binary_state = ( BinaryStates ) parent_state;
                final BinaryStates current_binary_state = ( BinaryStates ) traceback_states.get( current_node );
                if ( ( parent_binary_state == PRESENT ) && ( current_binary_state == ABSENT ) ) {
                    ++_total_losses;
                    setGainLossState( character_state_column, current_node, LOSS );
                }
                else if ( ( ( parent_binary_state == ABSENT ) || ( parent_binary_state == null ) )
                        && ( current_binary_state == PRESENT ) ) {
                    ++_total_gains;
                    setGainLossState( character_state_column, current_node, GAIN );
                }
                else {
                    ++_total_unchanged;
                    if ( current_binary_state == PRESENT ) {
                        setGainLossState( character_state_column, current_node, UNCHANGED_PRESENT );
                    }
                    else if ( current_binary_state == ABSENT ) {
                        setGainLossState( character_state_column, current_node, UNCHANGED_ABSENT );
                    }
                }
            }
            else if ( isReturnGainLossMatrix() && current_node.isRoot() ) {
                final BinaryStates current_binary_state = ( BinaryStates ) traceback_states.get( current_node );
                ++_total_unchanged; //new
                if ( current_binary_state == PRESENT ) {//new
                    setGainLossState( character_state_column, current_node, UNCHANGED_PRESENT );//new
                }//new
                else if ( current_binary_state == ABSENT ) {//new
                    setGainLossState( character_state_column, current_node, UNCHANGED_ABSENT );//new
                }//new
                // setGainLossState( character_state_column, current_node, UNKNOWN_GAIN_LOSS );
            }
        }
    }

    private void reset() {
        setCost( 0 );
        setTotalLosses( 0 );
        setTotalGains( 0 );
        setTotalUnchanged( 0 );
        setRandomGenerator( new Random( getRandomNumberSeed() ) );
    }

    private void setCost( final int cost ) {
        _cost = cost;
    }

    private void setGainLossMatrix( final CharacterStateMatrix<GainLossStates> gain_loss_matrix ) {
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

    private void setInternalNodeState( final Map<PhylogenyNode, STATE_TYPE> states,
                                       final int character_state_column,
                                       final PhylogenyNode node ) {
        getInternalStatesMatrix().setState( ForesterUtil.isEmpty( node.getNodeName() ) ? node.getNodeId() + "" : node
                                                    .getNodeName(),
                                            character_state_column,
                                            states.get( node ) );
    }

    private void setInternalNodeStatePriorToTraceback( final Map<PhylogenyNode, SortedSet<STATE_TYPE>> states,
                                                       final int character_state_column,
                                                       final PhylogenyNode node ) {
        getInternalStatesMatrixPriorToTraceback().setState( ForesterUtil.isEmpty( node.getNodeName() ) ? node
                                                                    .getNodeId()
                                                                    + "" : node.getNodeName(),
                                                            character_state_column,
                                                            toListSorted( states.get( node ) ) );
    }

    private void setInternalStatesMatrixPriorToTraceback( final CharacterStateMatrix<List<STATE_TYPE>> internal_states_matrix_prior_to_traceback ) {
        _internal_states_matrix_prior_to_traceback = internal_states_matrix_prior_to_traceback;
    }

    private void setInternalStatesMatrixTraceback( final CharacterStateMatrix<STATE_TYPE> internal_states_matrix_after_traceback ) {
        _internal_states_matrix_after_traceback = internal_states_matrix_after_traceback;
    }

    private void setRandomGenerator( final Random random_generator ) {
        _random_generator = random_generator;
    }

    public void setRandomize( final boolean randomize ) {
        if ( randomize && isUseLast() ) {
            throw new IllegalArgumentException( "attempt to allways use last state (ordered) if more than one choices and randomization at the same time" );
        }
        _randomize = randomize;
    }

    public void setRandomNumberSeed( final long random_number_seed ) {
        if ( !isRandomize() ) {
            throw new IllegalArgumentException( "attempt to set random number generator seed without randomization enabled" );
        }
        _random_number_seed = random_number_seed;
    }

    public void setReturnGainLossMatrix( final boolean return_gain_loss ) {
        _return_gain_loss = return_gain_loss;
    }

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

    /**
     * This sets whether to use the first or last state in the sorted
     * states at the undecided internal nodes.
     * For randomized choices set randomize to true (and this to false).
     * 
     * Note. It might be advisable to set this to false
     * for BinaryStates if absence at the root is preferred
     * (given the enum BinaryStates sorts in the following order: 
     * ABSENT, UNKNOWN, PRESENT).
     * 
     * 
     * @param use_last
     */
    public void setUseLast( final boolean use_last ) {
        if ( use_last && isRandomize() ) {
            throw new IllegalArgumentException( "attempt to allways use last state (ordered) if more than one choices and randomization at the same time" );
        }
        _use_last = use_last;
    }

    private List<STATE_TYPE> toListSorted( final SortedSet<STATE_TYPE> states ) {
        final List<STATE_TYPE> l = new ArrayList<STATE_TYPE>( states.size() );
        for( final STATE_TYPE state : states ) {
            l.add( state );
        }
        return l;
    }
}
