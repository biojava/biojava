// $Id: DomainParsimonyCalculator.java,v 1.24 2009/10/26 23:29:40 cmzmasek Exp $
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

package org.forester.surfacing;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.BinaryCharacters;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogenyinference.BasicCharacterStateMatrix;
import org.forester.phylogenyinference.CharacterStateMatrix;
import org.forester.phylogenyinference.DolloParsimony;
import org.forester.phylogenyinference.FitchParsimony;
import org.forester.phylogenyinference.CharacterStateMatrix.BinaryStates;
import org.forester.phylogenyinference.CharacterStateMatrix.GainLossStates;
import org.forester.surfacing.BinaryDomainCombination.DomainCombinationType;
import org.forester.util.ForesterUtil;

public final class DomainParsimonyCalculator {

    private static final String                     TYPE_FORBINARY_CHARACTERS = "parsimony inferred";
    private CharacterStateMatrix<GainLossStates>    _gain_loss_matrix;
    private CharacterStateMatrix<BinaryStates>      _binary_internal_states_matrix;
    private final List<GenomeWideCombinableDomains> _gwcd_list;
    private final Phylogeny                         _phylogeny;
    private int                                     _total_losses;
    private int                                     _total_gains;
    private int                                     _total_unchanged;
    private int                                     _cost;
    private Map<DomainId, Set<String>>              _domain_id_to_secondary_features_map;
    private SortedSet<DomainId>                     _positive_filter;

    private DomainParsimonyCalculator( final Phylogeny phylogeny ) {
        init();
        _phylogeny = phylogeny;
        _gwcd_list = null;
    }

    private DomainParsimonyCalculator( final Phylogeny phylogeny, final List<GenomeWideCombinableDomains> gwcd_list ) {
        init();
        _phylogeny = phylogeny;
        _gwcd_list = gwcd_list;
    }

    private DomainParsimonyCalculator( final Phylogeny phylogeny,
                                       final List<GenomeWideCombinableDomains> gwcd_list,
                                       final Map<DomainId, Set<String>> domain_id_to_secondary_features_map ) {
        init();
        _phylogeny = phylogeny;
        _gwcd_list = gwcd_list;
        setDomainIdToSecondaryFeaturesMap( domain_id_to_secondary_features_map );
    }

    int calculateNumberOfBinaryDomainCombination() {
        if ( getGenomeWideCombinableDomainsList().isEmpty() ) {
            throw new IllegalArgumentException( "genome wide combinable domains list is empty" );
        }
        final Set<BinaryDomainCombination> all_binary_combinations = new HashSet<BinaryDomainCombination>();
        for( final GenomeWideCombinableDomains gwcd : getGenomeWideCombinableDomainsList() ) {
            for( final BinaryDomainCombination bc : gwcd.toBinaryDomainCombinations() ) {
                all_binary_combinations.add( bc );
            }
        }
        return all_binary_combinations.size();
    }

    CharacterStateMatrix<BinaryStates> createMatrixOfBinaryDomainCombinationPresenceOrAbsence() {
        return createMatrixOfBinaryDomainCombinationPresenceOrAbsence( getGenomeWideCombinableDomainsList() );
    }

    CharacterStateMatrix<BinaryStates> createMatrixOfDomainPresenceOrAbsence() {
        return createMatrixOfDomainPresenceOrAbsence( getGenomeWideCombinableDomainsList(), getPositiveFilter() );
    }

    CharacterStateMatrix<BinaryStates> createMatrixOfSecondaryFeaturePresenceOrAbsence( final Map<Species, MappingResults> mapping_results_map ) {
        return createMatrixOfSecondaryFeaturePresenceOrAbsence( getGenomeWideCombinableDomainsList(),
                                                                getDomainIdToSecondaryFeaturesMap(),
                                                                mapping_results_map );
    }

    Phylogeny decoratePhylogenyWithDomains( final Phylogeny phylogeny ) {
        for( final PhylogenyNodeIterator it = phylogeny.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            final String node_identifier = node.getNodeName();
            final BinaryCharacters bc = new BinaryCharacters( getUnitsOnNode( node_identifier ),
                                                              getUnitsGainedOnNode( node_identifier ),
                                                              getUnitsLostOnNode( node_identifier ),
                                                              TYPE_FORBINARY_CHARACTERS,
                                                              getSumOfPresentOnNode( node_identifier ),
                                                              getSumOfGainsOnNode( node_identifier ),
                                                              getSumOfLossesOnNode( node_identifier ) );
            node.getNodeData().setBinaryCharacters( bc );
        }
        return phylogeny;
    }

    private void executeDolloParsimony( final boolean on_domain_presence ) {
        reset();
        final DolloParsimony dollo = DolloParsimony.createInstance();
        dollo.setReturnGainLossMatrix( true );
        dollo.setReturnInternalStates( true );
        CharacterStateMatrix<BinaryStates> states = null;
        if ( on_domain_presence ) {
            states = createMatrixOfDomainPresenceOrAbsence();
        }
        else {
            states = createMatrixOfBinaryDomainCombinationPresenceOrAbsence();
        }
        dollo.execute( getPhylogeny(), states );
        setGainLossMatrix( dollo.getGainLossMatrix() );
        setBinaryInternalStatesMatrix( dollo.getInternalStatesMatrix() );
        setCost( dollo.getCost() );
        setTotalGains( dollo.getTotalGains() );
        setTotalLosses( dollo.getTotalLosses() );
        setTotalUnchanged( dollo.getTotalUnchanged() );
    }

    public void executeDolloParsimonyOnBinaryDomainCombintionPresence() {
        executeDolloParsimony( false );
    }

    public void executeDolloParsimonyOnDomainPresence() {
        executeDolloParsimony( true );
    }

    public void executeDolloParsimonyOnDomainPresence( final SortedSet<DomainId> positive_filter ) {
        setPositiveFilter( positive_filter );
        executeDolloParsimony( true );
        setPositiveFilter( null );
    }

    public void executeDolloParsimonyOnSecondaryFeatures( final Map<Species, MappingResults> mapping_results_map ) {
        if ( getDomainIdToSecondaryFeaturesMap() == null ) {
            throw new IllegalStateException( "Domain id to secondary features map has apparently not been set" );
        }
        reset();
        final DolloParsimony dollo = DolloParsimony.createInstance();
        dollo.setReturnGainLossMatrix( true );
        dollo.setReturnInternalStates( true );
        final CharacterStateMatrix<BinaryStates> states = createMatrixOfSecondaryFeaturePresenceOrAbsence( mapping_results_map );
        dollo.execute( getPhylogeny(), states );
        setGainLossMatrix( dollo.getGainLossMatrix() );
        setBinaryInternalStatesMatrix( dollo.getInternalStatesMatrix() );
        setCost( dollo.getCost() );
        setTotalGains( dollo.getTotalGains() );
        setTotalLosses( dollo.getTotalLosses() );
        setTotalUnchanged( dollo.getTotalUnchanged() );
    }

    private void executeFitchParsimony( final boolean on_domain_presence,
                                        final boolean use_last,
                                        final boolean randomize,
                                        final long random_number_seed ) {
        reset();
        if ( use_last ) {
            System.out.println( "   Fitch parsimony: use_last = true" );
        }
        final FitchParsimony<BinaryStates> fitch = new FitchParsimony<BinaryStates>();
        fitch.setRandomize( randomize );
        if ( randomize ) {
            fitch.setRandomNumberSeed( random_number_seed );
        }
        fitch.setUseLast( use_last );
        fitch.setReturnGainLossMatrix( true );
        fitch.setReturnInternalStates( true );
        CharacterStateMatrix<BinaryStates> states = null;
        if ( on_domain_presence ) {
            states = createMatrixOfDomainPresenceOrAbsence( getGenomeWideCombinableDomainsList() );
        }
        else {
            states = createMatrixOfBinaryDomainCombinationPresenceOrAbsence( getGenomeWideCombinableDomainsList() );
        }
        fitch.execute( getPhylogeny(), states );
        setGainLossMatrix( fitch.getGainLossMatrix() );
        setBinaryInternalStatesMatrix( fitch.getInternalStatesMatrix() );
        setCost( fitch.getCost() );
        setTotalGains( fitch.getTotalGains() );
        setTotalLosses( fitch.getTotalLosses() );
        setTotalUnchanged( fitch.getTotalUnchanged() );
    }

    public void executeFitchParsimonyOnBinaryDomainCombintion( final boolean use_last ) {
        executeFitchParsimony( false, use_last, false, 0 );
    }

    public void executeFitchParsimonyOnBinaryDomainCombintion( final long random_number_seed ) {
        executeFitchParsimony( false, false, true, random_number_seed );
    }

    public void executeFitchParsimonyOnDomainPresence( final boolean use_last ) {
        executeFitchParsimony( true, use_last, false, 0 );
    }

    public void executeFitchParsimonyOnDomainPresence( final long random_number_seed ) {
        executeFitchParsimony( true, false, true, random_number_seed );
    }

    public void executeOnGivenBinaryStatesMatrix( final CharacterStateMatrix<BinaryStates> binary_states_matrix,
                                                  final String[] character_labels ) {
        reset();
        if ( binary_states_matrix.getNumberOfCharacters() != character_labels.length ) {
            throw new IllegalArgumentException( "binary states matrix number of characters is not equal to the number of character labels provided" );
        }
        if ( binary_states_matrix.getNumberOfIdentifiers() != getPhylogeny().getNumberOfBranches() ) {
            throw new IllegalArgumentException( "binary states matrix number of identifiers is not equal to the number of tree nodes provided" );
        }
        final CharacterStateMatrix<GainLossStates> gl_matrix = new BasicCharacterStateMatrix<GainLossStates>( binary_states_matrix
                                                                                                                      .getNumberOfIdentifiers(),
                                                                                                              binary_states_matrix
                                                                                                                      .getNumberOfCharacters() );
        int total_gains = 0;
        int total_losses = 0;
        int total_unchanged = 0;
        int i = 0;
        for( final PhylogenyNodeIterator it = getPhylogeny().iteratorPostorder(); it.hasNext(); ) {
            gl_matrix.setIdentifier( i++, it.next().getNodeName() );
        }
        for( int c = 0; c < character_labels.length; ++c ) {
            gl_matrix.setCharacter( c, character_labels[ c ] );
            final PhylogenyNodeIterator it = getPhylogeny().iteratorPostorder();
            while ( it.hasNext() ) {
                final PhylogenyNode node = it.next();
                final String name = node.getNodeName();
                final BinaryStates bin_state = binary_states_matrix.getState( binary_states_matrix
                        .getIdentifierIndex( name ), c );
                final PhylogenyNode parent_node = getPhylogeny().getNode( name ).getParent();
                GainLossStates gl_state = null;
                if ( node.isRoot() ) {
                    ++total_unchanged;
                    if ( bin_state == BinaryStates.ABSENT ) {
                        gl_state = GainLossStates.UNCHANGED_ABSENT;
                    }
                    else {
                        gl_state = GainLossStates.UNCHANGED_PRESENT;
                    }
                }
                else {
                    final BinaryStates parent_bin_state = binary_states_matrix.getState( binary_states_matrix
                            .getIdentifierIndex( parent_node.getNodeName() ), c );
                    if ( bin_state == BinaryStates.ABSENT ) {
                        if ( parent_bin_state == BinaryStates.ABSENT ) {
                            ++total_unchanged;
                            gl_state = GainLossStates.UNCHANGED_ABSENT;
                        }
                        else {
                            ++total_losses;
                            gl_state = GainLossStates.LOSS;
                        }
                    }
                    else {
                        if ( parent_bin_state == BinaryStates.ABSENT ) {
                            ++total_gains;
                            gl_state = GainLossStates.GAIN;
                        }
                        else {
                            ++total_unchanged;
                            gl_state = GainLossStates.UNCHANGED_PRESENT;
                        }
                    }
                }
                gl_matrix.setState( name, c, gl_state );
            }
        }
        setTotalGains( total_gains );
        setTotalLosses( total_losses );
        setTotalUnchanged( total_unchanged );
        setCost( total_gains + total_losses );
        setGainLossMatrix( gl_matrix );
    }

    public int getCost() {
        return _cost;
    }

    private Map<DomainId, Set<String>> getDomainIdToSecondaryFeaturesMap() {
        return _domain_id_to_secondary_features_map;
    }

    public CharacterStateMatrix<Integer> getGainLossCountsMatrix() {
        final CharacterStateMatrix<Integer> matrix = new BasicCharacterStateMatrix<Integer>( getGainLossMatrix()
                .getNumberOfIdentifiers(), 3 );
        for( int i = 0; i < getGainLossMatrix().getNumberOfIdentifiers(); ++i ) {
            matrix.setIdentifier( i, getGainLossMatrix().getIdentifier( i ) );
        }
        matrix.setCharacter( 0, "GAINS" );
        matrix.setCharacter( 1, "LOSSES" );
        matrix.setCharacter( 2, "NET" );
        for( int i = 0; i < getGainLossMatrix().getNumberOfIdentifiers(); ++i ) {
            int gains = 0;
            int losses = 0;
            for( int c = 0; c < getGainLossMatrix().getNumberOfCharacters(); ++c ) {
                final GainLossStates s = getGainLossMatrix().getState( i, c );
                if ( s == GainLossStates.GAIN ) {
                    ++gains;
                }
                else if ( s == GainLossStates.LOSS ) {
                    ++losses;
                }
            }
            matrix.setState( i, 0, gains );
            matrix.setState( i, 1, losses );
            matrix.setState( i, 2, gains - losses );
        }
        return matrix;
    }

    public CharacterStateMatrix<GainLossStates> getGainLossMatrix() {
        return _gain_loss_matrix;
    }

    private List<GenomeWideCombinableDomains> getGenomeWideCombinableDomainsList() {
        return _gwcd_list;
    }

    public CharacterStateMatrix<BinaryStates> getInternalStatesMatrix() {
        return _binary_internal_states_matrix;
    }

    public int getNetGainsOnNode( final String node_identifier ) {
        if ( getGainLossMatrix() == null ) {
            throw new IllegalStateException( "no gain loss matrix has been calculated" );
        }
        int net = 0;
        final int id_index = getGainLossMatrix().getIdentifierIndex( node_identifier );
        for( int c = 0; c < getGainLossMatrix().getNumberOfCharacters(); ++c ) {
            if ( getGainLossMatrix().getState( id_index, c ) == GainLossStates.GAIN ) {
                ++net;
            }
            else if ( getGainLossMatrix().getState( id_index, c ) == GainLossStates.LOSS ) {
                --net;
            }
        }
        return net;
    }

    private Phylogeny getPhylogeny() {
        return _phylogeny;
    }

    private SortedSet<DomainId> getPositiveFilter() {
        return _positive_filter;
    }

    public int getSumOfGainsOnNode( final String node_identifier ) {
        return getStateSumDeltaOnNode( node_identifier, getGainLossMatrix(), GainLossStates.GAIN );
    }

    public int getSumOfLossesOnNode( final String node_identifier ) {
        return getStateSumDeltaOnNode( node_identifier, getGainLossMatrix(), GainLossStates.LOSS );
    }

    public int getSumOfPresentOnNode( final String node_identifier ) {
        return getSumOfGainsOnNode( node_identifier ) + getSumOfUnchangedPresentOnNode( node_identifier );
    }

    int getSumOfUnchangedAbsentOnNode( final String node_identifier ) {
        return getStateSumDeltaOnNode( node_identifier, getGainLossMatrix(), GainLossStates.UNCHANGED_ABSENT );
    }

    int getSumOfUnchangedOnNode( final String node_identifier ) {
        return getSumOfUnchangedPresentOnNode( node_identifier ) + getSumOfUnchangedAbsentOnNode( node_identifier );
    }

    int getSumOfUnchangedPresentOnNode( final String node_identifier ) {
        return getStateSumDeltaOnNode( node_identifier, getGainLossMatrix(), GainLossStates.UNCHANGED_PRESENT );
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

    public SortedSet<String> getUnitsGainedOnNode( final String node_identifier ) {
        return getUnitsDeltaOnNode( node_identifier, getGainLossMatrix(), GainLossStates.GAIN );
    }

    public SortedSet<String> getUnitsLostOnNode( final String node_identifier ) {
        return getUnitsDeltaOnNode( node_identifier, getGainLossMatrix(), GainLossStates.LOSS );
    }

    public SortedSet<String> getUnitsOnNode( final String node_identifier ) {
        final SortedSet<String> present = getUnitsGainedOnNode( node_identifier );
        present.addAll( getUnitsUnchangedPresentOnNode( node_identifier ) );
        return present;
    }

    SortedSet<String> getUnitsUnchangedAbsentOnNode( final String node_identifier ) {
        return getUnitsDeltaOnNode( node_identifier, getGainLossMatrix(), GainLossStates.UNCHANGED_ABSENT );
    }

    SortedSet<String> getUnitsUnchangedPresentOnNode( final String node_identifier ) {
        return getUnitsDeltaOnNode( node_identifier, getGainLossMatrix(), GainLossStates.UNCHANGED_PRESENT );
    }

    private void init() {
        setDomainIdToSecondaryFeaturesMap( null );
        setPositiveFilter( null );
        reset();
    }

    private void reset() {
        setGainLossMatrix( null );
        setBinaryInternalStatesMatrix( null );
        setCost( -1 );
        setTotalGains( -1 );
        setTotalLosses( -1 );
        setTotalUnchanged( -1 );
    }

    private void setBinaryInternalStatesMatrix( final CharacterStateMatrix<BinaryStates> binary_states_matrix ) {
        _binary_internal_states_matrix = binary_states_matrix;
    }

    private void setCost( final int cost ) {
        _cost = cost;
    }

    private void setDomainIdToSecondaryFeaturesMap( final Map<DomainId, Set<String>> domain_id_to_secondary_features_map ) {
        _domain_id_to_secondary_features_map = domain_id_to_secondary_features_map;
    }

    private void setGainLossMatrix( final CharacterStateMatrix<GainLossStates> gain_loss_matrix ) {
        _gain_loss_matrix = gain_loss_matrix;
    }

    private void setPositiveFilter( final SortedSet<DomainId> positive_filter ) {
        _positive_filter = positive_filter;
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

    public static DomainParsimonyCalculator createInstance( final Phylogeny phylogeny ) {
        return new DomainParsimonyCalculator( phylogeny );
    }

    public static DomainParsimonyCalculator createInstance( final Phylogeny phylogeny,
                                                            final List<GenomeWideCombinableDomains> gwcd_list ) {
        if ( phylogeny.getNumberOfExternalNodes() != gwcd_list.size() ) {
            throw new IllegalArgumentException( "number of external nodes [" + phylogeny.getNumberOfExternalNodes()
                    + "] does not equal size of genome wide combinable domains list [" + gwcd_list.size() + "]" );
        }
        return new DomainParsimonyCalculator( phylogeny, gwcd_list );
    }

    public static DomainParsimonyCalculator createInstance( final Phylogeny phylogeny,
                                                            final List<GenomeWideCombinableDomains> gwcd_list,
                                                            final Map<DomainId, Set<String>> domain_id_to_secondary_features_map ) {
        if ( phylogeny.getNumberOfExternalNodes() != gwcd_list.size() ) {
            throw new IllegalArgumentException( "size of external nodes does not equal size of genome wide combinable domains list" );
        }
        return new DomainParsimonyCalculator( phylogeny, gwcd_list, domain_id_to_secondary_features_map );
    }

    public static CharacterStateMatrix<BinaryStates> createMatrixOfBinaryDomainCombinationPresenceOrAbsence( final List<GenomeWideCombinableDomains> gwcd_list ) {
        if ( gwcd_list.isEmpty() ) {
            throw new IllegalArgumentException( "genome wide combinable domains list is empty" );
        }
        final int number_of_identifiers = gwcd_list.size();
        final SortedSet<BinaryDomainCombination> all_binary_combinations = new TreeSet<BinaryDomainCombination>();
        final Set<BinaryDomainCombination>[] binary_combinations_per_genome = new HashSet[ number_of_identifiers ];
        int identifier_index = 0;
        for( final GenomeWideCombinableDomains gwcd : gwcd_list ) {
            binary_combinations_per_genome[ identifier_index ] = new HashSet<BinaryDomainCombination>();
            for( final BinaryDomainCombination bc : gwcd.toBinaryDomainCombinations() ) {
                all_binary_combinations.add( bc );
                binary_combinations_per_genome[ identifier_index ].add( bc );
            }
            ++identifier_index;
        }
        final int number_of_characters = all_binary_combinations.size();
        final CharacterStateMatrix<CharacterStateMatrix.BinaryStates> matrix = new BasicCharacterStateMatrix<CharacterStateMatrix.BinaryStates>( number_of_identifiers,
                                                                                                                                                 number_of_characters );
        int character_index = 0;
        for( final BinaryDomainCombination bc : all_binary_combinations ) {
            matrix.setCharacter( character_index++, bc.toString() );
        }
        identifier_index = 0;
        final Set<String> all_identifiers = new HashSet<String>();
        for( final GenomeWideCombinableDomains gwcd : gwcd_list ) {
            final String species_id = gwcd.getSpecies().getSpeciesId();
            if ( all_identifiers.contains( species_id ) ) {
                throw new AssertionError( "species [" + species_id + "] is not unique" );
            }
            all_identifiers.add( species_id );
            matrix.setIdentifier( identifier_index, species_id );
            for( int ci = 0; ci < matrix.getNumberOfCharacters(); ++ci ) {
                BinaryDomainCombination bc = null;
                if ( gwcd.getDomainCombinationType() == DomainCombinationType.DIRECTED_ADJACTANT ) {
                    bc = AdjactantDirectedBinaryDomainCombination.createInstance( matrix.getCharacter( ci ) );
                }
                else if ( gwcd.getDomainCombinationType() == DomainCombinationType.DIRECTED ) {
                    bc = DirectedBinaryDomainCombination.createInstance( matrix.getCharacter( ci ) );
                }
                else {
                    bc = BasicBinaryDomainCombination.createInstance( matrix.getCharacter( ci ) );
                }
                if ( binary_combinations_per_genome[ identifier_index ].contains( bc ) ) {
                    matrix.setState( identifier_index, ci, CharacterStateMatrix.BinaryStates.PRESENT );
                }
                else {
                    matrix.setState( identifier_index, ci, CharacterStateMatrix.BinaryStates.ABSENT );
                }
            }
            ++identifier_index;
        }
        return matrix;
    }

    static CharacterStateMatrix<BinaryStates> createMatrixOfDomainPresenceOrAbsence( final List<GenomeWideCombinableDomains> gwcd_list ) {
        return createMatrixOfDomainPresenceOrAbsence( gwcd_list, null );
    }

    public static CharacterStateMatrix<BinaryStates> createMatrixOfDomainPresenceOrAbsence( final List<GenomeWideCombinableDomains> gwcd_list,
                                                                                            final SortedSet<DomainId> positive_filter ) {
        if ( gwcd_list.isEmpty() ) {
            throw new IllegalArgumentException( "genome wide combinable domains list is empty" );
        }
        if ( ( positive_filter != null ) && ( positive_filter.size() < 1 ) ) {
            throw new IllegalArgumentException( "positive filter is empty" );
        }
        final int number_of_identifiers = gwcd_list.size();
        final SortedSet<DomainId> all_domain_ids = new TreeSet<DomainId>();
        for( final GenomeWideCombinableDomains gwcd : gwcd_list ) {
            for( final DomainId domain : gwcd.getAllDomainIds() ) {
                all_domain_ids.add( domain );
            }
        }
        int number_of_characters = all_domain_ids.size();
        if ( positive_filter != null ) {
            number_of_characters = positive_filter.size();
        }
        final CharacterStateMatrix<CharacterStateMatrix.BinaryStates> matrix = new BasicCharacterStateMatrix<CharacterStateMatrix.BinaryStates>( number_of_identifiers,
                                                                                                                                                 number_of_characters );
        int character_index = 0;
        for( final DomainId id : all_domain_ids ) {
            if ( positive_filter == null ) {
                matrix.setCharacter( character_index++, id.getId() );
            }
            else {
                if ( positive_filter.contains( id ) ) {
                    matrix.setCharacter( character_index++, id.getId() );
                }
            }
        }
        int identifier_index = 0;
        final Set<String> all_identifiers = new HashSet<String>();
        for( final GenomeWideCombinableDomains gwcd : gwcd_list ) {
            final String species_id = gwcd.getSpecies().getSpeciesId();
            if ( all_identifiers.contains( species_id ) ) {
                throw new IllegalArgumentException( "species [" + species_id + "] is not unique" );
            }
            all_identifiers.add( species_id );
            matrix.setIdentifier( identifier_index, species_id );
            for( int ci = 0; ci < matrix.getNumberOfCharacters(); ++ci ) {
                if ( ForesterUtil.isEmpty( matrix.getCharacter( ci ) ) ) {
                    throw new IllegalArgumentException( "problem with character #" + ci
                            + ", possibly using domain(s) in positive filter not present in the genomes analyzed" );
                }
                final DomainId id = new DomainId( matrix.getCharacter( ci ) );
                if ( gwcd.contains( id ) ) {
                    matrix.setState( identifier_index, ci, CharacterStateMatrix.BinaryStates.PRESENT );
                }
                else {
                    matrix.setState( identifier_index, ci, CharacterStateMatrix.BinaryStates.ABSENT );
                }
            }
            ++identifier_index;
        }
        return matrix;
    }

    /**
     * For folds instead of Pfam-domains, for example
     * 
     * 
     * @param gwcd_list
     * @return
     */
    static CharacterStateMatrix<BinaryStates> createMatrixOfSecondaryFeaturePresenceOrAbsence( final List<GenomeWideCombinableDomains> gwcd_list,
                                                                                               final Map<DomainId, Set<String>> domain_id_to_second_features_map,
                                                                                               final Map<Species, MappingResults> mapping_results_map ) {
        if ( gwcd_list.isEmpty() ) {
            throw new IllegalArgumentException( "genome wide combinable domains list is empty" );
        }
        if ( ( domain_id_to_second_features_map == null ) || domain_id_to_second_features_map.isEmpty() ) {
            throw new IllegalArgumentException( "domain id to secondary features map is null or empty" );
        }
        final int number_of_identifiers = gwcd_list.size();
        final SortedSet<String> all_secondary_features = new TreeSet<String>();
        for( final GenomeWideCombinableDomains gwcd : gwcd_list ) {
            int mapped = 0;
            int not_mapped = 0;
            for( final DomainId domain : gwcd.getAllDomainIds() ) {
                if ( domain_id_to_second_features_map.containsKey( domain ) ) {
                    all_secondary_features.addAll( domain_id_to_second_features_map.get( domain ) );
                    mapped++;
                }
                else {
                    not_mapped++;
                }
            }
            if ( mapping_results_map != null ) {
                final MappingResults mr = new MappingResults();
                mr.setDescription( gwcd.getSpecies().getSpeciesId() );
                mr.setSumOfSuccesses( mapped );
                mr.setSumOfFailures( not_mapped );
                mapping_results_map.put( gwcd.getSpecies(), mr );
            }
        }
        final int number_of_characters = all_secondary_features.size();
        final CharacterStateMatrix<CharacterStateMatrix.BinaryStates> matrix = new BasicCharacterStateMatrix<CharacterStateMatrix.BinaryStates>( number_of_identifiers,
                                                                                                                                                 number_of_characters );
        int character_index = 0;
        for( final String second_id : all_secondary_features ) {
            matrix.setCharacter( character_index++, second_id );
        }
        int identifier_index = 0;
        final Set<String> all_identifiers = new HashSet<String>();
        for( final GenomeWideCombinableDomains gwcd : gwcd_list ) {
            final String species_id = gwcd.getSpecies().getSpeciesId();
            if ( all_identifiers.contains( species_id ) ) {
                throw new IllegalArgumentException( "species [" + species_id + "] is not unique" );
            }
            all_identifiers.add( species_id );
            matrix.setIdentifier( identifier_index, species_id );
            final Set<String> all_second_per_gwcd = new HashSet<String>();
            for( final DomainId domain : gwcd.getAllDomainIds() ) {
                if ( domain_id_to_second_features_map.containsKey( domain ) ) {
                    all_second_per_gwcd.addAll( domain_id_to_second_features_map.get( domain ) );
                }
            }
            for( int ci = 0; ci < matrix.getNumberOfCharacters(); ++ci ) {
                if ( all_second_per_gwcd.contains( matrix.getCharacter( ci ) ) ) {
                    matrix.setState( identifier_index, ci, CharacterStateMatrix.BinaryStates.PRESENT );
                }
                else {
                    matrix.setState( identifier_index, ci, CharacterStateMatrix.BinaryStates.ABSENT );
                }
            }
            ++identifier_index;
        }
        return matrix;
    }

    private static int getStateSumDeltaOnNode( final String node_identifier,
                                               final CharacterStateMatrix<GainLossStates> gain_loss_matrix,
                                               final GainLossStates state ) {
        if ( gain_loss_matrix == null ) {
            throw new IllegalStateException( "no gain loss matrix has been calculated" );
        }
        if ( ForesterUtil.isEmpty( node_identifier ) ) {
            throw new IllegalArgumentException( "node identifier must not be empty" );
        }
        if ( gain_loss_matrix.isEmpty() ) {
            throw new IllegalStateException( "gain loss matrix is empty" );
        }
        int sum = 0;
        final int id_index = gain_loss_matrix.getIdentifierIndex( node_identifier );
        for( int c = 0; c < gain_loss_matrix.getNumberOfCharacters(); ++c ) {
            if ( gain_loss_matrix.getState( id_index, c ) == state ) {
                ++sum;
            }
        }
        return sum;
    }

    private static SortedSet<String> getUnitsDeltaOnNode( final String node_identifier,
                                                          final CharacterStateMatrix<GainLossStates> gain_loss_matrix,
                                                          final GainLossStates state ) {
        if ( gain_loss_matrix == null ) {
            throw new IllegalStateException( "no gain loss matrix has been calculated" );
        }
        if ( ForesterUtil.isEmpty( node_identifier ) ) {
            throw new IllegalArgumentException( "node identifier must not be empty" );
        }
        if ( gain_loss_matrix.isEmpty() ) {
            throw new IllegalStateException( "gain loss matrix is empty" );
        }
        final SortedSet<String> d = new TreeSet<String>();
        final int id_index = gain_loss_matrix.getIdentifierIndex( node_identifier );
        for( int c = 0; c < gain_loss_matrix.getNumberOfCharacters(); ++c ) {
            if ( gain_loss_matrix.getState( id_index, c ) == state ) {
                if ( d.contains( gain_loss_matrix.getCharacter( c ) ) ) {
                    throw new AssertionError( "this should not have happended: character ["
                            + gain_loss_matrix.getCharacter( c ) + "] already in set" );
                }
                d.add( gain_loss_matrix.getCharacter( c ) );
            }
        }
        return d;
    }
}
