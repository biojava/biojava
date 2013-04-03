// $Id: CombinableDomainsCollectionSimilarityCalculator.java,v 1.2 2007/11/01
// 19:38:35 cmzmasek Exp $
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
import java.util.Set;

public class DomainArchitectureBasedGenomeSimilarityCalculator {

    public static final double                MAX_SIMILARITY_SCORE = 1.0;
    public static final double                MIN_SIMILARITY_SCORE = 0.0;
    final private GenomeWideCombinableDomains _combinable_domains_genome_0;
    final private GenomeWideCombinableDomains _combinable_domains_genome_1;
    private Set<DomainId>                     _domain_ids_to_ignore;
    private boolean                           _allow_domains_to_be_ignored;
    private Set<DomainId>                     _all_domains;
    private Set<DomainId>                     _shared_domains;
    private Set<DomainId>                     _domains_specific_to_0;
    private Set<DomainId>                     _domains_specific_to_1;
    private Set<BinaryDomainCombination>      _all_binary_domain_combinations;
    private Set<BinaryDomainCombination>      _shared_binary_domain_combinations;
    private Set<BinaryDomainCombination>      _binary_domain_combinations_specific_to_0;
    private Set<BinaryDomainCombination>      _binary_domain_combinations_specific_to_1;

    public DomainArchitectureBasedGenomeSimilarityCalculator( final GenomeWideCombinableDomains combinable_domains_genome_0,
                                                              final GenomeWideCombinableDomains combinable_domains_genome_1 ) {
        if ( ( combinable_domains_genome_0 == null ) || ( combinable_domains_genome_0.getSize() < 1 )
                || ( combinable_domains_genome_1 == null ) || ( combinable_domains_genome_1.getSize() < 1 ) ) {
            throw new IllegalArgumentException( "attempt to compare null or empty combinable domains collection" );
        }
        if ( combinable_domains_genome_0.getSpecies().equals( combinable_domains_genome_1.getSpecies() ) ) {
            throw new IllegalArgumentException( "attempt to compare combinable domains collection from the same species" );
        }
        _combinable_domains_genome_0 = combinable_domains_genome_0;
        _combinable_domains_genome_1 = combinable_domains_genome_1;
        init();
        forceRecalculation();
    }

    public void addDomainIdToIgnore( final DomainId domain_id_to_ignore ) {
        forceRecalculation();
        getDomainIdsToIgnore().add( domain_id_to_ignore );
    }

    /**
     * This returns a score between 0.0 (no binary domain combination in common) 
     * and 1.0 (all binary domain combinations in common) measuring the similarity between two
     * genomes based on the number of shared binary domain combinations:
     *   
     * t: sum of (distinct) binary domain combinations
     * s: sum of shared (distinct) binary domain combinations
     *
     * 1 - ( ( t - s ) / t )
     *  
     * @return shared binary domain combinations based similarity score 
     */
    public double calculateSharedBinaryDomainCombinationBasedGenomeSimilarityScore() {
        final double t = getAllBinaryDomainCombinations().size();
        final double s = getSharedBinaryDomainCombinations().size();
        if ( t == 0.0 ) {
            return MIN_SIMILARITY_SCORE;
        }
        return ( MAX_SIMILARITY_SCORE - ( ( t - s ) / t ) );
    }

    /**
     * This returns a score between 0.0 (no domains in common) 
     * and 1.0 (all domains in common) measuring the similarity between two
     * genomes based on the number of shared domains:
     * 
     * t: sum of (distinct) domains
     * s: sum of shared (distinct) domains
     *
     * 1 - ( ( t - s ) / t )
     * 
     * @return shared domains based similarity score 
     */
    public double calculateSharedDomainsBasedGenomeSimilarityScore() {
        final double t = getAllDomains().size();
        final double s = getSharedDomains().size();
        if ( t == 0.0 ) {
            return MIN_SIMILARITY_SCORE;
        }
        return ( MAX_SIMILARITY_SCORE - ( ( t - s ) / t ) );
    }

    public void deleteAllDomainIdsToIgnore() {
        forceRecalculation();
        setDomainIdsToIgnore( new HashSet<DomainId>() );
    }

    private void forceRecalculation() {
        _all_domains = null;
        _shared_domains = null;
        _domains_specific_to_0 = null;
        _domains_specific_to_1 = null;
        _all_binary_domain_combinations = null;
        _shared_binary_domain_combinations = null;
        _binary_domain_combinations_specific_to_0 = null;
        _binary_domain_combinations_specific_to_1 = null;
    }

    /**
     * Does not return binary combinations which contain one or two domains
     * to be ignored -- if ignoring is allowed.
     * 
     * @return SortedSet<BinaryDomainCombination>
     */
    public Set<BinaryDomainCombination> getAllBinaryDomainCombinations() {
        if ( _all_binary_domain_combinations == null ) {
            final Set<BinaryDomainCombination> all = new HashSet<BinaryDomainCombination>();
            all.addAll( getCombinableDomainsGenome0().toBinaryDomainCombinations() );
            all.addAll( getCombinableDomainsGenome1().toBinaryDomainCombinations() );
            if ( isAllowDomainsToBeIgnored() && !getDomainIdsToIgnore().isEmpty() ) {
                _all_binary_domain_combinations = pruneBinaryCombinations( all );
            }
            else {
                _all_binary_domain_combinations = all;
            }
        }
        return _all_binary_domain_combinations;
    }

    /**
     * Does not return domains which are to be
     * ignored -- if ignoring is allowed.
     * 
     * 
     * @return
     */
    public Set<DomainId> getAllDomains() {
        if ( _all_domains == null ) {
            final Set<DomainId> all = new HashSet<DomainId>();
            all.addAll( getCombinableDomainsGenome0().getAllDomainIds() );
            all.addAll( getCombinableDomainsGenome1().getAllDomainIds() );
            if ( isAllowDomainsToBeIgnored() && !getDomainIdsToIgnore().isEmpty() ) {
                _all_domains = pruneDomains( all );
            }
            else {
                _all_domains = all;
            }
        }
        return _all_domains;
    }

    private Set<BinaryDomainCombination> getBinaryDomainCombinationsSpecificToGenome( final boolean specific_to_genome_0 ) {
        final Set<BinaryDomainCombination> specific = new HashSet<BinaryDomainCombination>();
        final Set<BinaryDomainCombination> bc0 = getCombinableDomainsGenome0().toBinaryDomainCombinations();
        final Set<BinaryDomainCombination> bc1 = getCombinableDomainsGenome1().toBinaryDomainCombinations();
        if ( specific_to_genome_0 ) {
            for( final BinaryDomainCombination binary_domain_combination0 : bc0 ) {
                if ( !bc1.contains( binary_domain_combination0 ) ) {
                    specific.add( binary_domain_combination0 );
                }
            }
        }
        else {
            for( final BinaryDomainCombination binary_domain_combination1 : bc1 ) {
                if ( !bc0.contains( binary_domain_combination1 ) ) {
                    specific.add( binary_domain_combination1 );
                }
            }
        }
        if ( isAllowDomainsToBeIgnored() && !getDomainIdsToIgnore().isEmpty() ) {
            return pruneBinaryCombinations( specific );
        }
        return specific;
    }

    public Set<BinaryDomainCombination> getBinaryDomainCombinationsSpecificToGenome0() {
        if ( _binary_domain_combinations_specific_to_0 == null ) {
            _binary_domain_combinations_specific_to_0 = getBinaryDomainCombinationsSpecificToGenome( true );
        }
        return _binary_domain_combinations_specific_to_0;
    }

    public Set<BinaryDomainCombination> getBinaryDomainCombinationsSpecificToGenome1() {
        if ( _binary_domain_combinations_specific_to_1 == null ) {
            _binary_domain_combinations_specific_to_1 = getBinaryDomainCombinationsSpecificToGenome( false );
        }
        return _binary_domain_combinations_specific_to_1;
    }

    private GenomeWideCombinableDomains getCombinableDomainsGenome0() {
        return _combinable_domains_genome_0;
    }

    private GenomeWideCombinableDomains getCombinableDomainsGenome1() {
        return _combinable_domains_genome_1;
    }

    private Set<DomainId> getDomainIdsToIgnore() {
        return _domain_ids_to_ignore;
    }

    private Set<DomainId> getDomainsSpecificToGenome( final boolean specific_to_genome_0 ) {
        final Set<DomainId> specific = new HashSet<DomainId>();
        final Set<DomainId> d0 = getCombinableDomainsGenome0().getAllDomainIds();
        final Set<DomainId> d1 = getCombinableDomainsGenome1().getAllDomainIds();
        if ( specific_to_genome_0 ) {
            for( final DomainId domain0 : d0 ) {
                if ( !d1.contains( domain0 ) ) {
                    specific.add( domain0 );
                }
            }
        }
        else {
            for( final DomainId domain1 : d1 ) {
                if ( !d0.contains( domain1 ) ) {
                    specific.add( domain1 );
                }
            }
        }
        if ( isAllowDomainsToBeIgnored() && !getDomainIdsToIgnore().isEmpty() ) {
            return pruneDomains( specific );
        }
        return specific;
    }

    public Set<DomainId> getDomainsSpecificToGenome0() {
        if ( _domains_specific_to_0 == null ) {
            _domains_specific_to_0 = getDomainsSpecificToGenome( true );
        }
        return _domains_specific_to_0;
    }

    public Set<DomainId> getDomainsSpecificToGenome1() {
        if ( _domains_specific_to_1 == null ) {
            _domains_specific_to_1 = getDomainsSpecificToGenome( false );
        }
        return _domains_specific_to_1;
    }

    public Set<BinaryDomainCombination> getSharedBinaryDomainCombinations() {
        if ( _shared_binary_domain_combinations == null ) {
            final Set<BinaryDomainCombination> shared = new HashSet<BinaryDomainCombination>();
            final Set<BinaryDomainCombination> bc0 = getCombinableDomainsGenome0().toBinaryDomainCombinations();
            final Set<BinaryDomainCombination> bc1 = getCombinableDomainsGenome1().toBinaryDomainCombinations();
            for( final BinaryDomainCombination binary_domain_combination0 : bc0 ) {
                if ( bc1.contains( binary_domain_combination0 ) ) {
                    shared.add( binary_domain_combination0 );
                }
            }
            _shared_binary_domain_combinations = shared;
            if ( isAllowDomainsToBeIgnored() && !getDomainIdsToIgnore().isEmpty() ) {
                _shared_binary_domain_combinations = pruneBinaryCombinations( shared );
            }
        }
        return _shared_binary_domain_combinations;
    }

    public Set<DomainId> getSharedDomains() {
        if ( _shared_domains == null ) {
            final Set<DomainId> shared = new HashSet<DomainId>();
            final Set<DomainId> d0 = getCombinableDomainsGenome0().getAllDomainIds();
            final Set<DomainId> d1 = getCombinableDomainsGenome1().getAllDomainIds();
            for( final DomainId domain0 : d0 ) {
                if ( d1.contains( domain0 ) ) {
                    shared.add( domain0 );
                }
            }
            _shared_domains = shared;
            if ( isAllowDomainsToBeIgnored() && !getDomainIdsToIgnore().isEmpty() ) {
                _shared_domains = pruneDomains( shared );
            }
        }
        return _shared_domains;
    }

    private void init() {
        deleteAllDomainIdsToIgnore();
        setAllowDomainsToBeIgnored( false );
    }

    private boolean isAllowDomainsToBeIgnored() {
        return _allow_domains_to_be_ignored;
    }

    private Set<BinaryDomainCombination> pruneBinaryCombinations( final Set<BinaryDomainCombination> all ) {
        final Set<BinaryDomainCombination> pruned = new HashSet<BinaryDomainCombination>();
        for( final BinaryDomainCombination bc : all ) {
            if ( ( !getDomainIdsToIgnore().contains( bc.getId0() ) )
                    && ( !getDomainIdsToIgnore().contains( bc.getId1() ) ) ) {
                pruned.add( bc );
            }
        }
        return pruned;
    }

    private Set<DomainId> pruneDomains( final Set<DomainId> all ) {
        final Set<DomainId> pruned = new HashSet<DomainId>();
        for( final DomainId d : all ) {
            if ( !getDomainIdsToIgnore().contains( d ) ) {
                pruned.add( d );
            }
        }
        return pruned;
    }

    public void setAllowDomainsToBeIgnored( final boolean allow_domains_to_be_ignored ) {
        forceRecalculation();
        _allow_domains_to_be_ignored = allow_domains_to_be_ignored;
    }

    void setDomainIdsToIgnore( final Set<DomainId> domain_ids_to_ignore ) {
        forceRecalculation();
        _domain_ids_to_ignore = domain_ids_to_ignore;
    }
}