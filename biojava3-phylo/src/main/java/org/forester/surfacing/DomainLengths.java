// $Id: DomainLengths.java,v 1.6 2009/10/26 23:29:40 cmzmasek Exp $
//
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2010 Christian M. Zmasek
// Copyright (C) 2008-2010 Burnham Institute for Medical Research
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

import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;

public class DomainLengths {

    final DomainId                                  _domain_id;
    final SortedMap<Species, DescriptiveStatistics> _length_statistics;

    public DomainLengths( final DomainId domain_id ) {
        _domain_id = domain_id;
        _length_statistics = new TreeMap<Species, DescriptiveStatistics>();
    }

    public void addLength( final Species species, final int domain_length ) {
        if ( !getLengthStatistics().containsKey( species ) ) {
            addLengthStatistics( species, new BasicDescriptiveStatistics() );
        }
        getLengthStatistic( species ).addValue( domain_length );
    }

    private void addLengthStatistics( final Species species, final DescriptiveStatistics length_statistic ) {
        if ( getLengthStatistics().containsKey( species ) ) {
            throw new IllegalArgumentException( "length statistics for [" + species.getSpeciesId() + "] already added" );
        }
        getLengthStatistics().put( species, length_statistic );
    }

    /**
     * Returns descriptive statistics based on the arithmetic means
     * for each species.  
     * 
     * 
     * @return
     */
    public DescriptiveStatistics calculateMeanBasedStatistics() {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( final DescriptiveStatistics s : getLengthStatisticsList() ) {
            stats.addValue( s.arithmeticMean() );
        }
        return stats;
    }

    /**
     * 
     * Note. This is not technically a Z-score since the distribution
     * of means is unknown (and not normal).
     * 
     * @param species
     * @return
     */
    public double calculateZScoreForSpecies( final Species species ) {
        final double species_mean = getLengthStatistic( species ).arithmeticMean();
        final DescriptiveStatistics domain_stats = calculateMeanBasedStatistics();
        final double population_sd = domain_stats.sampleStandardDeviation();
        final double population_mean = domain_stats.arithmeticMean();
        return ( species_mean - population_mean ) / population_sd;
    }

    public DomainId getDomainId() {
        return _domain_id;
    }

    public DescriptiveStatistics getLengthStatistic( final Species species ) {
        return getLengthStatistics().get( species );
    }

    private SortedMap<Species, DescriptiveStatistics> getLengthStatistics() {
        return _length_statistics;
    }

    public List<DescriptiveStatistics> getLengthStatisticsList() {
        final List<DescriptiveStatistics> list = new ArrayList<DescriptiveStatistics>();
        for( final DescriptiveStatistics stats : _length_statistics.values() ) {
            list.add( stats );
        }
        return list;
    }

    public List<Species> getMeanBasedOutlierSpecies( final double z_score_limit ) {
        final List<Species> species = new ArrayList<Species>();
        if ( getSpeciesList().size() > 1 ) {
            for( final Species s : getSpeciesList() ) {
                final double z = calculateZScoreForSpecies( s );
                if ( z_score_limit < 0 ) {
                    if ( z <= z_score_limit ) {
                        species.add( s );
                    }
                }
                else if ( z_score_limit > 0 ) {
                    if ( z >= z_score_limit ) {
                        species.add( s );
                    }
                }
            }
        }
        return species;
    }

    public List<Species> getSpeciesList() {
        final List<Species> list = new ArrayList<Species>();
        for( final Species s : _length_statistics.keySet() ) {
            list.add( s );
        }
        return list;
    }

    public boolean isHasLengthStatistic( final Species species ) {
        return getLengthStatistics().containsKey( species );
    }
}
