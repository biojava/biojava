// $Id: DomainLengthsTable.java,v 1.9 2009/11/11 01:36:50 cmzmasek Exp $
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

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class DomainLengthsTable {

    private final static DecimalFormat       DF = new DecimalFormat( "#.0" );
    final SortedMap<DomainId, DomainLengths> _domain_lengths;
    final List<Species>                      _species;

    public DomainLengthsTable() {
        _domain_lengths = new TreeMap<DomainId, DomainLengths>();
        _species = new ArrayList<Species>();
    }

    private void addDomainLengths( final DomainLengths domain_lengths ) {
        if ( getDomainLengths().containsKey( domain_lengths.getDomainId() ) ) {
            throw new IllegalArgumentException( "domain lengths for [" + domain_lengths.getDomainId()
                    + "] already added" );
        }
        getDomainLengths().put( domain_lengths.getDomainId(), domain_lengths );
    }

    private void addLength( final DomainId domain_id, final Species species, final int domain_length ) {
        if ( !getDomainLengths().containsKey( domain_id ) ) {
            addDomainLengths( new DomainLengths( domain_id ) );
        }
        getDomainLengths().get( domain_id ).addLength( species, domain_length );
    }

    public void addLengths( final List<Protein> protein_list ) {
        for( final Protein protein : protein_list ) {
            final Species species = protein.getSpecies();
            if ( !_species.contains( species ) ) {
                _species.add( species );
            }
            for( final Domain domain : protein.getProteinDomains() ) {
                addLength( domain.getDomainId(), species, ( domain.getTo() - domain.getFrom() ) + 1 );
            }
        }
    }

    public DescriptiveStatistics calculateMeanBasedStatisticsForAllSpecies() {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( final Species species : getSpecies() ) {
            final DescriptiveStatistics stats_per_species = calculateMeanBasedStatisticsForSpecies( species );
            stats.addValue( stats_per_species.arithmeticMean() );
        }
        return stats;
    }

    public DescriptiveStatistics calculateMeanBasedStatisticsForDomain( final DomainId domain_id ) {
        return getDomainLengths( domain_id ).calculateMeanBasedStatistics();
    }

    public DescriptiveStatistics calculateMeanBasedStatisticsForSpecies( final Species species ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( final DomainLengths l : getDomainLengths().values() ) {
            if ( l.isHasLengthStatistic( species ) ) {
                stats.addValue( l.getLengthStatistic( species ).arithmeticMean() );
            }
        }
        return stats;
    }

    public StringBuilder createMeanBasedStatisticsPerSpeciesTable() {
        final StringBuilder sb = new StringBuilder();
        sb.append( "SPECIES" );
        sb.append( "\t" );
        sb.append( "MEAN" );
        sb.append( "\t" );
        sb.append( "SD" );
        sb.append( "\t" );
        sb.append( "MIN" );
        sb.append( "\t" );
        sb.append( "MAX" );
        sb.append( "\t" );
        sb.append( "MEDIAN" );
        sb.append( ForesterUtil.LINE_SEPARATOR );
        for( final Species species : getSpecies() ) {
            final DescriptiveStatistics stats = calculateMeanBasedStatisticsForSpecies( species );
            sb.append( species );
            sb.append( "\t" );
            sb.append( DF.format( stats.arithmeticMean() ) );
            sb.append( "\t" );
            try {
                sb.append( DF.format( stats.sampleStandardDeviation() ) );
            }
            catch ( final ArithmeticException e ) {
                sb.append( "" );
            }
            sb.append( "\t" );
            sb.append( DF.format( stats.getMin() ) );
            sb.append( "\t" );
            sb.append( DF.format( stats.getMax() ) );
            sb.append( "\t" );
            try {
                sb.append( DF.format( stats.median() ) );
            }
            catch ( final ArithmeticException e ) {
                sb.append( "" );
            }
            sb.append( ForesterUtil.LINE_SEPARATOR );
        }
        return sb;
    }

    private SortedMap<DomainId, DomainLengths> getDomainLengths() {
        return _domain_lengths;
    }

    public DomainLengths getDomainLengths( final DomainId domain_id ) {
        return getDomainLengths().get( domain_id );
    }

    public List<DomainLengths> getDomainLengthsList() {
        final List<DomainLengths> list = new ArrayList<DomainLengths>();
        for( final DomainLengths l : getDomainLengths().values() ) {
            list.add( l );
        }
        return list;
    }

    public DescriptiveStatistics getLengthStatistic( final DomainId domain_id, final Species species ) {
        return getDomainLengths( domain_id ).getLengthStatistic( species );
    }

    public List<Species> getSpecies() {
        return _species;
    }
}
