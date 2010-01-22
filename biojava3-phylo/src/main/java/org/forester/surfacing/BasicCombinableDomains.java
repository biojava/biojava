// $Id: BasicCombinableDomains.java,v 1.14 2009/10/26 23:29:40 cmzmasek Exp $
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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.util.DescriptiveStatistics;

public class BasicCombinableDomains implements CombinableDomains {

    final private DomainId                   _key_domain;
    private int                              _key_domain_count;
    private int                              _key_domain_proteins_count;
    final private Species                    _species;
    final private TreeMap<DomainId, Integer> _combining_domains;
    private DescriptiveStatistics            _key_domain_confidence_statistics;

    public BasicCombinableDomains( final DomainId key_domain, final Species species ) {
        _key_domain = key_domain;
        _species = species;
        _combining_domains = new TreeMap<DomainId, Integer>();
        init();
    }

    public void addCombinableDomain( final DomainId protein_domain ) {
        if ( getCombiningDomains().containsKey( protein_domain ) ) {
            getCombiningDomains().put( protein_domain, getCombiningDomains().get( protein_domain ) + 1 );
        }
        else {
            getCombiningDomains().put( protein_domain, 1 );
        }
    }

    public List<DomainId> getAllDomains() {
        final List<DomainId> domains = getCombinableDomains();
        if ( !domains.contains( getKeyDomain() ) ) {
            domains.add( getKeyDomain() );
        }
        return domains;
    }

    public List<DomainId> getCombinableDomains() {
        final List<DomainId> domains = new ArrayList<DomainId>( getNumberOfCombinableDomains() );
        for( final DomainId domain : getCombiningDomains().keySet() ) {
            domains.add( domain );
        }
        return domains;
    }

    public SortedMap<DomainId, Integer> getCombinableDomainsIds() {
        final SortedMap<DomainId, Integer> ids = new TreeMap<DomainId, Integer>();
        for( final DomainId domain : getCombiningDomains().keySet() ) {
            final DomainId pd = domain;
            ids.put( pd, getCombiningDomains().get( pd ) );
        }
        return ids;
    }

    public StringBuilder getCombiningDomainIdsAsStringBuilder() {
        final StringBuilder sb = new StringBuilder();
        for( final Iterator<DomainId> iter = getCombiningDomains().keySet().iterator(); iter.hasNext(); ) {
            final DomainId key = iter.next();
            sb.append( key.toString() );
            sb.append( " [" );
            final int count = getCombiningDomains().get( key );
            sb.append( count );
            sb.append( "]" );
            if ( iter.hasNext() ) {
                sb.append( ", " );
            }
        }
        return sb;
    }

    protected TreeMap<DomainId, Integer> getCombiningDomains() {
        return _combining_domains;
    }

    public DomainId getKeyDomain() {
        return _key_domain;
    }

    public DescriptiveStatistics getKeyDomainConfidenceDescriptiveStatistics() {
        return _key_domain_confidence_statistics;
    }

    public int getKeyDomainCount() {
        return _key_domain_count;
    }

    public int getKeyDomainProteinsCount() {
        return _key_domain_proteins_count;
    }

    public int getNumberOfCombinableDomains() {
        return _combining_domains.size();
    }

    public int getNumberOfProteinsExhibitingCombination( final DomainId protein_domain ) {
        if ( getCombiningDomains().containsKey( protein_domain ) ) {
            return getCombiningDomains().get( protein_domain );
        }
        else {
            return 0;
        }
    }

    public Species getSpecies() {
        return _species;
    }

    private void init() {
        _key_domain_count = 0;
        _key_domain_proteins_count = 0;
        _key_domain_confidence_statistics = null;
    }

    public boolean isCombinable( final DomainId protein_domain ) {
        return getCombiningDomains().containsKey( protein_domain );
    }

    public void setKeyDomainConfidenceDescriptiveStatistics( final DescriptiveStatistics key_domain_confidence_statistics ) {
        _key_domain_confidence_statistics = key_domain_confidence_statistics;
    }

    public void setKeyDomainCount( final int key_domain_count ) {
        _key_domain_count = key_domain_count;
    }

    public void setKeyDomainProteinsCount( final int key_domain_proteins_count ) {
        _key_domain_proteins_count = key_domain_proteins_count;
    }

    @Override
    public List<BinaryDomainCombination> toBinaryDomainCombinations() {
        final List<BinaryDomainCombination> binary_combinations = new ArrayList<BinaryDomainCombination>( getNumberOfCombinableDomains() );
        for( final DomainId domain : getCombiningDomains().keySet() ) {
            binary_combinations.add( new BasicBinaryDomainCombination( getKeyDomain(), domain ) );
        }
        return binary_combinations;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append( getKeyDomain() );
        sb.append( " [" );
        sb.append( getKeyDomainCount() );
        sb.append( ", " );
        sb.append( getKeyDomainProteinsCount() );
        sb.append( ", " );
        sb.append( getNumberOfCombinableDomains() );
        sb.append( "]: " );
        sb.append( getCombiningDomainIdsAsStringBuilder() );
        return sb.toString();
    }
}
