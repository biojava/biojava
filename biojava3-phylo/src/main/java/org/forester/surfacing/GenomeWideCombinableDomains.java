// $Id: CombinableDomainsCollection.java,v 1.1 2007/10/10 22:43:36 cmzmasek Exp
// $
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

import java.util.SortedMap;
import java.util.SortedSet;

import org.forester.surfacing.BinaryDomainCombination.DomainCombinationType;
import org.forester.util.DescriptiveStatistics;

public interface GenomeWideCombinableDomains {

    public boolean contains( DomainId key_id );

    public CombinableDomains get( DomainId key_id );

    public SortedMap<DomainId, CombinableDomains> getAllCombinableDomainsIds();

    /**
     * This should return all domains ids present in the genome.
     * 
     * @return a sorted set of domains ids
     */
    public SortedSet<DomainId> getAllDomainIds();

    public DomainCombinationType getDomainCombinationType();

    SortedSet<DomainId> getMostPromiscuosDomain();

    /**
     * This should return a statistic for per domain 
     * promiscuity in a genome.
     * 
     * @return descriptive statistics for per domain promiscuity in a genome
     */
    public DescriptiveStatistics getPerGenomeDomainPromiscuityStatistics();

    public int getSize();

    public Species getSpecies();

    /**
     * This should return all binary domain combinations present in the genome.
     * 
     * @return a sorted set of binary domain combinations
     */
    public SortedSet<BinaryDomainCombination> toBinaryDomainCombinations();

    public StringBuilder toStringBuilder( GenomeWideCombinableDomainsSortOrder order );

    public static enum GenomeWideCombinableDomainsSortOrder {
        ALPHABETICAL_KEY_ID, KEY_DOMAIN_PROTEINS_COUNT, KEY_DOMAIN_COUNT, COMBINATIONS_COUNT
    }
}
