// $Id: CombinableDomains.java,v 1.7 2009/10/26 23:29:40 cmzmasek Exp $
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

import java.util.List;
import java.util.SortedMap;

import org.forester.util.DescriptiveStatistics;

public interface CombinableDomains {

    /**
     * To add a new combinable domain.
     * 
     * @param protein_domain
     */
    public void addCombinableDomain( final DomainId protein_domain );

    /**
     * 
     * This must return all domains in this set of combinable domains (i.e.
     * the key domain and all domains which can combine with the key domain).
     * 
     *  @return all domains
     */
    List<DomainId> getAllDomains();

    List<DomainId> getCombinableDomains();

    /**
     * Returns the combinable domain identifiers sorted in alphabetical manner: -
     * keys are the combinable domain identifiers - values are the counts of
     * proteins exhibiting a particular combination
     * 
     * @return combining domain identifiers sorted in alphabetical manner
     */
    public SortedMap<DomainId, Integer> getCombinableDomainsIds();

    public StringBuilder getCombiningDomainIdsAsStringBuilder();

    /**
     * Returns the domain whose combinable domains are in stored in this
     * combinable domains.
     * 
     * @return the domain identifier
     */
    public DomainId getKeyDomain();

    /**
     * Gets descriptive statistics for the confidence (i.e. E-values) of the key
     * domain.
     * 
     * 
     * @return descriptive statistics for the confidence of the key domain
     */
    public DescriptiveStatistics getKeyDomainConfidenceDescriptiveStatistics();

    /**
     * Returns how many times the key domain is present in a given species
     * genome.
     * 
     * @return key domain count in species
     */
    public int getKeyDomainCount();

    /**
     * Returns how many proteins with the key domain are present in a given
     * species genome.
     * 
     * @return key domain proteins count in species
     */
    public int getKeyDomainProteinsCount();

    public int getNumberOfCombinableDomains();

    public int getNumberOfProteinsExhibitingCombination( final DomainId protein_domain );

    /**
     * Returns the species of this combinable domains.
     * 
     * @return the species
     */
    public Species getSpecies();

    public boolean isCombinable( final DomainId protein_domain );

    /**
     * This is to set descriptive statistics for the confidence (i.e. E-values)
     * of the key domain.
     * 
     * 
     * @param statistics
     */
    void setKeyDomainConfidenceDescriptiveStatistics( final DescriptiveStatistics statistics );

    /**
     * Sets how many times the key domain is present in a given species genome.
     * 
     * @param key_domain_count
     *            key domain count in species
     */
    void setKeyDomainCount( final int key_domain_count );

    /**
     * Sets how many proteins with the key domain are present in a given species
     * genome.
     * 
     * @param key_domain_proteins_count
     *            key domain protein count in species
     */
    void setKeyDomainProteinsCount( final int key_domain_proteins_count );

    public List<BinaryDomainCombination> toBinaryDomainCombinations();
}