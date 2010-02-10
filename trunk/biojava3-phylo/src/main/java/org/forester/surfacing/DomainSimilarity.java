// $Id: DomainSimilarity.java,v 1.8 2009/01/13 19:49:32 cmzmasek Exp $
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

/*
 * This is to represent a measure of similarity between two or more domains from
 * different genomes.
 */
public interface DomainSimilarity extends Comparable<DomainSimilarity> {

    public SortedSet<DomainId> getCombinableDomainIds( final Species species_of_combinable_domain );;

    public DomainId getDomainId();

    /**
     * For pairwise similarities, this should return the "difference"; for example the difference in counts
     * for copy number based features (the same as getMaximalDifferenceInCounts(), or the number
     * of actually different domain combinations. 
     * For pairwise similarities, this should return the difference,
     * while for comparisons of more than two domains, this should return the maximal difference
     * 
     * 
     * 
     * @return
     */
    public int getMaximalDifference();

    /**
     * For pairwise similarities, this should return the difference in counts,
     * while for comparisons of more than two domains, this should return the maximal difference
     * in counts
     * 
     * 
     * @return the (maximal) difference in counts
     */
    public int getMaximalDifferenceInCounts();

    public double getMaximalSimilarityScore();

    public double getMeanSimilarityScore();

    public double getMinimalSimilarityScore();

    /**
     * This should return the number of pairwise distances used to calculate
     * this similarity score
     * 
     * @return the number of pairwise distances
     */
    public int getN();

    public SortedSet<Species> getSpecies();

    /**
     * This should return a map, which maps species names to
     * SpeciesSpecificDomainSimilariyData
     * 
     * 
     * @return SortedMap<String, SpeciesSpecificDomainSimilariyData>
     */
    public SortedMap<Species, SpeciesSpecificDomainSimilariyData> getSpeciesData();

    public double getStandardDeviationOfSimilarityScore();

    public StringBuffer toStringBuffer( final PrintableDomainSimilarity.PRINT_OPTION print_option );

    static public enum DomainSimilarityScoring {
        DOMAINS, PROTEINS, COMBINATIONS;
    }

    public static enum DomainSimilaritySortField {
        MIN, MAX, SD, MEAN, ABS_MAX_COUNTS_DIFFERENCE, MAX_COUNTS_DIFFERENCE, MAX_DIFFERENCE, SPECIES_COUNT, DOMAIN_ID,
    }
}
