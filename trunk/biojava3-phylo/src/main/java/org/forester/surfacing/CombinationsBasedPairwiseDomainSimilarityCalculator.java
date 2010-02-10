// $Id: CombinationsBasedPairwiseSimilarityCalculator.java,v 1.4 2007/10/10
// 22:43:35 cmzmasek Exp $
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

public class CombinationsBasedPairwiseDomainSimilarityCalculator implements PairwiseDomainSimilarityCalculator {

    public PairwiseDomainSimilarity calculateSimilarity( final CombinableDomains domains_1,
                                                         final CombinableDomains domains_2 ) {
        if ( !domains_1.getKeyDomain().equals( domains_2.getKeyDomain() ) ) {
            throw new IllegalArgumentException( "attempt to calculate similarity between domain collection with different keys" );
        }
        final List<DomainId> d1 = domains_1.getCombinableDomains();
        final List<DomainId> d2 = domains_2.getCombinableDomains();
        int same = 0;
        int different = 0;
        for( final DomainId domain : d1 ) {
            if ( d2.contains( domain ) ) {
                same++;
            }
            else {
                different++;
            }
        }
        for( final DomainId domain : d2 ) {
            if ( !( d1.contains( domain ) ) ) {
                different++;
            }
        }
        final int difference = domains_1.getNumberOfCombinableDomains() - domains_2.getNumberOfCombinableDomains();
        return new CombinationsBasedPairwiseDomainSimilarity( same, different, difference );
    }
}
