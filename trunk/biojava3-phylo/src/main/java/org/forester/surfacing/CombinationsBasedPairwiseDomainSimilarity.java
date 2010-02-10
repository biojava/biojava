// $Id: CombinationsBasedPairwiseSimilarity.java,v 1.1 2007/10/01 23:02:37
// cmzmasek Exp $
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

public class CombinationsBasedPairwiseDomainSimilarity implements PairwiseDomainSimilarity {

    private final int    _same_domains;
    private final int    _different_domains;
    private final int    _difference_in_counts;
    private final double _score;

    public CombinationsBasedPairwiseDomainSimilarity( final int same_domains,
                                                      final int different_domains,
                                                      final int difference_in_counts ) {
        if ( ( same_domains < 0 ) || ( different_domains < 0 ) ) {
            throw new IllegalArgumentException( "attempt to use domain counts less than 0" );
        }
        _difference_in_counts = difference_in_counts;
        _same_domains = same_domains;
        _different_domains = different_domains;
        if ( _different_domains == 0 ) {
            _score = 1.0;
        }
        else {
            _score = ( double ) _same_domains / ( _different_domains + _same_domains );
        }
    }

    @Override
    public int getDifferenceInCounts() {
        return _difference_in_counts;
    }

    public int getNumberOfDifferentDomains() {
        return _different_domains;
    }

    public int getNumberOfSameDomains() {
        return _same_domains;
    }

    public double getSimilarityScore() {
        return _score;
    }
}
