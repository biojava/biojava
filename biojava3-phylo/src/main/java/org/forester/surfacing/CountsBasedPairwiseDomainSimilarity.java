// $Id: CopyNumberBasedPairwiseSimilarity.java,v 1.2 2007/10/01 23:57:45
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

public class CountsBasedPairwiseDomainSimilarity implements PairwiseDomainSimilarity {

    private final double _score;
    private final int    _copy_number_difference;

    /**
     * counts_difference: (counts for domain 1) minus (counts for domain 2).
     * 
     * 
     * @param counts_difference value of domain_1 minus value of domain_2
     * @param counts_sum
     */
    public CountsBasedPairwiseDomainSimilarity( final int counts_difference, final int counts_sum ) {
        if ( counts_sum <= 0 ) {
            throw new IllegalArgumentException( "attempt to use copy sum of less than or equal to 0: " + counts_sum );
        }
        _copy_number_difference = counts_difference;
        final int abs_copy_number_difference = Math.abs( counts_difference );
        if ( abs_copy_number_difference > counts_sum ) {
            throw new IllegalArgumentException( "attempt to use absolute copy number difference larger than copy number sum" );
        }
        _score = 1.0 - ( double ) abs_copy_number_difference / counts_sum;
    }

    /**
     * Returns (counts for domain 1) minus (counts for domain 2).
     * 
     */
    public int getDifferenceInCounts() {
        return _copy_number_difference;
    }

    public double getSimilarityScore() {
        return _score;
    }
}
