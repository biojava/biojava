// $Id: ExternalNodeBasedCoverage.java,v 1.5 2009/01/13 19:49:31 cmzmasek Exp $
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

package org.forester.pccx;

import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

/*
 * @author Christian M. Zmasek
 */
public class ExternalNodeBasedCoverage implements Coverage {

    private final double _av_normalized_score;
    private final double _av_raw_score;
    private final int    _n;
    private final double _sd;
    private final double _max;
    private final double _min;

    public ExternalNodeBasedCoverage( final DescriptiveStatistics stats,
                                      final double average_raw_score,
                                      final CoverageCalculationOptions options ) {
        _av_normalized_score = stats.arithmeticMean();
        _av_raw_score = average_raw_score;
        _n = stats.getN();
        if ( _n > 1 ) {
            _sd = stats.sampleStandardDeviation();
        }
        else {
            _sd = 0.0;
        }
        _max = stats.getMax();
        _min = stats.getMin();
    }

    public String asString() {
        final StringBuffer sb = new StringBuffer();
        if ( getN() == 1 ) {
            sb.append( "Normalized score: " + getScore() + ForesterUtil.getLineSeparator() );
            sb.append( "Raw score       : " + getAvarageRawScore() );
        }
        else {
            sb.append( "Avarage normalized score: " + getScore() + " [sd=" + getSD() + " min=" + getMin() + " max="
                    + getMax() + " n=" + getN() + "]" + ForesterUtil.getLineSeparator() );
            sb.append( "Avarage raw score       : " + getAvarageRawScore() );
        }
        return sb.toString();
    }

    public double getAvarageNormalizedScore() {
        return _av_normalized_score;
    }

    public double getAvarageRawScore() {
        return _av_raw_score;
    }

    public double getMax() {
        return _max;
    }

    public double getMin() {
        return _min;
    }

    public int getN() {
        return _n;
    }

    public double getScore() {
        return getAvarageNormalizedScore();
    }

    public double getSD() {
        return _sd;
    }
}
