// $Id: AsciiHistogram.java,v 1.11 2009/10/26 23:29:40 cmzmasek Exp $
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

package org.forester.util;

public class AsciiHistogram {

    final private DescriptiveStatistics _stats;
    final private String                _title;

    public AsciiHistogram( final DescriptiveStatistics stats ) {
        _stats = stats;
        _title = "";
    }

    public AsciiHistogram( final DescriptiveStatistics stats, final String title ) {
        _stats = stats;
        _title = title;
    }

    private void drawToStringBuffer( final double min,
                                     final char symbol,
                                     final int size,
                                     final int digits,
                                     final StringBuffer sb,
                                     final int[] bins,
                                     final int max_count,
                                     final int under,
                                     final int over,
                                     final double binning_factor ) {
        final double draw_factor = ( double ) max_count / size;
        final int counts_size = ForesterUtil.roundToInt( Math.log10( max_count ) ) + 1;
        if ( !ForesterUtil.isEmpty( getTitle() ) ) {
            sb.append( getTitle() );
            sb.append( ForesterUtil.LINE_SEPARATOR );
            sb.append( ForesterUtil.LINE_SEPARATOR );
        }
        if ( under > 0 ) {
            sb.append( "[" + under + "] " );
            sb.append( ForesterUtil.LINE_SEPARATOR );
        }
        for( int i = 0; i < bins.length; ++i ) {
            final int count = bins[ i ];
            final double label = ForesterUtil.round( ( min + i * ( 1.0 / binning_factor ) ), digits );
            sb.append( ForesterUtil.pad( label + "", digits, '0', false ) );
            sb.append( " [" + ForesterUtil.pad( count + "", counts_size, ' ', true ) + "] " );
            final int s = ForesterUtil.roundToInt( count / draw_factor );
            for( int j = 0; j < s; ++j ) {
                sb.append( symbol );
            }
            sb.append( ForesterUtil.LINE_SEPARATOR );
        }
        if ( over > 0 ) {
            sb.append( "[" + over + "] " );
            sb.append( ForesterUtil.LINE_SEPARATOR );
        }
    }

    private DescriptiveStatistics getDescriptiveStatistics() {
        return _stats;
    }

    private String getTitle() {
        return _title;
    }

    public StringBuffer toStringBuffer( final double min,
                                        final double max,
                                        final int number_of_bins,
                                        final char symbol,
                                        final int size,
                                        final int digits ) {
        if ( min >= max ) {
            throw new IllegalArgumentException( "min [" + min + "] is larger than or equal to max [" + max + "]" );
        }
        if ( number_of_bins < 3 ) {
            throw new IllegalArgumentException( "number of bins is smaller than 3" );
        }
        if ( size < 2 ) {
            throw new IllegalArgumentException( "size is smaller than 2" );
        }
        final StringBuffer sb = new StringBuffer();
        int max_count = 0;
        final double binning_factor = number_of_bins / ( max - min );
        final int[] bins = BasicDescriptiveStatistics
                .performBinning( getDescriptiveStatistics().getDataAsDoubleArray(), min, max, number_of_bins );
        for( final int bin : bins ) {
            if ( bin > max_count ) {
                max_count = bin;
            }
        }
        drawToStringBuffer( min, symbol, size, digits, sb, bins, max_count, 0, 0, binning_factor );
        return sb;
    }

    public StringBuffer toStringBuffer( final int bins, final char symbol, final int size, final int digits ) {
        return toStringBuffer( getDescriptiveStatistics().getMin(),
                               getDescriptiveStatistics().getMax(),
                               bins,
                               symbol,
                               size,
                               digits );
    }
}
