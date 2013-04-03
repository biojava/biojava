// $Id: DescriptiveStatistics.java,v 1.12 2009/01/13 19:49:31 cmzmasek Exp $
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

public interface DescriptiveStatistics {

    public final static String PLUS_MINUS = "" + ( char ) 177;

    public abstract void addValue( final double d );

    public abstract double arithmeticMean();

    public abstract String asSummary();

    /**
     * Computes the coefficient of variation. Used to express standard deviation
     * independent of units of measure.
     * 
     * @return
     */
    public abstract double coefficientOfVariation();

    public abstract double[] getDataAsDoubleArray();

    public abstract double getMax();

    public abstract double getMin();

    public abstract int getN();

    public abstract double getSum();

    public abstract String getSummaryAsString();

    public abstract double getValue( final int index );

    public abstract double median();

    public abstract double midrange();

    /**
     * Determines relationship between the mean and the median. This reflects
     * how the data differs from the normal bell shaped distribution.
     * 
     * @return
     */
    public abstract double pearsonianSkewness();

    public abstract double sampleStandardDeviation();

    public abstract double sampleStandardUnit( final double value );

    public abstract double sampleVariance();

    public abstract double standardErrorOfMean();

    public abstract double sumDeviations();

    public abstract String toString();
}