/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

package org.biojava.bio.symbol;

import java.util.AbstractList;
import java.util.Collections;
import java.util.Iterator;

import org.biojava.bio.BioError;

/**
 * <code>FuzzyPointLocation</code> represents two types of EMBL-style
 * partially-defined locations. These are the '(123.567)' type, which
 * represent a single residue somewhere between these coordinates and
 * the '<123' or '>123' type, which represent an unbounded location,
 * not including the residue at that coordinate.
 *
 * @author <a href="mailto:kdj@sanger.ac.uk">Keith James</a>
 * @author Greg Cox
 */
public class FuzzyPointLocation extends AbstractLocation
{
    // Use the minimum value
    public final static PointResolver RESOLVE_MIN;

    // Use the maximum value
    public final static PointResolver RESOLVE_MAX;

    // Use the arithmetic mean of the two values, unless they are
    // unbounded, in which case Integer.MIN_VALUE or Integer.MAX_VALUE
    // is returned
    public final static PointResolver RESOLVE_AVERAGE;

    static
    {
	RESOLVE_MIN     = new MinPointResolver();
	RESOLVE_MAX     = new MaxPointResolver();
	RESOLVE_AVERAGE = new AveragePointResolver();
    }

    private int           min;
    private int           max;
    private PointResolver resolver;

    /**
     * Creates a new <code>FuzzyPointLocation</code> object.  If the minimum
     * value is equal to the maximum value, it is interperted as a swissprot
     * point location such as "?24".  If the minimum is Integer.MIN_VALUE and
     * the maximum is Integer.MAX_VALUE, it is interperted as a swissprot
     * location like '?'.
     *
     * @param min an <code>int</code> value for the minimum boundary
     * of the location, Integer.MIN_VALUE if unbounded.
     * @param max an <code>int</code> value for the minimum boundary
     * of the location, Integer.MAX_VALUE if unbounded.
     * @param resolver a <code>PointResolver</code> which defines the
     * policy used to calculate the location's min and max
     * properties.
     */
    public FuzzyPointLocation(int min, int max, PointResolver resolver)
    {
	this.min      = min;
	this.max      = max;
	this.resolver = resolver;
    }

    public PointResolver getResolver()
    {
	return resolver;
    }

    public int getMin()
    {
	return min;
    }

    public int getMax()
    {
	return max;
    }

    public boolean hasBoundedMin()
    {
	return min != Integer.MIN_VALUE;
    }

    public boolean hasBoundedMax()
    {
	return max != Integer.MAX_VALUE;
    }

    public boolean overlaps(Location loc)
    {
	return loc.contains(this);
    }

    public boolean contains(Location loc)
    {
	// If the location is unbounded, it is not certain that it
	// contains any other specific location
	return (hasBoundedMin() && hasBoundedMax()) &&
	    (resolver.resolve(this) == loc.getMin()) &&
	    (resolver.resolve(this) == loc.getMax());
    }

    public boolean contains(int point)
    {
	// If the location is unbounded, it is not certain that it
	// contains any other specific coordinate
	return (hasBoundedMin() && hasBoundedMax()) &&
	    resolver.resolve(this) == point;
    }

    public boolean equals(Location loc)
    {
	return this.contains(loc) && loc.contains(this);
    }

    public int hashCode()
    {
	return getMin();
    }

    public Location intersection(Location loc)
    {
	return loc.contains(this)
	    ? this
	    : Location.empty;
    }


    public SymbolList symbols(SymbolList slist)
    {
	final Symbol sym = slist.symbolAt(resolver.resolve(this));
	try
	{
	    return new SimpleSymbolList(slist.getAlphabet(), new AbstractList()
		{
		    public Object get(int index)
			throws IndexOutOfBoundsException
		    {
			if (index == 0)
			{
			    return sym;
			}

			throw new IndexOutOfBoundsException("Index " + index + " greater than 0");
		    }

		    public int size()
		    {
			return 1;
		    }
		});
	}
	catch (IllegalSymbolException ise)
	{
	    throw new BioError(ise);
	}
    }

    public boolean isContiguous()
    {
	return true;
    }

    public Iterator blockIterator()
    {
	return Collections.singleton(this).iterator();
    }

    public Location translate(int dist)
    {
	if (dist == 0)
	    return this;

	return new FuzzyPointLocation(this.min + dist,
				      this.max + dist,
				      this.resolver);
    }

    public String toString()
    {
	if (hasBoundedMin() && hasBoundedMax())
	{
	    if (getMin() == getMax())
	    {
	    	return("?" + getMin());
	    }
	    else
	    {
		    return "["
			+ Integer.toString(getMin())
			+ "."
			+ Integer.toString(getMax());
		}
	}
	else if (hasBoundedMin())
	{
	    return "[>"
		+ Integer.toString(getMin())
		+ "]";
	}
	else if (hasBoundedMax())
	{
	    return "[<"
		+ Integer.toString(getMax())
		+ "]";
	}
	else
	{
            return "?";
	}
    }

     
    /**
     * Determines how a <code>FuzzyPointLocation</code> should be treated when used
     * as a normal <code>Location</code>.
     *
     * Use one of the implementations of this interface when creating a <code>FuzzyPointLocation</code>
     * to specify how the fuzzy (inner/outer) properties are translated into the standard
     * Location min and max properties.
     *
     * It is possible to write custom implementations of this to create <code>FuzzyLocations</code>
     * with exotic behaviour.
     */
    
    public static interface PointResolver
    {
        /**
         * Return the actual point that the specified location should claim to
         * occupy.
         */
        
        public int resolve(FuzzyPointLocation loc);
    }

    private static class MinPointResolver implements PointResolver
    {
	public int resolve(FuzzyPointLocation loc)
	{
	    if (loc.hasBoundedMin())
		return loc.getMin();
	    else
		return Integer.MIN_VALUE;
	}
    }

    private static class MaxPointResolver implements PointResolver
    {
	public int resolve(FuzzyPointLocation loc)
	{
	    if (loc.hasBoundedMax())
		return loc.getMax();
	    else
		return Integer.MAX_VALUE;
	}
    }

    private static class AveragePointResolver implements PointResolver
    {
	public int resolve(FuzzyPointLocation loc)
	{
	    // Range of form: (123.567)
	    if (loc.hasBoundedMin() && loc.hasBoundedMax())
	    {
		return (loc.getMin() + loc.getMax()) / 2;
	    }
		// Swissprot range ?
		else if((loc.hasBoundedMin() == false) &&
				(loc.hasBoundedMax() == false))
		{
			return 0;
		}
	    // Range of form: <123 or >123
	    else
	    {
		return loc.hasBoundedMin() ? Integer.MAX_VALUE : Integer.MIN_VALUE;
	    }
	}
    }
}
