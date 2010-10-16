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

//import [statements];

/**
 * Helper class for location operations involving between locations.
 *
 * @author Greg Cox
 */
final class BetweenLocationTools
{
// Static variables

// Member variables

// Constructors and initialization

// Interface implementations

// Public methods

// Protected methods
	protected static Location union(Location locationA, Location locationB)
	{
		if(BetweenLocationTools.isBetween(locationA) && (BetweenLocationTools.isBetween(locationB)))
		{
			Location tempLocationA = BetweenLocationTools.unwrap(locationA);
			Location tempLocationB = BetweenLocationTools.unwrap(locationB);
			Location mergeLocation = LocationTools.union(tempLocationA, tempLocationB);
			mergeLocation = new BetweenLocation(mergeLocation);
			return mergeLocation;
		}
		else
		{
			java.util.List locList = new java.util.ArrayList();
			locList.add(locationA);
			locList.add(locationB);
			return new CompoundLocation(locList);
		}
	}

	protected static Location intersection(Location locationA, Location locationB)
	{
		Location tempLocationA = BetweenLocationTools.unwrap(locationA);
		Location tempLocationB = BetweenLocationTools.unwrap(locationB);

		Location tempIntersection = LocationTools.intersection(tempLocationA, tempLocationB);
		// If the result is one nucleotide, there is a case like
		// 1..10 intersect 10^11.  This reduces to the empty location
		if(tempIntersection.getMin() == tempIntersection.getMax())
		{
			tempIntersection = Location.empty;
		}

		if(!LocationTools.areEqual(tempIntersection, Location.empty))
		{
			tempIntersection = new BetweenLocation(tempIntersection);
		}
		return tempIntersection;
	}

	protected static boolean overlaps(Location locationA, Location locationB)
	{
		boolean toReturn = true;

		Location tempLocation = BetweenLocationTools.intersection(locationA, locationB);
		if(LocationTools.areEqual(tempLocation, Location.empty))
		{
			toReturn = false;
		}
		return toReturn;
	}

	protected static boolean contains(Location locationA, Location locationB)
	{
		Location tempLocationA = BetweenLocationTools.unwrap(locationA);
		Location tempLocationB = BetweenLocationTools.unwrap(locationB);

		return LocationTools.contains(tempLocationA, tempLocationB);
	}

	protected static boolean areEqual(Location locationA, Location locationB)
	{
		boolean toReturn = false;
		// Check that both are between types
		if(BetweenLocationTools.isBetween(locationA) && BetweenLocationTools.isBetween(locationB))
		{
			Location tempLocationA = BetweenLocationTools.unwrap(locationA);
			Location tempLocationB = BetweenLocationTools.unwrap(locationB);
			toReturn = LocationTools.areEqual(tempLocationA, tempLocationB);
		}

		return toReturn;
	}

	protected static Location unwrap(Location theLoc)
	{
		Location toReturn = theLoc;
		try
		{
			if(BetweenLocationTools.isBetween(theLoc))
			{
				toReturn = ((BetweenLocation)theLoc).getWrapped();
			}
		}
		catch(ClassCastException cce)
		{
			throw new org.biojava.bio.BioError("Unwrapping not possible: decorated by other than between decorator");
		}

		return toReturn;
	}

	/**
	 * Tests if the location passed in is an instance of a between location
	 * or not.
	 *
	 * @param theLocation Location to test
	 * @return True if it is a between location, false otherwise
	 */
	protected static boolean isBetween(Location theLocation)
	{
		boolean toReturn = false;
	    try
	    {
		    if(theLocation.getDecorator(Class.forName("org.biojava.bio.symbol.BetweenLocation")) != null)
	    	{
	    		toReturn = true;
		    }
		}
		catch(Exception e)
		{
			throw new org.biojava.bio.BioError("class org.biojava.bio.symbol.BetweenLocation could not be loaded");
		}
    	return toReturn;
    }

// Private methods
}
