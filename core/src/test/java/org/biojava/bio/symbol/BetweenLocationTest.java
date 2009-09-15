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

import junit.framework.TestCase;

/**
 * Automated tests for BetweenLocation operations.
 *
 * @author Greg Cox
 */
public class BetweenLocationTest extends TestCase
{
// Static variables

// Member variables
	private Location betweenFiveSix;
	private Location betweenFiveSeven;
	private Location betweenFiveTen;
	private Location betweenSevenTwelve;
	private Location betweenNineEleven;
	private Location betweenTenEleven;
	private Location rangeOneTen;
	private Location rangeFiveSix;

// Constructors and initialization
	public BetweenLocationTest(String theName)
	{
		super(theName);
	}

// Interface implementations

// Public methods
	/**
	 * Runs the unit tests defined here.
	 */
	public static void main(String args[])
	{
		junit.textui.TestRunner.run(BetweenLocationTest.class);
	}

	// Test intersection operator
	public void testIntersectionOneTenBetweenFiveSix()
	{
		Location expected = new BetweenLocation(new RangeLocation(5, 6));

		Location result = LocationTools.intersection(rangeOneTen, betweenFiveSix);

		this.assertEqual(expected, result);
	}

	public void testIntersectionOneTenBetweenFiveSeven()
	{
		Location expected = new BetweenLocation(new RangeLocation(5, 7));

		Location result = LocationTools.intersection(rangeOneTen, betweenFiveSeven);

		this.assertEqual(expected, result);
	}

	public void testIntersectionOneTenBetweenTenEleven()
	{
		Location expected = Location.empty;

		Location result = LocationTools.intersection(rangeOneTen, betweenTenEleven);

		this.assertEqual(expected, result);
	}

	public void testIntersectionOneTenBetweenNineEleven()
	{
		Location expected = new BetweenLocation(new RangeLocation(9, 10));

		Location result = LocationTools.intersection(rangeOneTen, betweenNineEleven);

		this.assertEqual(expected, result);
	}

	public void testIntersectionBetweenFiveSixBetweenFiveSix()
	{
		Location expected = new BetweenLocation(new RangeLocation(5, 6));

		Location result = LocationTools.intersection(betweenFiveSix, betweenFiveSix);

		this.assertEqual(expected, result);
	}

	public void testIntersectionBetweenFiveTenBetweenFiveSix()
	{
		Location expected = new BetweenLocation(new RangeLocation(5, 6));

		Location result = LocationTools.intersection(betweenFiveTen, betweenFiveSix);

		this.assertEqual(expected, result);
	}

	public void testIntersectionBetweenFiveTenBetweenSevenTwelve()
	{
		Location expected = new BetweenLocation(new RangeLocation(7, 10));

		Location result = LocationTools.intersection(betweenFiveTen, betweenSevenTwelve);

		this.assertEqual(expected, result);
	}

	// Test overlaps operator
	public void testOverlapsOneTenBetweenFiveSix()
	{
		assertTrue(LocationTools.overlaps(rangeOneTen, betweenFiveSix));
	}

	public void testOverlapsOneTenBetweenFiveSeven()
	{
		assertTrue(LocationTools.overlaps(rangeOneTen, betweenFiveSeven));
	}

	public void testOverlapsOneTenBetweenTenEleven()
	{
		assertTrue(!LocationTools.overlaps(rangeOneTen, betweenTenEleven));
	}

	public void testOverlapsOneTenBetweenNineEleven()
	{
		assertTrue(LocationTools.overlaps(rangeOneTen, betweenNineEleven));
	}

	public void testOverlapsBetweenFiveSixBetweenFiveSix()
	{
		assertTrue(LocationTools.overlaps(betweenFiveSix, betweenFiveSix));
	}

	public void testOverlapsBetweenFiveTenBetweenFiveSix()
	{
		assertTrue(LocationTools.overlaps(betweenFiveTen, betweenFiveSix));
	}

	public void testOverlapsBetweenFiveTenBetweenSevenTwelve()
	{
		assertTrue(LocationTools.overlaps(betweenFiveTen, betweenSevenTwelve));
	}

	// Test equals operator
	public void testAreEqualRangeFiveSixBetweenFiveSix()
	{
		assertTrue(!LocationTools.areEqual(rangeFiveSix, betweenFiveSix));
	}

	public void testAreEqualBetweenFiveSixBetweenFiveSix()
	{
		assertTrue(LocationTools.areEqual(betweenFiveSix, betweenFiveSix));
	}

	public void testAreEqualBetweenFiveSixBetweenFiveTen()
	{
		assertTrue(!LocationTools.areEqual(betweenFiveSix, betweenFiveTen));
	}

	// Test contains operator
	public void testContainsOneTenBetweenFiveSix()
	{
		assertTrue(LocationTools.contains(rangeOneTen, betweenFiveSix));
	}

	public void testContainsBetweenFiveSixOneTen()
	{
		assertTrue(!LocationTools.contains(betweenFiveSix, rangeOneTen));
	}

	public void testContainsOneTenBetweenFiveSeven()
	{
		assertTrue(LocationTools.contains(rangeOneTen, betweenFiveSeven));
	}

	public void testContainsOneTenBetweenTenEleven()
	{
		assertTrue(!LocationTools.contains(rangeOneTen, betweenTenEleven));
	}

	public void testContainsOneTenBetweenNineEleven()
	{
		assertTrue(!LocationTools.contains(rangeOneTen, betweenNineEleven));
	}

	public void testContainsBetweenFiveSixBetweenFiveSix()
	{
		assertTrue(LocationTools.contains(betweenFiveSix, betweenFiveSix));
	}

	public void testContainsBetweenFiveTenBetweenFiveSix()
	{
		assertTrue(LocationTools.contains(betweenFiveTen, betweenFiveSix));
	}

	public void testContainsBetweenFiveSixBetweenFiveTen()
	{
		assertTrue(!LocationTools.contains(betweenFiveSix, betweenFiveTen));
	}

	public void testContainsBetweenFiveTenBetweenSevenTwelve()
	{
		assertTrue(!LocationTools.contains(betweenFiveTen, betweenSevenTwelve));
	}

	public void testContainsBetweenSevenTwelveBetweenFiveTen()
	{
		assertTrue(!LocationTools.contains(betweenSevenTwelve, betweenFiveTen));
	}

	// Test union operator
	public void testUnionRangeOneTenBetweenFiveSix()
	{
	}

	public void testUnionRangeOneTenBetweenFiveSeven()
	{
	}

	public void testUnionRangeOneTenBetweenTenEleven()
	{
	}

	public void testUnionRangeOneTenBetweenNineEleven()
	{
	}

	public void testUnionBetweenFiveSixBetweenFiveSix()
	{
	}

	public void testUnionBetweenFiveTenBetweenFiveSix()
	{
	}

	public void testUnionBetweenFiveTenBetweenSevenTwelve()
	{
	}

// Protected methods
	protected void setUp()
	{
		betweenFiveSix = new BetweenLocation(new RangeLocation(5, 6));
		betweenFiveSeven = new BetweenLocation(new RangeLocation(5, 7));
		betweenFiveTen = new BetweenLocation(new RangeLocation(5, 10));
		betweenSevenTwelve = new BetweenLocation(new RangeLocation(7, 12));
		betweenNineEleven = new BetweenLocation(new RangeLocation(9, 11));
		betweenTenEleven = new BetweenLocation(new RangeLocation(10, 11));
		rangeOneTen = new RangeLocation(1, 10);
		rangeFiveSix = new RangeLocation(5, 6);
	}

// Private methods
	/**
	 * This method assumes that the only decorator to deal with is between
	 * location.  Because of sheer lazyness, Between(5,10) and Circular(5,10)
	 * will return as equal from this method.
	 */
	private void assertEqual(Location theLocationA, Location theLocationB)
	{
		boolean locationAWrapped = false;
		boolean locationBWrapped = false;
		try
		{
			// Unwrap decorations
			if(LocationTools.isDecorated(theLocationA))
			{
				locationAWrapped = true;
				theLocationA = ((AbstractLocationDecorator)theLocationA).getWrapped();
			}

			if(LocationTools.isDecorated(theLocationB))
			{
				locationBWrapped = true;
				theLocationB = ((AbstractLocationDecorator)theLocationB).getWrapped();
			}
		}
		catch(ClassCastException cce)
		{
			throw new org.biojava.bio.BioError("Not possible: decorated by other than abstract location decorator");
		}
		assertTrue("The locations do not have the same decorators", locationAWrapped == locationBWrapped);

		// Delegate to LocationTools
		boolean result = LocationTools.areEqual(theLocationA, theLocationB);
		if(result == false)
		{
			fail("Locations not equal: " + theLocationA + " and " + theLocationB);
		}
	}
}
