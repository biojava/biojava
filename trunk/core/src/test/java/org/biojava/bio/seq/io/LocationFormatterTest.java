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

package org.biojava.bio.seq.io;

import junit.framework.TestCase;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.BetweenLocation;
import org.biojava.bio.symbol.CircularLocation;
import org.biojava.bio.symbol.FuzzyLocation;
import org.biojava.bio.symbol.FuzzyPointLocation;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.RangeLocation;

/**
 * This test suite is built on JUnit (www.junit.org).  It tests the
 * Genbank (and EMBL, as they use the same code) Location formatter
 * against various locations. The formatting code being tested is
 * located in AbstractGenEmblFileFormer.
 *
 * @author Greg Cox
 * @author Keith James
 */
public class LocationFormatterTest extends TestCase
{
    // Static variables
    protected static final StrandedFeature.Strand  PLUS = StrandedFeature.POSITIVE;
    protected static final StrandedFeature.Strand MINUS = StrandedFeature.NEGATIVE;

    // Member variables
    protected SeqFileFormer fileFormer;

    protected StringBuffer sb;

    protected Location mSimplePointLocation;
    // 467
    protected Location mRangedPointLocation;
    // (102.110)
    protected Location mBoundedPointLocation;
    // >12

    protected Location mSimpleRangedLocation;
    // 340..565
    protected Location mBoundedRangedLocation;
    // <345..500
    protected Location mBoundedPointRangedPointRangedLocation;
    // <94..(200.250)
    protected Location mRangedPointSimplePointRangedLocation;
    // (23.45)..600
    protected Location mRangedPointRangedPointRangedLocation;
    // (122.133)..(204.221)

    protected Location mCircularSimpleRangedLocation;
    // 340..565
    protected Location mCircularBoundedRangedLocation;
    // <345..500
    protected Location mCircularBoundedPointRangedPointRangedLocation;
    // <94..(200.250)
    protected Location mCircularRangedPointSimplePointRangedLocation;
    // (23.45)..600
    protected Location mCircularRangedPointRangedPointRangedLocation;
    // (122.133)..(204.221)

    //	protected Location mJoinCompoundLocation;
    // join(>12,(122.133)..(204.221),340..565,<345..500)
    //	protected Location mOrderCompoundLocation;
    // order(>12,(122.133)..(204.221),340..565,<345..500)

    protected Location mAdjacentBetweenLocation;
    // 20^21
    protected Location mFuzzyBetweenLocation;
    // 20^30

    /**
     * Creates a LocationFormatterTest object.
     * @param name The name of the test case
     * @see TestCase
     */
    public LocationFormatterTest(String theString)
    {
        super(theString);
        sb = new StringBuffer(512);
    }

    // Interface methods

    // Public methods

    /**
     * Runs the unit tests defined here.
     */
    public static void main(String args[])
    {
        junit.textui.TestRunner.run(LocationFormatterTest.class);
    }

    /**
     * Checks that mSimplePointLocation is formatted correctly under genbank.
     */
    public void testSimplePointLocationGenbank()
    {
        String expected = "467";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mSimplePointLocation, PLUS).toString();
        assertEquals(expected, actual);
    }

    /**
     * Checks that mRangedPointLocation is formatted correctly under
     * genbank.
     */
    public void testRangedPointLocationGenbank()
    {
        String expected = "(102.110)";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mRangedPointLocation, PLUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that mBoundedPointLocation is formatted correctly under genbank.
     */
    public void testBoundedPointLocationGenbank()
    {
        String expected = ">12";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mBoundedPointLocation, PLUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that mSimpleRangedLocation is formatted correctly under
     * genbank.
     */
    public void testSimpleRangedLocationGenbank()
    {
        String expected = "340..565";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mSimpleRangedLocation, PLUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that mBoundedRangedLocation is formatted correctly under genbank.
     */
    public void testBoundedRangedLocationGenbank()
    {
        String expected = "<345..500";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mBoundedRangedLocation, PLUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that mBoundedPointRangedPointRangedLocation is formatted
     * correctly under genbank.
     */
    public void testBoundedPointRangedPointRangedLocationGenbank()
    {
        String expected = "<94..(200.250)";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mBoundedPointRangedPointRangedLocation, PLUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that mRangedPointSimplePointRangedLocation is formatted
     * correctly under genbank.
     */
    public void testRangedPointSimplePointRangedLocationGenbank()
    {
        String expected = "(23.45)..600";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mRangedPointSimplePointRangedLocation, PLUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that mRangedPointRangedPointRangedLocation is formatted
     * correctly under genbank.
     */
    public void testRangedPointRangedPointRangedLocationGenbank()
    {
        String expected = "(122.133)..(204.221)";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mRangedPointRangedPointRangedLocation, PLUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that mCircularSimpleRangedLocation is formatted
     * correctly under genbank.
     */
    public void testCircularSimpleRangedLocationGenbank()
    {
        String expected = "340..565";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mCircularSimpleRangedLocation, PLUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that mCircularBoundedRangedLocation is formatted
     * correctly under genbank.
     */
    public void testCircularBoundedRangedLocationGenbank()
    {
        String expected = "<345..500";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mCircularBoundedRangedLocation, PLUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that mCircularBoundedPointRangedPointRangedLocation is
     * formatted correctly under genbank.
     */
    public void testCircularBoundedPointRangedPointRangedLocationGenbank()
    {
        String expected = "<94..(200.250)";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mCircularBoundedPointRangedPointRangedLocation, PLUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that mCircularRangedPointSimplePointRangedLocation is
     * formatted correctly under genbank.
     */
    public void testCircularRangedPointSimplePointRangedLocationGenbank()
    {
        String expected = "(23.45)..600";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mCircularRangedPointSimplePointRangedLocation, PLUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that mCircularRangedPointRangedPointRangedLocation is
     * formatted correctly under genbank.
     */
    public void testCircularRangedPointRangedPointRangedLocationGenbank()
    {
        String expected = "(122.133)..(204.221)";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mCircularRangedPointRangedPointRangedLocation, PLUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that mJoinCompoundLocation is formatted correctly under genbank.
     */
    //	public void testJoinCompoundLocationGenbank()
    //	{
    //		String expected = "join(>12,(122.133)..(204.221),340..565,<345..500)";
    //
    //		String actual = fileFormer.formatLocation(mJoinCompoundLocation, PLUS);
    //
    //		this.assertEquals(expected, actual);
    //	}
    //
    //	/**
    //	 * Checks that mOrderCompoundLocation is formatted correctly under genbank.
    //	 */
    //	public void testOrderCompoundLocationGenbank()
    //	{
    //		String expected = "order(>12,(122.133)..(204.221),340..565,<345..500)";
    //
    //		String actual = fileFormer.formatLocation(mOrderCompoundLocation, PLUS);
    //
    //		this.assertEquals(expected, actual);
    //	}
    //
    //	/**
    //	 * Checks that a complement mJoinCompoundLocation is formatted correctly under genbank.
    //	 */
    //	public void testComplementJoinCompoundLocationGenbank()
    //	{
    //		String expected = "complement(join(>12,(122.133)..(204.221),340..565,<345..500))";
    //
    //		String actual = fileFormer.formatLocation(mJoinCompoundLocation, MINUS);
    //
    //		this.assertEquals(expected, actual);
    //	}
    //
    //	/**
    //	 * Checks that a complement mOrderCompoundLocation is formatted correctly under genbank.
    //	 */
    //	public void testComplementOrderCompoundLocationGenbank()
    //	{
    //		String expected = "complement(order(>12,(122.133)..(204.221),340..565,<345..500))";
    //
    //		String actual = fileFormer.formatLocation(mOrderCompoundLocation, MINUS);
    //
    //		this.assertEquals(expected, actual);
    //	}

    /**
     * Checks that a complement mSimplePointLocation is formatted
     * correctly under genbank.
     */
    public void testComplementSimplePointLocationGenbank()
    {
        String expected = "complement(467)";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mSimplePointLocation, MINUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that a complement mRangedPointLocation is formatted
     * correctly under genbank.
     */
    public void testComplementRangedPointLocationGenbank()
    {
        String expected = "complement((102.110))";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mRangedPointLocation, MINUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that a complement mBoundedPointLocation is formatted
     * correctly under genbank.
     */
    public void testComplementBoundedPointLocationGenbank()
    {
        String expected = "complement(>12)";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mBoundedPointLocation, MINUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that a complement mSimpleRangedLocation is formatted
     * correctly under genbank.
     */
    public void testComplementSimpleRangedLocationGenbank()
    {
        String expected = "complement(340..565)";

        sb.delete(0, sb.length() );
        String actual = fileFormer.formatLocation(sb, mSimpleRangedLocation, MINUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that a complement mBoundedRangedLocation is formatted
     * correctly under genbank.
     */
    public void testComplementBoundedRangedLocationGenbank()
    {
        String expected = "complement(<345..500)";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mBoundedRangedLocation, MINUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that a complement mBoundedPointRangedPointRangedLocation
     * is formatted correctly under genbank.
     */
    public void testComplementBoundedPointRangedPointRangedLocationGenbank()
    {
        String expected = "complement(<94..(200.250))";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mBoundedPointRangedPointRangedLocation, MINUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that a complement mRangedPointSimplePointRangedLocation
     * is formatted correctly under genbank.
     */
    public void testComplementRangedPointSimplePointRangedLocationGenbank()
    {
        String expected = "complement((23.45)..600)";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mRangedPointSimplePointRangedLocation, MINUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that a complement mRangedPointRangedPointRangedLocation
     * is formatted correctly under genbank.
     */
    public void testComplementRangedPointRangedPointRangedLocationGenbank()
    {
        String expected = "complement((122.133)..(204.221))";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mRangedPointRangedPointRangedLocation, MINUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that a complement mCircularSimpleRangedLocation is
     * formatted correctly under genbank.
     */
    public void testComplementCircularSimpleRangedLocationGenbank()
    {
        String expected = "complement(340..565)";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mCircularSimpleRangedLocation, MINUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that a complement mCircularBoundedRangedLocation is
     * formatted correctly under genbank.
     */
    public void testComplementCircularBoundedRangedLocationGenbank()
    {
        String expected = "complement(<345..500)";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mCircularBoundedRangedLocation, MINUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that a complement
     * mCircularBoundedPointRangedPointRangedLocation is formatted
     * correctly under genbank.
     */
    public void testComplementCircularBoundedPointRangedPointRangedLocationGenbank()
    {
        String expected = "complement(<94..(200.250))";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mCircularBoundedPointRangedPointRangedLocation, MINUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that a complement
     * mCircularRangedPointSimplePointRangedLocation is formatted
     * correctly under genbank.
     */
    public void testComplementCircularRangedPointSimplePointRangedLocationGenbank()
    {
        String expected = "complement((23.45)..600)";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mCircularRangedPointSimplePointRangedLocation, MINUS).toString();

        assertEquals(expected, actual);
    }

    /**
     * Checks that a complement
     * mCircularRangedPointRangedPointRangedLocation is formatted
     * correctly under genbank.
     */
    public void testComplementCircularRangedPointRangedPointRangedLocationGenbank()
    {
        String expected = "complement((122.133)..(204.221))";

        sb.delete(0, sb.length());
        String actual = fileFormer.formatLocation(sb, mCircularRangedPointRangedPointRangedLocation, MINUS).toString();

        assertEquals(expected, actual);
    }

    //	/**
    //	 * Checks that a complement mJoinCompoundLocation is formatted correctly under genbank.
    //	 */
    //	public void testJoinComplementCompoundLocationGenbank()
    //	{
    //		String expected = "join(complement(>12),complement((122.133)..(204.221)),complement(340..565),complement(<345..500))";
    //
    //		String actual = fileFormer.formatLocation(mJoinCompoundLocation, MINUS);
    //
    //		assertEquals(expected, actual);
    //	}
    //
    //	/**
    //	 * Checks that a complement mOrderCompoundLocation is formatted correctly under genbank.
    //	 */
    //	public void testOrderComplementCompoundLocationGenbank()
    //	{
    //		String expected = "order(complement(>12),complement((122.133)..(204.221)),complement(340..565),complement(<345..500))";
    //
    //		String actual = fileFormer.formatLocation(mOrderCompoundLocation, MINUS);
    //
    //		assertEquals(expected, actual);
    //	}

    /**
     * Declares member points and locations for use in testing the formatter.
     */
    protected void setUp()
    {
        fileFormer = new GenbankFileFormer();

        mSimplePointLocation = new PointLocation(467);
        mRangedPointLocation = new FuzzyPointLocation(102, 110, FuzzyPointLocation.RESOLVE_AVERAGE);
        mBoundedPointLocation = new FuzzyPointLocation(12, Integer.MAX_VALUE, FuzzyPointLocation.RESOLVE_AVERAGE);
        mSimpleRangedLocation = new RangeLocation(340, 565);
        mBoundedRangedLocation = new FuzzyLocation(Integer.MIN_VALUE, 500, 345, 500, FuzzyLocation.RESOLVE_INNER);
        mBoundedPointRangedPointRangedLocation = new FuzzyLocation(Integer.MIN_VALUE, 250, 94, 200, FuzzyLocation.RESOLVE_INNER);
        mRangedPointSimplePointRangedLocation = new FuzzyLocation(23, 600, 45, 600, FuzzyLocation.RESOLVE_AVERAGE);
        mRangedPointRangedPointRangedLocation = new FuzzyLocation(122, 221, 133, 204, FuzzyLocation.RESOLVE_AVERAGE);

        mCircularSimpleRangedLocation = new CircularLocation(mSimpleRangedLocation, 1000);
        mCircularBoundedRangedLocation = new CircularLocation(mBoundedRangedLocation, 1000);
        mCircularBoundedPointRangedPointRangedLocation = new CircularLocation(mBoundedPointRangedPointRangedLocation, 1000);
        mCircularRangedPointSimplePointRangedLocation = new CircularLocation(mRangedPointSimplePointRangedLocation, 1000);
        mCircularRangedPointRangedPointRangedLocation = new CircularLocation(mRangedPointRangedPointRangedLocation, 1000);
        mCircularRangedPointRangedPointRangedLocation = new CircularLocation(mRangedPointRangedPointRangedLocation, 1000);

        mAdjacentBetweenLocation = new BetweenLocation(new RangeLocation(20, 21));
        mFuzzyBetweenLocation = new BetweenLocation(new RangeLocation(20, 30));
    }
    
    // Private methods
    /**
     * Asserts that theString1 is equal to theString2.
     *
     * @param theString1 The first string to compare
     * @param theString2 The second string to compare
     */
    //	private void assertEquals(String theString1, String theString2)
    //	{
    //		if(theString1.equals(theString2) == false)
    //		{
    //			fail("\n\t" + theString1 + "\n\t\t != \n\t" + theString2);
    //		}
    //	}
}
