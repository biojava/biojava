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

import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;

/**
 * <code>CompoundLocationTest</code> tests the behaviour of
 * <code>CompoundLocation</code> by itself and combined with
 * <code>LocationTools</code>.
 *
 * @author <a href="mailto:kdj@sanger.ac.uk">Keith James</a>
 * @author Thomas Down
 * @author Francois Pepin
 */
public class CompoundLocationTest extends TestCase
{
    protected Location r1;
    protected Location r2;
    protected Location r3;
    protected Location r4;
    protected Location r5;
    protected Location r6;
    protected Location r7;
    protected Location r8;
    protected Location r9;


    private   List     locs1;
    private   List     locs2;
    private   List     locs3;
    private   List     locs4;
    private   List     locs5;
    private   List     locs6;
    private   List     locs7;
    private   List     locs8;
    private   List     locs9;



    public CompoundLocationTest(String name)
    {
	super(name);
    }

  /**
   * Runs the unit tests defined here.
   */
  public static void main(String args[])
  {
    junit.textui.TestRunner.run(CompoundLocationTest.class);
  }

    protected void setUp() throws Exception
    {
	locs1 = new ArrayList();
	locs1.add(new RangeLocation(1, 100));
	locs1.add(new RangeLocation(150, 200));

	locs2 = new ArrayList();
	locs2.add(new RangeLocation(150, 160));
	locs2.add(new RangeLocation(170, 190));

	locs3 = new ArrayList();
	locs3.add(new RangeLocation(250, 300));
	locs3.add(new RangeLocation(350, 400));

	locs4 = new ArrayList();
	locs5 = new ArrayList();
	locs6 = new ArrayList();
	locs7 = new ArrayList();
	for (int i = 0; i < 173; ++i) {
	    int j = i * 17;
	    locs4.add(new RangeLocation(j, j + 6));
	    locs5.add(new RangeLocation(j + 5, j + 8));
	    locs6.add(new RangeLocation(j + 5, j + 6));
	    locs7.add(new RangeLocation(j + 10, j + 11));
	}
        locs8 = new ArrayList();
        locs8.add(new RangeLocation(250, 255));
        locs8.add(new RangeLocation(350, 400));

        locs9 = new ArrayList();
        locs9.add(new RangeLocation(256, 300));
        locs9.add(new RangeLocation(350, 400));


	r1 = LocationTools.buildLoc(locs1);
	r2 = LocationTools.buildLoc(locs2);
	r3 = LocationTools.buildLoc(locs3);
	r4 = LocationTools.buildLoc(locs4);
	r5 = LocationTools.buildLoc(locs5);
	r6 = LocationTools.buildLoc(locs6);
	r7 = LocationTools.buildLoc(locs7);
        r8 = LocationTools.buildLoc(locs8);
        r9 = LocationTools.buildLoc(locs9);
    }

    /**
     * <code>testEquals</code> tests equality directly.
     *
     */
    public void testEquals()
    {
  	assertEquals(r1, r1);
  	assertEquals(r1, LocationTools.buildLoc(locs1));
    }

    /**
     * <code>testAreEqual</code> tests equality via
     * <code>LocationTools</code>.
     *
     */
    public void testAreEqual()
    {
	assertTrue(LocationTools.areEqual(r1, r1));
	assertEquals(r1, LocationTools.buildLoc(locs1));
    }

    /**
     * <code>testOverlaps</code> tests overlaps via
     * <code>LocationTools</code>.
     *
     */
    public void testOverlaps()
    {
	assertTrue(LocationTools.overlaps(r1, r2));
	assertTrue(! LocationTools.overlaps(r1, r3));
    }

    /**
     * <code>testContains</code> tests contains via
     * <code>LocationTools</code>.
     *
     */
    public void testContains()
    {
	assertTrue(LocationTools.contains(r1, r2));
	assertTrue(! LocationTools.contains(r1, r3));
	assertTrue(r1.contains(1));
	assertTrue(r1.contains(100));
	assertTrue(r1.contains(150));
	assertTrue(r1.contains(200));

	// Between contained ranges
	assertTrue(! r1.contains(101));
	// Outside contained ranges
	assertTrue(! r1.contains(201));
    }

    /**
     * <code>testIntersection</code> tests intersection via
     * <code>LocationTools</code>.
     *
     */
    public void testIntersection()
    {
	assertEquals(LocationTools.intersection(r1, r2),
		     LocationTools.intersection(r2, r1));
	assertEquals(LocationTools.intersection(r1, r3),
		     Location.empty);
	assertEquals(LocationTools.intersection(r3, r1),
		     Location.empty);
	assertEquals(LocationTools.intersection(r1, new RangeLocation(1, 1000)),
		     r1);
	assertEquals(LocationTools.intersection(new RangeLocation(1, 1000), r1),
		     r1);
    }

    /**
     * <code>testUnion</code> tests union via
     * <code>LocationTools</code>.
     *
     */
    public void testUnion()
    {
	assertEquals(r1, LocationTools.union(r1, r1));
	assertEquals(r2, LocationTools.union(r2, r2));
	assertEquals(LocationTools.union(r1, r2),
		     LocationTools.union(r2, r1));
	assertEquals(LocationTools.union(r1, r3),
                     LocationTools.union(r3, r1));
        assertEquals(r3, LocationTools.union(r8,r9));

    }

    /**
     * <code>testIsContiguous</code> tests contiguous.
     *
     */
    public void testIsContiguous()
    {
	assertTrue(! r1.isContiguous());

	List single = new ArrayList();
	single.add(new RangeLocation(1, 100));
	Location contig = LocationTools.buildLoc(single);

	assertTrue(contig.isContiguous());
    }

    public void testOverlapsBig()
    {
	assertTrue(LocationTools.overlaps(r4, r5));
	assertTrue(!LocationTools.overlaps(r4, r7));
    }

    public void testContainsBig()
    {
	assertTrue(LocationTools.contains(r4, r6));
	assertTrue(!LocationTools.contains(r4, r5));
	assertTrue(!LocationTools.contains(r4, r7));
    }

    public void testIntersectionBig() {
	assertEquals(LocationTools.intersection(r4, r5), r6);
	assertEquals(LocationTools.intersection(r4, r7), Location.empty);
    }

    public void testSubtract() {
        Location l1 = new RangeLocation(100, 200);
        Location l2 = LocationTools.union(
            new RangeLocation(100, 119),
            new RangeLocation(180, 200)
        );
        assertEquals(
            LocationTools.subtract(
                l1,
                new RangeLocation(120, 179)
            ),
            l2
        );
        assertEquals(
            LocationTools.subtract(
                l1,
                new RangeLocation(80, 119)
            ),
            new RangeLocation(120, 200)
        );
        assertEquals(
            LocationTools.subtract(
                l1,
                new RangeLocation(180, 1019)
            ),
            new RangeLocation(100, 179)
        );
        assertEquals(
            LocationTools.subtract(
                l2,
                Location.empty
            ),
            l2
        );
    }
}
