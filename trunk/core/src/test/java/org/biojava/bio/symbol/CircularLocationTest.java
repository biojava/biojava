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
import java.util.Iterator;

import junit.framework.TestCase;

/**
 * <code>CircularLocationTest</code> tests the behaviour of
 * <code>CircularLocation</code> by itself and combined with
 * <code>CircularLocationTools</code>.
 *
 * @author <a href="mailto:mark.schreiber@agresearch.co.nz">Mark Schreiber</a>
 * @author Francois Pepin
 * @since 1.3
 */
public class CircularLocationTest extends TestCase
{
    protected CircularLocation r1;
    protected CircularLocation r2;
    protected CircularLocation r3;
    protected CircularLocation r4;
    protected CircularLocation r5;
    protected CircularLocation r6;
    protected CircularLocation r7;
    protected CircularLocation r8;
    protected CircularLocation r9;


    public CircularLocationTest(String name)
    {
        super(name);
    }

  /**
   * Runs the unit tests defined here.
   */
  public static void main(String args[])
  {
    junit.textui.TestRunner.run(CircularLocationTest.class);
  }

    protected void setUp() throws Exception
    {
        r1 = LocationTools.makeCircularLocation(1, 100, 200);
        r2 = LocationTools.makeCircularLocation(1, 200, 200);
        r3 = LocationTools.makeCircularLocation(105, 300, 200);
        r4 = LocationTools.makeCircularLocation(-50, 20, 200);
        r5 = LocationTools.makeCircularLocation(90, 110, 200);

        r6 = LocationTools.makeCircularLocation(18,24,20);
        r7 = LocationTools.makeCircularLocation(2,8,20);
        r8 = LocationTools.makeCircularLocation(2,7,20);
        r9 = LocationTools.makeCircularLocation(8,8,20);
    }

     /**
     * <code>testConstruction</code> tests construction of Locations.
     *
     */
    public void testConstruction()
    {
        assertTrue(r1 != null);
        assertTrue(r2 != null);
        assertTrue(r3 != null);
        assertTrue(r4 != null);
        assertTrue(r5 != null);
    }
    /**
     * <code>testEquals</code> tests equality directly.
     *
     */
    public void testEquals()
    {
        assertEquals(r1, r1);
        assertEquals(r1, LocationTools.makeCircularLocation(1, 100, 200));
    }

    public void testToString(){
      assertEquals(r1.toString(), "[1,100] (circular)");
      assertEquals(r2.toString(), "[1,200] (circular)");
      assertEquals(r3.toString(), "105, 100 {([105,200]), ([1,100])}  (circular)");
      assertEquals(r4.toString(), "150, 20 {([150,200]), ([1,20])}  (circular)");
      assertEquals(r5.toString(), "[90,110] (circular)");
      assertEquals(r6.toString(), "18, 4 {([18,20]), ([1,4])}  (circular)");
      assertEquals(r7.toString(),"[2,8] (circular)");
    }

    /**
     * <code>testAreEqual</code> tests equality via
     * <code>LocationTools</code>.
     *
     */
    public void testAreEqual()
    {
        assertTrue(LocationTools.areEqual(r1, r1));
        assertTrue(LocationTools.areEqual(r1, LocationTools.makeCircularLocation(1, 100, 200)));
    }

    /**
     * <code>testOverlaps</code> tests overlaps via
     * <code>LocationTools</code>.
     *
     */
    public void testOverlaps()
    {
        assertTrue(LocationTools.overlaps(r1, r1));
        assertTrue(LocationTools.overlaps(r2, r2));
        assertTrue(LocationTools.overlaps(r1, r2));
        assertTrue(LocationTools.overlaps(r2, r1));
        assertTrue(LocationTools.overlaps(r4, r2));
        assertTrue(LocationTools.overlaps(r2, r4));
        assertTrue(LocationTools.overlaps(r5, r3));
        assertTrue(LocationTools.overlaps(r3 ,r5));
        assertTrue(LocationTools.overlaps(r6 ,r7));
        assertTrue(LocationTools.overlaps(r7 ,r6));

        assertTrue(! LocationTools.overlaps(r5, r4));
    }

    /**
     * <code>testContains</code> tests contains via
     * <code>LocationTools</code>.
     *
     */
    public void testContains()
    {
        assertTrue(LocationTools.contains(r2, r1));
        assertTrue(LocationTools.contains(r2, r3));
        assertTrue(LocationTools.contains(r2, r4));
        assertTrue(LocationTools.contains(r2, r5));
        assertTrue(LocationTools.contains(r3, r4));

        assertTrue(! LocationTools.contains(r4, r5));
        assertTrue(! LocationTools.contains(r4, r3));
        assertTrue(! LocationTools.contains(r1, r2));

        assertTrue(r1.contains(1));
        assertTrue(r1.contains(100));
        assertTrue(r1.contains(202));
        assertTrue(! r1.contains(101));
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
        assertEquals(LocationTools.intersection(r5, r4),
                     Location.empty);
        assertEquals(LocationTools.intersection(r4, r5),
                     Location.empty);
    }

    /**
     * <code>testUnion</code> tests union via
     * <code>LocationTools</code>.
     *
     */
    public void testUnion()
    {
        assertEquals(r1, LocationTools.union(r1, r1));
        assertEquals(LocationTools.union(r1, r2),
                     LocationTools.union(r2, r1));

        CircularLocation[] locs = new CircularLocation[3];
        locs[0] = LocationTools.makeCircularLocation(13,14,20);
        locs[1] = LocationTools.makeCircularLocation(18,4,20);
        locs[2] = (CircularLocation)LocationTools.union(locs[0], locs[1]);

        assertTrue(locs[2].get5PrimeEnd() == 13);
        assertEquals(locs[2].toString(),
                     "13, 4 {([13,14]), ([18,20]), ([1,4])}  (circular)");


        //test a more complex union
        ArrayList al = new ArrayList();
        al.add(LocationTools.makeCircularLocation(5,8,20));
        al.add(LocationTools.makeCircularLocation(7,16,20));
        al.add(LocationTools.makeCircularLocation(19,2, 20));

        CircularLocation loc = (CircularLocation)LocationTools.union(al);
        assertTrue(loc.get5PrimeEnd() == 5);
        assertTrue(loc.overlapsOrigin());


        //check the components
        Iterator blocki = loc.blockIterator();
        Location locA = (Location)blocki.next();
        assertTrue(locA.getMin() == 1);
        assertTrue(locA.getMax() == 2);
        Location locB = (Location)blocki.next();
        assertTrue(locB.getMin() == 5);
        assertTrue(locB.getMax() == 16);
        Location locC = (Location)blocki.next();
        assertTrue(locC.getMin() == 19);
        assertTrue(locC.getMax() == 20);
        //shouldn't be more
        assertTrue(blocki.hasNext() == false);
        //test unions when locations are touching
        assertEquals(r7, LocationTools.union(r8,r9));
    }

    /**
     * <code>testIsContiguous</code> tests contiguous.
     *
     */
    public void testIsContiguous()
    {
        assertTrue(r1.isContiguous());
        assertTrue(r2.isContiguous());
        assertTrue(! r3.isContiguous());
        assertTrue(! r4.isContiguous());
        assertTrue(r5.isContiguous());

        Location l = LocationTools.union(new RangeLocation(1,4),new RangeLocation(7,10));
        CircularLocation cl = new CircularLocation(l,200);
        assertTrue(! cl.isContiguous());
    }
}
