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
 * <code>PointLocationTest</code> tests the behaviour of
 * <code>PointLocation</code> by itself and combined with
 * <code>LocationTools</code>.
 *
 * @author <a href="mailto:kdj@sanger.ac.uk">Keith James</a>
 * @author Francois Pepin
 * @since 1.2
 */
public class PointLocationTest extends TestCase
{
    protected Location r1;
    protected Location r2;
    protected Location r3;
    protected Location r4;
  
    public PointLocationTest(String name)
    {
	super(name);
    }

  /**
   * Runs the unit tests defined here.
   */
  public static void main(String args[])
  {
    junit.textui.TestRunner.run(PointLocationTest.class);
  }
    protected void setUp() throws Exception
    {
	r1 = new PointLocation(1);
	r2 = new PointLocation(10);
         r3 = new PointLocation(11);
         r4 = LocationTools.makeLocation(10,11);
         
    }

    /**
     * <code>testEquals</code> tests equality directly.
     *
     */
    public void testEquals()
    {
  	assertEquals(r1, r1);
  	assertEquals(r1, new PointLocation(1));
    }

    /**
     * <code>testAreEqual</code> tests equality via
     * <code>LocationTools</code>.
     *
     */
    public void testAreEqual()
    {
	assertTrue(LocationTools.areEqual(r1, r1));
	assertTrue(LocationTools.areEqual(r1, new PointLocation(1)));
    }

    /**
     * <code>testOverlaps</code> tests overlaps via
     * <code>LocationTools</code>.
     *
     */
    public void testOverlaps()
    {
	assertTrue(LocationTools.overlaps(r1, r1));
	assertTrue(! LocationTools.overlaps(r1, r2));
    }

    /**
     * <code>testContains</code> tests contains via
     * <code>LocationTools</code>.
     *
     */
    public void testContains()
    {
	assertTrue(LocationTools.contains(r1, new PointLocation(1)));
	assertTrue(! LocationTools.contains(r1, r2));

	assertTrue(r1.contains(1));
	assertTrue(! r1.contains(2));
    }

    /**
     * <code>testIntersection</code> tests intersection via
     * <code>LocationTools</code>.
     *
     */
    public void testIntersection()
    {
	assertEquals(LocationTools.intersection(r1, r1),
		     new PointLocation(1));
	assertEquals(LocationTools.intersection(r1, r2),
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
        assertEquals(r4, LocationTools.union(r2,r3));
   	assertTrue(LocationTools.union(r2,r3).isContiguous());     
    }
    /**
     * <code>testIsContiguous</code> tests contiguous.
     *
     */
    public void testIsContiguous()
    {
	assertTrue(r1.isContiguous());
    } 
}
