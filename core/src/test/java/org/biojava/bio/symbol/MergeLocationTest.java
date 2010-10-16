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
 */


package org.biojava.bio.symbol;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import junit.framework.TestCase;

/**
 * <p>Title: MergeLocationTest</p>
 * <p>Description: Tests the MergeLocation class</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: AgResearch</p>
 * @author Mark Schreiber
 * @author Francois Pepin
 * @version 1.0
 */

public class MergeLocationTest extends TestCase {
  Location a;
  Location b;
  Location c;
  Location d;
  Location e;
  Location f;
  Location g;
  MergeLocation ml;

  List la;
  List lb;

  List nestl;

  Location abc;
  Location abcd;
  MergeLocation fg;
  Location nested;



  public MergeLocationTest(String name){
    super(name);
  }

  /**
   * Runs the unit tests defined here.
   */
  public static void main(String args[])
  {
    junit.textui.TestRunner.run(MergeLocationTest.class);
  }

  protected void setUp() throws java.lang.Exception {
    super.setUp();

    a = new RangeLocation(1, 10);
    b = new RangeLocation(30,40);
    c = new RangeLocation(35,65);
    d = new RangeLocation(80,90);
    e = new RangeLocation(60,70);
    f = new RangeLocation(60,64);
    g = new RangeLocation(65,70);

    la = new ArrayList();
    lb = new ArrayList();
    nestl = new ArrayList();

    la.add(a);
    la.add(b);
    la.add(c);

    lb.add(a);
    lb.add(b);
    lb.add(c);
    lb.add(d);

    nestl.add(a);
    nestl.add(d);
    nestl.add(e);

    ml = MergeLocation.mergeLocations(b,c);

    nestl.add(ml);

    abc = LocationTools.union(la);
    abcd = LocationTools.union(lb);
    fg = MergeLocation.mergeLocations(f,g);
    nested = LocationTools.union(nestl);

  }
  protected void tearDown() throws java.lang.Exception {
    super.tearDown();
  }

  public void testInstantiation()throws Exception{
    assertNotNull(ml);
    assertNotNull(abc);
    assertNotNull(abcd);
    assertNotNull(nested);
  }

  public void testBlocks(){
    assertNotNull(ml.getComponentList(false));

    assertTrue(fg.getComponentList(false).size() == 2);
    assertTrue(ml.getComponentList(false).size() == 2);
    assertTrue(ml.getComponentList(true).size() == 2);
  }

  public void testTouchingMerge(){
    assertTrue(fg.getComponentList(true).size() == 2);
    assertTrue(fg.isContiguous());
  }


  public void testCompoundABC(){
    int blocks = 0;
    MergeLocation loc = null;

    for (Iterator i = abc.blockIterator(); i.hasNext(); ) {
      Object item = i.next();
      if(item instanceof MergeLocation){
        loc = (MergeLocation)item;
      }
      blocks++;
    }

    assertTrue(blocks == 2);
    assertNotNull(loc);
    assertTrue(abc.getMax() == 65);
    assertTrue(abc.getMin() == 1);

    assertTrue(loc.getComponentList(false).size() == 2);
    assertTrue(loc.getComponentList(true).size() == 2);
  }

  public void testCompoundABCD(){
    int blocks = 0;
    MergeLocation loc = null;

    for (Iterator i = abcd.blockIterator(); i.hasNext(); ) {
      Object item = i.next();
      if(item instanceof MergeLocation){
        loc = (MergeLocation)item;
      }
      blocks++;
    }
    assertTrue(blocks == 3);
    assertTrue(abcd.getMin() == 1);
    assertTrue(abcd.getMax() == 90);

    assertTrue(loc.getComponentList(true).size() == 2);
    assertTrue(loc.getComponentList(false).size() == 2);
  }

  public void testNested(){
    int blocks = 0;
    MergeLocation loc = null;

    for (Iterator i = nested.blockIterator(); i.hasNext(); ) {
      Object item = i.next();
      if(item instanceof MergeLocation){
        loc = (MergeLocation)item;
      }
      blocks++;
    }

    assertTrue(blocks == 3);
    assertTrue(nested.getMax() == 90);
    assertTrue(nested.getMin()==1);

    assertTrue(loc.getComponentList(false).size()==2);
    assertTrue(loc.getComponentList(true).size()==3);
  }
}
