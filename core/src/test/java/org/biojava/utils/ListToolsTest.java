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

package org.biojava.utils;

import junit.framework.TestCase;

/**
 * @author Matthew Pocock
 */
public class ListToolsTest
extends TestCase {
  public ListToolsTest(String name) {
    super(name);
  }
  
  public void testDoublet() {
    ListTools.Doublet doublet = new ListTools.Doublet("foo", "bar");
    assertTrue(doublet.size() == 2);

    assertTrue(doublet.getA() == "foo");
    assertTrue(doublet.getB() == "bar");
    
    assertTrue(doublet.get(0) == "foo");
    assertTrue(doublet.get(1) == "bar");
  }
  
  public void testTriplet() {
    ListTools.Triplet triplet = new ListTools.Triplet("foo", "bar", "baz");
    assertTrue(triplet.size() == 3);

    assertTrue(triplet.getA() == "foo");
    assertTrue(triplet.getB() == "bar");
    assertTrue(triplet.getC() == "baz");
   
    assertTrue(triplet.get(0) == "foo");
    assertTrue(triplet.get(1) == "bar");
    assertTrue(triplet.get(2) == "baz");
  }
  
}
