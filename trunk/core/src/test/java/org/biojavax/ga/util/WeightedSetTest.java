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
package org.biojavax.ga.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;

import junit.framework.TestCase;


/**
 * @author Mark Schreiber
 */
public class WeightedSetTest
    extends TestCase {
  private WeightedSet weightedSet = null;

  public WeightedSetTest(String name) {
    super(name);
  }

  protected void setUp() throws Exception {
    super.setUp();

    weightedSet = new WeightedSet();
  }

  protected void tearDown() throws Exception {
    weightedSet = null;
    super.tearDown();
  }


  public void testAdd() {
    Object o = new Object();
    boolean expectedReturn = true;
    boolean actualReturn = weightedSet.add(o);
    assertEquals("return value", expectedReturn, actualReturn);
    assertTrue(weightedSet.getWeight(o) == 0.0);
    assertTrue(weightedSet.size() == 1);

    actualReturn = weightedSet.add(o);
    expectedReturn = false;
    assertEquals("return value", expectedReturn, actualReturn);
    assertTrue(weightedSet.size() == 1);
  }

  public void testAddAll() {
    List c = new ArrayList();
    Object o = new Object();
    String s = "";

    c.add(o); c.add(s);

    boolean expectedReturn = true;
    boolean actualReturn = weightedSet.addAll(c);
    assertEquals("return value", expectedReturn, actualReturn);
    assertTrue(weightedSet.size() == 2);
    assertTrue(weightedSet.getWeight(o) == 0.0);
    assertTrue(weightedSet.getWeight(s) == 0.0);
  }

  public void testAsMap() {
    weightedSet.add("one");

    Map m = weightedSet.asMap();
    assertTrue(m.containsKey("one"));
    Double expectedReturn = new Double(0.0);
    Double actualReturn = (Double)m.get("one");
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testClear() {
    weightedSet.setWeight("one", 0.5);
    weightedSet.setWeight("two", 0.5);
    weightedSet.clear();
    assertTrue(weightedSet.getTotalWeight() == 0.0);
    assertTrue(weightedSet.contains("one") == false);
    assertTrue(weightedSet.contains("two") == false);
  }

  public void testContains() {
    Object o = new Object();
    boolean expectedReturn = false;
    boolean actualReturn = weightedSet.contains(o);
    assertEquals("return value", expectedReturn, actualReturn);

    weightedSet.add(o);
    expectedReturn = true;
    actualReturn = weightedSet.contains(o);
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testContainsAll() {
    List c = new ArrayList();
    Object o = new Object();
    String s = "";

    c.add(o);
    c.add(s);


    boolean expectedReturn = false;
    boolean actualReturn = weightedSet.containsAll(c);
    assertEquals("return value", expectedReturn, actualReturn);

    weightedSet.addAll(c);
    expectedReturn = true;
    actualReturn = weightedSet.containsAll(c);
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testGetTotalWeight() {
    Double expectedReturn = new Double(1.5);
    weightedSet.setWeight("one", 0.5);
    weightedSet.setWeight("two", 0.5);
    weightedSet.setWeight("three", 0.5);
    Double actualReturn = new Double(weightedSet.getTotalWeight());
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testGetWeight() throws NoSuchElementException {
    weightedSet.setWeight("one", 0.5);
    weightedSet.setWeight("two", 0.5);
    weightedSet.setWeight("three", 0.5);
    weightedSet.setWeight("four", 0.5);

    Double expectedReturn = new Double(0.25);
    Double actualReturn = new Double(weightedSet.getWeight("one"));
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testIsEmpty() {
    weightedSet.setWeight("four", 0.5);
    boolean expectedReturn = false;
    boolean actualReturn = weightedSet.isEmpty();
    assertEquals("return value", expectedReturn, actualReturn);
  }


  public void testRemove() {
    weightedSet.setWeight("one", 0.5);
    weightedSet.setWeight("two", 0.5);
    weightedSet.setWeight("three", 0.5);
    weightedSet.setWeight("four", 0.5);
    weightedSet.setWeight("five", 0.5);

    weightedSet.remove("one");
    assertTrue(weightedSet.getTotalWeight() == 2.0);
    assertTrue(weightedSet.getWeight("five") == 0.25);
  }

  public void testRetainAll() {
    weightedSet.setWeight("one", 0.5);
    weightedSet.setWeight("two", 0.5);
    weightedSet.setWeight("three", 0.5);
    weightedSet.setWeight("four", 0.5);
    weightedSet.setWeight("five", 0.5);

    Collection c = new ArrayList();
    c.add("one"); c.add("two");

    boolean expectedReturn = true;
    boolean actualReturn = weightedSet.retainAll(c);
    assertEquals("return value", expectedReturn, actualReturn);

    assertTrue(weightedSet.contains("one"));
    assertTrue(weightedSet.containsAll(c));
    assertTrue(! weightedSet.contains("three"));
  }

  public void testSample() {
    weightedSet.setWeight("one", 0.5);
    Object expectedReturn = "one";
    Object actualReturn = weightedSet.sample();
    assertEquals("return value", expectedReturn, actualReturn);

    weightedSet.setWeight("two", 4.5);
  }

  public void testSetWeight() {
    Object o = "one";
    double w = 0.3;
    weightedSet.setWeight(o, w);

    assertTrue(weightedSet.getTotalWeight() == 0.3);
    assertTrue(weightedSet.getWeight(o) == 1.0);
  }

  public void testSetWeight2(){
    Object o = "one";
    double w = 2.5;
    weightedSet.setWeight(o, w);

    assertTrue(weightedSet.getTotalWeight() == 2.5);
    assertTrue(weightedSet.getWeight(o) == 1.0);
  }

  public void testSetWeight3(){
    Object o = "one";
    double w = 2.5;
    Object p = "two";
    double x = 2.5;

    weightedSet.setWeight(o, w);
    assertTrue(weightedSet.getTotalWeight() == 2.5);
    assertTrue(weightedSet.getWeight(o) == 1.0);
    weightedSet.setWeight(p, x);
    assertTrue(weightedSet.getTotalWeight() == 5.0);
    assertTrue(weightedSet.getWeight(o) == 0.5);
    assertTrue(weightedSet.getWeight(p) == 0.5);
  }

  public void testSize() {
    int expectedReturn = 0;
    int actualReturn = weightedSet.size();
    assertEquals("return value", expectedReturn, actualReturn);

    weightedSet.setWeight("one", 0.5);
    weightedSet.setWeight("two", 0.5);
    weightedSet.setWeight("three", 0.5);
    weightedSet.setWeight("four", 0.5);
    weightedSet.setWeight("five", 0.5);

    expectedReturn = 5;
    actualReturn = weightedSet.size();
    assertEquals("return value", expectedReturn, actualReturn);

  }
}