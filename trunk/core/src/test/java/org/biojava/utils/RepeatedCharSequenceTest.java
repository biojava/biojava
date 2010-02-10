package org.biojava.utils;

import junit.framework.TestCase;

/**
 *
 *
 * @author Matthew Pocock
 */
public class RepeatedCharSequenceTest
        extends TestCase
{
  public void testDefaultConstructor()
  {
    RepeatedCharSequence rcs = new RepeatedCharSequence();
    String rcsS = rcs.toString();

    assertEquals("Empty String", "", rcsS);
  }

  public void testCustomConstructor()
  {
    RepeatedCharSequence rcs = new RepeatedCharSequence(5, 'x');
    String rcsS = rcs.toString();

    assertEquals("Five x", "xxxxx", rcsS);
  }


  public void testLengthIncrease()
  {
    RepeatedCharSequence rcs = new RepeatedCharSequence(4, '?');
    String rcsS1 = rcs.toString();
    rcs.setLength(8);
    String rcsS2 = rcs.toString();

    assertEquals("Original length 4", 4, rcsS1.length());
    assertEquals("Current length 8", 8, rcsS2.length());
    assertEquals("Eight ?", "????????", rcsS2);
  }

  public void testLengthDecrease()
  {
    RepeatedCharSequence rcs = new RepeatedCharSequence(8, '*');
    String rcsS1 = rcs.toString();
    rcs.setLength(4);
    String rcsS2 = rcs.toString();

    assertEquals("Original length 8", 8, rcsS1.length());
    assertEquals("Current length 4", 4, rcsS2.length());
    assertEquals("Four *", "****", rcsS2);
  }

  public void testCharChange()
  {
    RepeatedCharSequence rcs = new RepeatedCharSequence(4, '$');
    String rcsS1 = rcs.toString();
    rcs.setCharacter('@');
    String rcsS2 = rcs.toString();
System.out.println("JAM");
    assertEquals("Length not changed", rcsS1.length(), rcsS2.length());
    assertEquals("Now a string of @", "@@@@", rcsS2);
  }

  public void testSubSequence()
  {
    RepeatedCharSequence rcs = new RepeatedCharSequence(1000, 'k');
    CharSequence rcs2 = rcs.subSequence(100, 110);
    CharSequence rcs3 = rcs2.subSequence(3, 3);

    assertEquals("10 ks", "kkkkkkkkkk", rcs2.toString());
    assertEquals("empty", "", rcs3.toString());
  }
}