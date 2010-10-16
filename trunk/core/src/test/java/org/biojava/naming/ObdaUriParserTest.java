package org.biojava.naming;

import javax.naming.Name;
import javax.naming.NamingException;

import junit.framework.TestCase;

/**
 *
 *
 * @author Matthew Pocock
 */
public class ObdaUriParserTest
        extends TestCase
{
  public void testEmpty()
          throws NamingException
  {
    Name name = ObdaUriParser.getInstance().parse("");
    System.out.println(name);
    assertEquals("Empty string parses to no names", 0, name.size());
  }

  public void testSingle()
          throws NamingException
  {
    Name name = ObdaUriParser.getInstance().parse("oneName");
    System.out.println(name);
    assertEquals("Single part parses to single name", 1, name.size());
  }

  public void testLeading()
          throws NamingException
  {
    Name name = ObdaUriParser.getInstance().parse(":trail");
    System.out.println(name);
    assertEquals("leading part parses double name", 2, name.size());
    assertEquals("leader is empty", "", name.get(0));
    assertEquals("trailer is there", "trail", name.get(1));
  }

  public void testTrailing()
          throws NamingException
  {
    Name name = ObdaUriParser.getInstance().parse("lead:");
    System.out.println(name);
    assertEquals("trailing part parses double name", 2, name.size());
    assertEquals("leader is there", "lead", name.get(0));
    assertEquals("trailer is empty", "", name.get(1));
  }

  public void testLots()
          throws NamingException
  {
    Name name = ObdaUriParser.getInstance().parse("urn:obda.org:format:embl/ac");
    System.out.println(name);
    assertEquals("Splits into 4 parts", 4, name.size());
    assertEquals("part0=urn", "urn", name.get(0));
    assertEquals("part1=obda.org", "obda.org", name.get(1));
    assertEquals("part2=format", "format", name.get(2));
    assertEquals("part3=embl/ac", "embl/ac", name.get(3));
  }
}