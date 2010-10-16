package org.biojava.naming;

import java.util.Properties;

import javax.naming.Binding;
import javax.naming.Context;
import javax.naming.NamingEnumeration;
import javax.naming.directory.Attributes;
import javax.naming.directory.DirContext;
import javax.naming.directory.InitialDirContext;

import junit.framework.TestCase;

/**
 *
 *
 * @author Matthew Pocock
 */
public class ObdaInitialContextFactoryTest
        extends TestCase
{
  public void testLookup()
          throws Exception
  {
    Properties env = new Properties();
    env.put(Context.INITIAL_CONTEXT_FACTORY,
            "org.biojava.naming.ObdaInitialContextFactory");
    DirContext context = new InitialDirContext(env);

    DirContext embl = (DirContext) context.lookup("urn:open-bio.org:format:embl");
    assertNotNull("Fetched embl name", embl);

    Attributes desc = embl.getAttributes("", new String[] { "description" });
    assertEquals("Got one description attribute", 1, desc.size());
    assertNotNull("Description is not null", desc.get("description"));
  }

  public void testWalk()
          throws Exception
  {
    Properties env = new Properties();
    env.put(Context.INITIAL_CONTEXT_FACTORY,
            "org.biojava.naming.ObdaInitialContextFactory");
    DirContext context = new InitialDirContext(env);

    walk(context);
  }

  private void walk(DirContext context)
          throws Exception
  {
    System.out.println("Reached " + context.getNameInNamespace());
    System.out.println("  Attributes: " + context.getAttributes(""));

    NamingEnumeration ne = context.listBindings("");
    while(ne.hasMore()) {
      Binding b = (Binding) ne.nextElement();
      System.out.println("Binding: " + b.getName() + " -> " + b.getObject());
      walk((DirContext) b.getObject());
    }
  }
}