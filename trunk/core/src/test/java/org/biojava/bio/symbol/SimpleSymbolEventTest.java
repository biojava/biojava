package org.biojava.bio.symbol;

import org.biojava.bio.Annotation;


/**
 * @author Matthew Pocock
 */
public class SimpleSymbolEventTest extends AbstractSymbolEventTest
{
  protected Symbol createSymbol(Annotation ann)
          throws Exception
  {
    return new SimpleSymbol(ann);
  }
}

