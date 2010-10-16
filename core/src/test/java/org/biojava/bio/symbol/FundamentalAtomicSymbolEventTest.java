package org.biojava.bio.symbol;

import org.biojava.bio.Annotation;
import org.biojava.utils.ChangeSupport;

/**
 * @author Matthew Pocock
 */
public class FundamentalAtomicSymbolEventTest extends AbstractSymbolEventTest
{
  protected Symbol createSymbol(Annotation ann)
  {
    return new FundamentalAtomicSymbol("me", ann)
    {
      protected ChangeSupport generateChangeSupport()
      {
        return new ChangeSupport();
      }
    };
  }
}
