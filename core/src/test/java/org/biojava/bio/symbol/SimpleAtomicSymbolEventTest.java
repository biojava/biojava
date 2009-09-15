package org.biojava.bio.symbol;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.DNATools;
import org.biojava.utils.ListTools;

/**
 * @author Matthew Pocock
 */
public class SimpleAtomicSymbolEventTest
        extends AbstractSymbolEventTest
{
  protected Symbol createSymbol(Annotation ann)
          throws Exception
  {
    return new SimpleAtomicSymbol(ann, new ListTools.Doublet(DNATools.c(), DNATools.t()));
  }
}
