package org.biojava.bio.symbol;

import java.io.InvalidObjectException;
import java.io.ObjectStreamException;
import java.io.Serializable;
import java.util.List;

import org.biojava.bio.Annotation;

/**
 * Provides custom serialization hooks.
 *
 * @since 1.4
 * @author Matthew Pocock
 */
class SimpleBasisSymbol
        extends AbstractSimpleBasisSymbol
        implements Serializable
{
  public SimpleBasisSymbol(Annotation annotation, List symbols)
          throws IllegalSymbolException
  {
    super(annotation, symbols);
  }

  public SimpleBasisSymbol(Annotation annotation)
  {
    super(annotation);
  }

  public SimpleBasisSymbol(Annotation annotation, Alphabet matches)
  {
    super(annotation, matches);
  }

  public SimpleBasisSymbol(Annotation annotation, List symbols, Alphabet matches)
          throws IllegalSymbolException
  {
    super(annotation, symbols, matches);
  }

  private Object writeReplace()
  {
    return new SBSH(this);
  }

  private static class SBSH
          implements Serializable
  {
    private List syms;
    private Annotation ann;

    public SBSH(SimpleBasisSymbol sym)
    {
      syms = sym.getSymbols();
      ann = sym.getAnnotation();
    }

    public Object readResolve()
            throws ObjectStreamException
    {
      try {
        return AlphabetManager.createSymbol(ann, syms, null);
      } catch (IllegalSymbolException ex) {
        throw new InvalidObjectException("Couldn't resolve symbol:" + syms);
      }
    }
  }
}
