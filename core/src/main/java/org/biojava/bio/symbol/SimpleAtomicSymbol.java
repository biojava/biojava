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


package org.biojava.bio.symbol;

import java.io.InvalidObjectException;
import java.io.ObjectStreamException;
import java.io.Serializable;
import java.util.List;

import org.biojava.bio.Annotation;

/**
 * A basic implementation of AtomicSymbol.
 *
 * If you wish to construct new Symbols, you should normally do so via utility methods
 * on <code>AlphabetManager</code>.
 *
 * This may be a useful base class for custom implementations.
 *
 * @author Matthew Pocock
 */
public class SimpleAtomicSymbol
        extends AbstractSimpleBasisSymbol
        implements AtomicSymbol, Serializable
{
  protected SimpleAtomicSymbol(
    Annotation annotation, List syms
  ) throws IllegalSymbolException {
    super(annotation, syms);
  }

  protected Alphabet createMatches() {
    return new SingletonAlphabet(this);
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

    public SBSH(SimpleAtomicSymbol sym)
    {
      syms = sym.getSymbols();
      ann = sym.getAnnotation();
    }

    public Object readResolve()
            throws ObjectStreamException
    {
      try {
        Symbol sym = AlphabetManager.createSymbol(ann, syms, null);
        return sym;
      } catch (IllegalSymbolException ex) {
        throw new InvalidObjectException("Couldn't resolve symbol:" + syms);
      }
    }
  }
}
