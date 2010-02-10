/**
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

package org.biojava.bio.seq.io;

import java.io.NotSerializableException;
import java.io.ObjectStreamException;
import java.io.Serializable;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.PackedSymbolListFactory;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.StaticMemberPlaceHolder;

/**
 * Basic SequenceBuilder implementation which accumulates all
 * notified information and chooses a sequence implementation
 * suited to the size of the sequence. This may or may not bit-encode
 * the symbols (using PackedSymbolList), and may or may not store the
 * symbols in multiple fixed-length chunks (using ChunkedSymbolList).
 *
 * <p>More functionality is offered by {@link org.biojavax.bio.seq.io.SimpleRichSequenceBuilder SimpleRichSequenceBuilder},
 * Use of this class is prefered.</p>
 *
 * @author David Huen
 * @author Matthew Pocock
 */

public class SmartSequenceBuilder extends SequenceBuilderBase {
  public final static SequenceBuilderFactory FACTORY = new SSBFactory(-1);
  public final static SequenceBuilderFactory BIT_PACKED = new SSBFactory(0);

  private static class SSBFactory implements SequenceBuilderFactory, Serializable {
    private final int threshold;

    private SSBFactory(int threshold) {
      this.threshold = threshold;
    }

    public SequenceBuilder makeSequenceBuilder() {
      return new SmartSequenceBuilder(threshold);
    }

    private Object writeReplace() throws ObjectStreamException {
      try {
        return new StaticMemberPlaceHolder(SimpleSequenceBuilder.class.getField("FACTORY"));
      } catch (NoSuchFieldException nsfe) {
        throw new NotSerializableException(nsfe.getMessage());
      }
    }
  }

  private ChunkedSymbolListFactory slFactory;

  private SmartSequenceBuilder(int threshold) {
    slFactory = new ChunkedSymbolListFactory(
            new PackedSymbolListFactory(), threshold);
  }

  //
  // SeqIOListener
  //

  public void addSymbols(Alphabet alpha, Symbol[] syms, int pos, int len)
          throws IllegalAlphabetException
  {
    slFactory.addSymbols(alpha, syms, pos, len);
  }


  public Sequence makeSequence()
          throws BioException
  {
    SymbolList symbols;
    try {
      symbols = slFactory.makeSymbolList();
      seq = new SimpleSequence(symbols, uri, name, annotation);
    }
    catch (IllegalAlphabetException iae) {
      // this shouldn't happen!!!
    }

    return super.makeSequence();
  }
}
