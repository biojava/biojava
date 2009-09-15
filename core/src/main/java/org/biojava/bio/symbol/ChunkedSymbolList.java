package org.biojava.bio.symbol;

import java.io.IOException;
import java.io.Serializable;

import org.biojava.utils.ChangeListener;

/**
 * SymbolList implementation using constant-size chunks. Each chunk provides
 * the same number of symbols (except the last one, which may be shorter). When
 * a request for symbols comes in, firstly the apropreate chunk is located, and
 * then the symbols are extracted. This implementation is used in the IO package
 * to avoid allocating and re-allocating memory when the total length of the
 * symbol list is unknown. It can also be used when chunks are to be lazily
 * fetched from some high-latency stoorage by implementing a single lazy-fetch
 * SymbolList class and populating a ChunkedSymbolList with a complete
 * tile-path of them.
 *
 * @author David Huen
 * @author Matthew Pocock
 */
public class ChunkedSymbolList
        extends AbstractSymbolList
        implements Serializable
{
  // state
  private SymbolList [] chunks;
  private final int chunkSize;
  private final Alphabet alpha;
  private final int length;

  // cached info for speedups
  private transient int currentMin = Integer.MAX_VALUE;
  private transient int currentMax = Integer.MIN_VALUE;
  private transient SymbolList currentChunk = null;

  private void readObject(java.io.ObjectInputStream stream)
          throws IOException, ClassNotFoundException
  {
    stream.defaultReadObject();

    currentMin = Integer.MAX_VALUE;
    currentMax = Integer.MIN_VALUE;
    currentChunk = null;
  }

  protected void finalize() throws Throwable {
    super.finalize();
    alpha.removeChangeListener(ChangeListener.ALWAYS_VETO, Alphabet.SYMBOLS);
  }

  public ChunkedSymbolList(SymbolList [] chunks,
                           int chunkSize,
                           int length,
                           Alphabet alpha)
  {
    this.chunks = chunks;
    this.chunkSize = chunkSize;
    this.length = length;
    this.alpha = alpha;
    alpha.addChangeListener(ChangeListener.ALWAYS_VETO, Alphabet.SYMBOLS);
  }

  public Alphabet getAlphabet() {
    return alpha;
  }

  public int length() {
    return length;
  }

  public Symbol symbolAt(int pos) {
    int offset;
    --pos;

    synchronized(this) {
      if ((pos < currentMin) || (pos > currentMax)) {
        // bad - we are outside the current chunk
        int chnk = pos / chunkSize;
        offset =  pos % chunkSize;

        currentMin = pos - offset;
        currentMax = currentMin + chunkSize - 1;
        currentChunk = chunks[chnk];
      }
      else {
        offset = pos - currentMin;
      }
    }

    return currentChunk.symbolAt(offset + 1);
  }

  public SymbolList subList(int start, int end) {
    if (start < 1 || end > length()) {
      throw new IndexOutOfBoundsException(
              "Sublist index out of bounds " + length() + ":" + start + "," + end
      );
    }

    if (end < start) {
      throw new IllegalArgumentException(
              "end must not be lower than start: start=" + start + ", end=" + end
      );
    }

    //
    // Mildly optimized for case where from and to are within
    // the same chunk.
    //

    int afrom = start - 1;
    int ato = end - 1;
    int cfrom = afrom / chunkSize;
    if (ato / chunkSize == cfrom) {
      return chunks[cfrom].subList((afrom % chunkSize) + 1, (ato % chunkSize) + 1);
    } else {
      return super.subList(start, end);
    }
  }
}
