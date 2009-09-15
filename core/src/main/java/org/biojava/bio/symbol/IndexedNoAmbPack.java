package org.biojava.bio.symbol;

/**
 * Packing that uses the ordering in an AlphabetIndex to map between bytes and
 * symbols.
 *
 * @since 1.4 
 * @author Matthew Pocock
 */
class IndexedNoAmbPack
        implements Packing, java.io.Serializable
{
  private final AlphabetIndex index;
  private final byte wordSize;

  IndexedNoAmbPack(AlphabetIndex index)
          throws IllegalAlphabetException
  {
    this.index = index;

    // calculate the number of bits needed to encode the alphabet without
    // ambiguity
    // ceil(log_2(size))
    int size = (int) Math.ceil(Math.log(index.getAlphabet().size())
                               / Math.log(2));
    if(size > 8) {
      throw new IllegalAlphabetException("Alphabet too big to pack into a byte");
    }

    this.wordSize = (byte) size;
  }

  public AlphabetIndex getIndex()
  {
    return index;
  }

  public FiniteAlphabet getAlphabet()
  {
    return index.getAlphabet();
  }

  public byte wordSize()
  {
    return wordSize;
  }

  public boolean handlesAmbiguity()
  {
    return false;
  }

  public byte pack(Symbol sym)
          throws IllegalSymbolException
  {
    return (byte) index.indexForSymbol(sym);
  }

  public Symbol unpack(byte packed)
          throws IllegalSymbolException
  {
    return index.symbolForIndex(packed);
  }

  public boolean equals(Object o)
  {
    if (this == o) {
      return true;
    }
    if (!(o instanceof IndexedNoAmbPack)) {
      return false;
    }

    final IndexedNoAmbPack indexedNoAmbPack = (IndexedNoAmbPack) o;

    if (!index.equals(indexedNoAmbPack.index)) {
      return false;
    }

    return true;
  }

  public int hashCode()
  {
    return index.hashCode();
  }
}
