package org.biojava.bio.symbol;

import java.util.Iterator;

import junit.framework.TestCase;

import org.biojava.bio.seq.ProteinTools;

/**
 * @author Matthew Pocock
 */
public class IndexedNoAmbPackTest
        extends TestCase
{
  public void testSymToIndexToSym()
          throws IllegalAlphabetException, IllegalSymbolException
  {
    FiniteAlphabet alpha = ProteinTools.getAlphabet();
    AlphabetIndex index = AlphabetManager.getAlphabetIndex(alpha);
    Packing packing = new IndexedNoAmbPack(index);

    for(Iterator i = alpha.iterator(); i.hasNext(); ) {
      Symbol sym = (Symbol) i.next();
      byte packed = packing.pack(sym);
      Symbol unpacked = packing.unpack(packed);

      assertEquals("Packing and unpacking round-trips", sym, unpacked);
    }
  }

  public void testIndexToSymToIndex()
          throws IllegalAlphabetException, IllegalSymbolException
  {
    FiniteAlphabet alpha = ProteinTools.getAlphabet();
    AlphabetIndex index = AlphabetManager.getAlphabetIndex(alpha);
    Packing packing = new IndexedNoAmbPack(index);

    for (byte i = 0; i < alpha.size(); i++) {
      Symbol unpacked = packing.unpack(i);
      byte repacked = packing.pack(unpacked);

      assertEquals("Unpacking and packing round-trips", i, repacked);
    }
  }
}