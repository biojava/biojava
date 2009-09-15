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

import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.seq.DNATools;
import org.biojava.utils.AssertionFailure;
import org.biojava.utils.SmallSet;

/**
 * Packing utility class for DNA.  Also represents
 * ambiguity.
 *
 * @author Matthew Pocock
 * @author David Huen (bugfix)
 */  
public class DNAAmbPack
  implements
    Packing,
    java.io.Serializable
{
  private final Symbol[] syms;
  
  public DNAAmbPack() {
    this.syms = new Symbol[16];
    for(byte i = 0; i < 16; i++) {
      syms[i] = _unpack(i);
    }
  }
  
  public FiniteAlphabet getAlphabet() {
    return DNATools.getDNA();
  }
  
  public byte pack(Symbol sym) {
    if(false) {
    } else if(sym == DNATools.a()) {
      return 1;
    } else if(sym == DNATools.g()) {
      return 2;
    } else if(sym == DNATools.c()) {
      return 4;
    } else if(sym == DNATools.t()) {
      return 8;
    } else if(sym == DNATools.n()) {
      return 15;
    }
    
    byte p = 0;
    for(Iterator i = ((FiniteAlphabet) sym.getMatches()).iterator(); i.hasNext(); ) {
      p |= pack((AtomicSymbol) i.next());
    }
    return p;
  }
  
  public Symbol unpack(byte b) {
    return syms[b];
  }
  
  private Symbol _unpack(byte b) {
    Set syms = new SmallSet();
    if(0 != (b & 1)) {
      syms.add(DNATools.a());
    }
    if(0 != (b & 2)) {
      syms.add(DNATools.g());
    }
    if(0 != (b & 4)) {
      syms.add(DNATools.c());
    }
    if(0 != (b & 8)) {
      syms.add(DNATools.t());
    }
    try {
      // do incorporate the gap symbol when appropriate
      if (b != 0)
          return DNATools.getDNA().getAmbiguity(syms);
      else
          return DNATools.getDNA().getGapSymbol();

    } catch (IllegalSymbolException ise) {
      throw new AssertionFailure("Assertion failure: couldn't get DNA ambiguity from DNA: " + syms, ise);
    }
  }
  
  public byte wordSize() {
    return 4;
  }
  
  public boolean handlesAmbiguity() {
    return true;
  }
}

