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

import java.io.NotSerializableException;
import java.io.ObjectStreamException;
import java.io.Serializable;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.StaticMemberPlaceHolder;
import org.biojava.utils.Unchangeable;

/**
 * @author Matthew Pocock
 */
class EmptyAlphabet
  extends
    Unchangeable
  implements
    FiniteAlphabet,
    Serializable
{
    public String getName() {
	return "Empty Alphabet";
    }

    public Annotation getAnnotation() {
	return Annotation.EMPTY_ANNOTATION;
    }

    public boolean contains(Symbol s) {
	return s == AlphabetManager.getGapSymbol();
    }

    public void validate(Symbol sym) throws IllegalSymbolException {
	throw new IllegalSymbolException(
					 "The empty alphabet does not contain symbol " + sym.getName());
    }

    public SymbolTokenization getTokenization(String name) throws NoSuchElementException {
      throw new NoSuchElementException("There is no parser for the empty alphabet. Attempted to retrieve " + name);
    }

    public int size() {
      return 0;
    }

    public List getAlphabets() {
	return Collections.EMPTY_LIST;
    }

    public Symbol getSymbol(List syms) throws IllegalSymbolException {
	if(syms.size() != 0) {
	    throw new IllegalSymbolException("The empty alphabet contains nothing");
	}
	return AlphabetManager.getGapSymbol();
    }

    public Symbol getAmbiguity(Set syms) throws IllegalSymbolException {
	for(Iterator i = syms.iterator(); i.hasNext(); ) {
	    this.validate((Symbol) i.next());
	}
	return AlphabetManager.getGapSymbol();
    }

    public Symbol getGapSymbol()
    {
        return AlphabetManager.getGapSymbol();
    }

    public Iterator iterator() {
	return SymbolList.EMPTY_LIST.iterator();
    }

    public void addSymbol(Symbol sym) throws IllegalSymbolException {
      throw new IllegalSymbolException(
        "Can't add symbols to alphabet: " + sym.getName() +
        " in " + getName()
      );
    }

    public void removeSymbol(Symbol sym) throws IllegalSymbolException {
      throw new IllegalSymbolException(
        "Can't remove symbols from alphabet: " + sym.getName() +
        " in " + getName()
      );
    }

    private Object writeReplace() throws ObjectStreamException {
      try {
        return new StaticMemberPlaceHolder(Alphabet.class.getField("EMPTY_ALPHABET"));
      } catch (NoSuchFieldException nsfe) {
        throw new NotSerializableException(nsfe.getMessage());
      }
    }
}

