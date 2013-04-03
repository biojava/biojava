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


package org.biojava.bio.dp;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.alignment.SimpleAlignment;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;

/**
 * A no-frills implementation of StatePath.
 *
 * @author Matthew Pocock
 * @author Nimesh Singh
 */
public class SimpleStatePath
  extends
    Unchangeable
  implements
    StatePath,
    Serializable
{
  private final double score;
  private final Alignment delegate;

  public double getScore() {
    return score;
  }

  public SimpleStatePath(
    double score,
    SymbolList sequence,
    SymbolList states,
    SymbolList scores
  ) throws IllegalArgumentException {
    this.score = score;
    Map map = new HashMap();
    map.put(StatePath.SEQUENCE, sequence);
    map.put(StatePath.STATES, states);
    map.put(StatePath.SCORES, scores);
    this.delegate = new SimpleAlignment(map);
  }

  public Alphabet getAlphabet() {
    return delegate.getAlphabet();
  }

  public List getLabels() {
    return delegate.getLabels();
  }

  public int length() {
    return delegate.length();
  }

  public Alignment subAlignment(Set labels, Location loc)
  throws NoSuchElementException {
    return delegate.subAlignment(labels, loc);
  }

  public Symbol symbolAt(int col)
  throws IndexOutOfBoundsException {
    return delegate.symbolAt(col);
  }

  public Symbol symbolAt(String label, int col)
  throws IndexOutOfBoundsException, NoSuchElementException {
    return delegate.symbolAt(label, col);
  }

  public SymbolList symbolListForLabel(String label)
  throws NoSuchElementException {
    return delegate.symbolListForLabel(label);
  }

  public Iterator<Symbol> iterator() {
    return delegate.iterator();
  }

  public SymbolList subList(int start, int end) {
    return delegate.subList(start, end);
  }

  public List<Symbol> toList() {
    return delegate.toList();
  }

  public String seqString() {
    return delegate.seqString();
  }

  public String subStr(int start, int end)
  throws IndexOutOfBoundsException {
    return delegate.subStr(start, end);
  }

  public void edit(Edit edit)
  throws IllegalAlphabetException, IndexOutOfBoundsException, ChangeVetoException {
    throw new ChangeVetoException("Can't edit SimpleStatePath");
  }

  public Iterator<SymbolList> symbolListIterator() {
    return new Alignment.SymbolListIterator(this);
  }
}
