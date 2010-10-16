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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.alignment.Alignment;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;

/**
 * An alignment that relabels another alignment.
 *
 * @author Matthew Pocock
 * @author Nimesh Singh
 */
public class RelabeledAlignment
  extends
    Unchangeable
  implements
    Alignment
{
  private Alignment delegate;
  private Map<String, String> labelMap = new HashMap<String, String>();

  public RelabeledAlignment(Alignment delegate) {
    this.delegate = delegate;
    for(Iterator<String> i = delegate.getLabels().iterator(); i.hasNext(); ) {
      String label = i.next();
      labelMap.put(label, label);
    }
  }

  public List<String> getLabels() {
    return new ArrayList<String>(labelMap.keySet());
  }

  public Symbol symbolAt(String label, int column)
  throws NoSuchElementException {
    return delegate.symbolAt(map(label), column);
  }

  public SymbolList symbolListForLabel(String label)
  throws NoSuchElementException {
    return delegate.symbolListForLabel(map(label));
  }

  public Alignment subAlignment(Set<String> labels, Location loc)
  throws NoSuchElementException {
    return new RelabeledAlignment(delegate.subAlignment(map(labels), loc));
  }

  public String seqString() {
    return delegate.seqString();
  }

  public String subStr(int min, int max) {
    return delegate.subStr(min, max);
  }

  public Alphabet getAlphabet() {
    return delegate.getAlphabet();
  }

  public Iterator iterator() {
    return delegate.iterator();
  }

  public int length() {
    return delegate.length();
  }

  public List toList() {
    return delegate.toList();
  }

  public SymbolList subList(int min, int max) {
    return delegate.subList(min, max);
  }

  public Symbol symbolAt(int pos) {
    return delegate.symbolAt(pos);
  }

  public void edit(Edit edit)
  throws IllegalAlphabetException, ChangeVetoException {
    delegate.edit(edit);
  }

  public Iterator symbolListIterator() {
    return new Alignment.SymbolListIterator(this);
  }
  
  protected Set<String> map(Set<String> labels) {
    Set<String> set = new HashSet<String>();
    for(Iterator<String> i = labels.iterator(); i.hasNext(); ) {
      set.add(map(i.next()));
    }
    return set;
  }

  protected String map(Object label) {
    return labelMap.get(label);
  }
}

