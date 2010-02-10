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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Suffix tree implementation.  The interface is a bit strange, as it
 * needed to be as space-efficient as possible. More work could be
 * done on the space issue.
 *
 * <p>
 * A suffix tree is an efficient method for encoding the frequencies
 * of motifs in a sequence.  They are sometimes used to quickly screen
 * for similar sequences.  For instance, all motifs of length up to
 * 2 in the sequence <code>AAGT</code> could be encoded as:
 * </p>
 *
 * <pre>
 * root(4)
 * |
 * A(2)--------G(1)-----T(1)
 * |           |
 * A(1)--G(1)  T(1)
 * </pre>
 *
 * <p>
 * A possible method of comparing SuffixTrees is provided as a kernel
 * function as <code>org.biojava.stats.svm.tools.SuffixTreeKernel</code>.
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down (documentation and other updates) 
 */

public class SuffixTree implements Serializable {
  private FiniteAlphabet alphabet;
  private SuffixNode root;
  private AlphabetIndex indexer;
  private List counts;
  
    /**
     * Return the Alphabet containing all Symbols which might be found in
     * this SuffixTree.
     */

  public FiniteAlphabet getAlphabet() {
    return alphabet;
  }

    /**
     * Return the node object which is the root of this suffix tree.
     * This represents the set of all motifs found in this tree.
     */

  public SuffixNode getRoot() {
    return root;
  }
  
    /**
     * Get a child of a SuffixTree.SuffixNode, constructing a new
     * one if need be.  This method is here due to memory optimisations.
     */

  public SuffixNode getChild(SuffixNode node, Symbol s)
  throws IllegalSymbolException {
    if(!getAlphabet().contains(s)) {
      return null;
    }
    int index = indexer.indexForSymbol(s);
    return getChild(node, index);
  }
  
    /**
     * Get the n'th child of a node.
     */

  public SuffixNode getChild(SuffixNode node, int i) {
    if(!node.hasChild(i)) {
      node.addChild(this, i, new SimpleNode(alphabet.size()));
    }
    return node.getChild(i);
  }
  
    /**
     * Add a count for all motifs with length of up to <code>window</code>
     * to this tree.
     *
     * @param sList a SymbolList whose motifs should be added to the
     *              tree.
     * @param window The maximum motif length to count.
     */

  public void addSymbols(SymbolList sList, int window)
  throws IllegalSymbolException {
    SuffixNode [] buf = new SuffixNode[window];
    int [] counts = new int[window];
    for(int i = 0; i < window; i++) {
      buf[i] = getRoot();
    }
    
    for(int p = 1; p <= sList.length(); p++) {
      Symbol s = sList.symbolAt(p);
      buf[p % window] = getRoot();
      for(int i = 0; i < window; i++) {
        int pi = (p + i) % window;
        if(buf[pi] != null) {
          buf[pi] = getChild(buf[pi], s);
          if(buf[pi] != null) {
            counts[i]++;
            buf[pi].setNumber(buf[pi].getNumber() + 1.0f);
          }
        }
      }
    }
    
    for(int i = 0; i < window; i++) {
      incCounts(i+1, counts[i]);
    } 
  }
  
  protected void incCounts(int i, int c) {
    if(i < counts.size()) {
      Integer oldC = (Integer) counts.get(i-1);
      Integer newC = new Integer(oldC.intValue() + c);
      counts.set(i-1, newC);
    } else {
      counts.add(new Integer(c));
    }
  }
  
    /**
     * Return the length of the longest motif represented in this
     * SuffixTree
     */

  public int maxLength() {
    return counts.size();
  }
  
    /**
     * Return the number of motifs of a given length encoded
     * in this SuffixTree.
     */

  public int frequency(int length) {
    return ((Integer) counts.get(length - 1)).intValue();
  }
  
    /**
     * Construct a new SuffixTree to contain motifs over the
     * specified alphabet.
     *
     * @param alphabet The alphabet of this SuffixTree (must be
     *                 finite).
     */

  public SuffixTree(FiniteAlphabet alphabet) {
    this.alphabet = alphabet;
    this.indexer = AlphabetManager.getAlphabetIndex(alphabet);
    this.counts = new ArrayList();
    this.root = new SimpleNode(alphabet.size());
  }
  
  /**
   * A node in the suffix tree.
   * <p>
   * This class is realy stupid & delegates most work off to a SuffixTree so
   * that it is as small (in memory-per-object terms) as possible.
   *
   * @author Matthew Pocock
   */
  public static abstract class SuffixNode implements Serializable {
      /**
       * Determine is this node is terminal (has no children).
       *
       * @return <code>true</code> if and only if this node has no children.
       */

    abstract public boolean isTerminal();
      
      /**
       * Determine if this node has a child corresponding to a given
       * index number.
       */

    abstract public boolean hasChild(int i);

      /**
       * Return a number (usually, but not always, a motif count)
       * associated with this node of the tree.
       */

    abstract public float getNumber();
      
      /**
       * Set the number associated with this node.
       */

    abstract public void setNumber(float n);

    abstract SuffixNode getChild(int i);
    abstract void addChild(SuffixTree tree, int i, SuffixNode n);
  }

  private static class SimpleNode extends SuffixNode {
    private float number = 0.0f;
    private SuffixNode [] child;
    
    private SuffixNode [] childArray(SuffixTree tree) {
      if(child == null)
        child = new SuffixNode[tree.getAlphabet().size()];
      return child;
    }
    
    public boolean isTerminal() {
      return false;
    }
    
    public boolean hasChild(int i) {
      return child != null && child[i] != null;
    }
    
    public float getNumber() {
      return number;
    }
    
    SuffixNode getChild(int i) {
      if(hasChild(i))
        return child[i];
      return null;
    }
    
    void addChild(SuffixTree tree, int i, SuffixNode n) {
      childArray(tree)[i] = n;
    }
    
    public void setNumber(float n) {
      number = n;
    }
    
    SimpleNode(int c) {
      child = new SuffixNode[c];
    }
  }
}
