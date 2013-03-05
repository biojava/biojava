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
 * Created on 5 Mar 2013
 * Created by Andreas Prlic
 *
 * @since 3.0.6
 */
package org.biojava.bio.structure.math;


import java.io.Serializable;
import java.util.Iterator;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 *  Sorted symbol table implementation using a java.util.TreeMap.
 *  Does not allow duplicate keys.
 * 
 *  This class represents an ordered symbol table. It assumes that
 *  the elements are <tt>Comparable</tt>.
 *  It supports the usual <em>put</em>, <em>get</em>, <em>contains</em>,
 *  and <em>delete</em> methods.
 *  It also provides ordered methods for finding the <em>minimum</em>,
 *  <em>maximum</em>, <em>floor</em>, and <em>ceiling</em>.
 *  <p>
 *  The class uses the convention that values cannot be null. Setting the
 *  value associated with a key to null is equivalent to removing the key.
 *  <p>
 *  This implementation uses a balanced binary search tree.
 *  The <em>add</em>, <em>contains</em>, <em>delete</em>, <em>minimum</em>,
 *  <em>maximum</em>, <em>ceiling</em>, and <em>floor</em> methods take
 *  logarithmic time.
 *  
 *  Derived from http://introcs.cs.princeton.edu/java/44st/ST.java.html
 *  
 *  <p>
 *  For additional documentation, see <a href="http://introcs.cs.princeton.edu/44st">Section 4.4</a> of
 *  <i>Introduction to Programming in Java: An Interdisciplinary Approach</i> by Robert Sedgewick and Kevin Wayne. 
 *
 */
public class SymbolTable<Key extends Comparable<Key>, Value> implements Iterable<Key>, Serializable {
	
    /**
	 * 
	 */
	private static final long serialVersionUID = -4417561575046471931L;
	
	private TreeMap<Key, Value> st;

    /**
     * Create an empty symbol table.
     */
    public SymbolTable() {
        st = new TreeMap<Key, Value>();
    }

    /**
     * Put key-value pair into the symbol table. Remove key from table if
     * value is null.
     */
    public void put(Key key, Value val) {
        if (val == null) st.remove(key);
        else             st.put(key, val);
    }

    /**
     * Return the value paired with given key; null if key is not in table.
     */
    public Value get(Key key) {
        return st.get(key);
    }

    /**
     * Delete the key (and paired value) from table.
     * Return the value paired with given key; null if key is not in table.
     */
    public Value delete(Key key) {
        return st.remove(key);
    }

    /**
     * Is the key in the table?
     */
    public boolean contains(Key key) {
        return st.containsKey(key);
    }

    /**
     * How many keys are in the table?
     */
    public int size() {
        return st.size();
    }

    /**
     * Return an <tt>Iterator</tt> for the keys in the table.
     * To iterate over all of the keys in the symbol table <tt>st</tt>, use the
     * foreach notation: <tt>for (Key key : st)</tt>.
     */ 
    public Iterator<Key> iterator() {
        return st.keySet().iterator();
    }

    /**
     * Return an <tt>Iterable</tt> for the keys in the table.
     * To iterate over all of the keys in the symbol table <tt>st</tt>, use the
     * foreach notation: <tt>for (Key key : st.keys())</tt>.
     */ 
    public Iterable<Key> keys() {
        return st.keySet();
    }

    /**
     * Return the smallest key in the table.
     */ 
    public Key min() {
        return st.firstKey();
    }

    /**
     * Return the largest key in the table.
     */ 
    public Key max() {
        return st.lastKey();
    }

    /**
     * Return the smallest key in the table >= k.
     */ 
    public Key ceil(Key k) {
        SortedMap<Key, Value> tail = st.tailMap(k);
        if (tail.isEmpty()) return null;
        else return tail.firstKey();
    }

    /**
     * Return the largest key in the table <= k.
     */ 
    public Key floor(Key k) {
        if (st.containsKey(k)) return k;

        // does not include key if present (!)
        SortedMap<Key, Value> head = st.headMap(k);
        if (head.isEmpty()) return null;
        else return head.lastKey();
    }

}
