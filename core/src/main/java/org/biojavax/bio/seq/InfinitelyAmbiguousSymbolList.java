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

package org.biojavax.bio.seq;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * A symbol list that is <code>Integer.MAX_VALUE</code>long, never gives index out of
 * bounds and always returns ambiguity symbols for everything.
 * @author Richard Holland
 * @author MarkSchreiber
 * @since 1.5
 */
public class InfinitelyAmbiguousSymbolList implements SymbolList {
    
    private FiniteAlphabet fa;
    private Symbol sym;
    
    /**
     * 
     * Creates a new instance of InfinitelyAmbiguousSymbolList 
     * @param fa the finite alphabet to return ambiguous symbols from.
     */
    public InfinitelyAmbiguousSymbolList(FiniteAlphabet fa) {
        this.fa = fa;
        this.sym = AlphabetManager.getAllAmbiguitySymbol(fa);
    }

    /**
     * {@inheritDoc}
     * IGNORED
     */
    public void addChangeListener(ChangeListener cl) {}

    /**
     * {@inheritDoc}
     * IGNORED
     */
    public void addChangeListener(ChangeListener cl, ChangeType ct) {}

    /**
     * {@inheritDoc}
     * IGNORED
     */
    public void edit(Edit edit) throws IndexOutOfBoundsException, IllegalAlphabetException, ChangeVetoException {}

    /**
     * {@inheritDoc}
     */
    public Alphabet getAlphabet() { return this.fa; }

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS TRUE
     */
    public boolean isUnchanging(ChangeType ct) { return true; }

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS AN ITERATOR OVER A SINGLE AMBIGUITY SYMBOL
     */
    public Iterator iterator() { return Collections.singletonList(this.sym).iterator(); }

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS <code>Integer.MAX_VALUE</code>
     */
    public int length() { return Integer.MAX_VALUE; }

    /**
     * {@inheritDoc}
     * IGNORED
     */
    public void removeChangeListener(ChangeListener cl) {}

    /**
     * {@inheritDoc}
     * IGNORED
     */
    public void removeChangeListener(ChangeListener cl, ChangeType ct) {}

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE AMBIGUITY SYMBOL REPRESENTED AS A STRING
     */
    public String seqString() { return this.sym.toString(); }

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS SELF
     */
    public SymbolList subList(int start, int end) throws IndexOutOfBoundsException { return this; }

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE CORRECT LENGTH STRING MADE UP OF AMBIGUITY SYMBOLS
     */
    public String subStr(int start, int end) throws IndexOutOfBoundsException {
        int min = Math.min(start,end);
        int max = Math.max(start,end);
        StringBuffer sb = new StringBuffer();
        String symStr = this.sym.toString();
        while (max-->=min) sb.append(symStr);
        return sb.toString();
    }

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS THE AMBIGUITY SYMBOL
     */
    public Symbol symbolAt(int index) throws IndexOutOfBoundsException { return this.sym; }

    /**
     * {@inheritDoc}
     * ALWAYS RETURNS A LIST CONTAINING THE AMBIGUITY SYMBOL
     */
    public List toList() { return Collections.singletonList(this.sym); }
    
}
