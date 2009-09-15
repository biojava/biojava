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
package org.biojava3.core.symbol;

import java.io.Serializable;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import org.biojava3.core.symbol.SymbolCoder.DefaultSymbolCoder;

/**
 * An alphabet is a plain collection of symbol objects, none of which have
 * to share the same object type. The alphabet maintains ambiguity mappings and
 * allows symbols to be compared for equality based on ambiguity or case match.
 * The alphabet also stores a gap symbol which does not appear directly in the
 * set of symbols reported by the alphabet's {@link Collection} methods. The
 * alphabet is key when it comes to writing or reading lists of symbols, by 
 * providing an encoder and decoder to parse and write them. The default encoder
 * and decoder instances used are {@link SymbolEncoder#DEFAULT} and 
 * {@link SymbolDecoder#DEFAULT}. Alphabets may have names.
 * @author Richard Holland
 * @since 3.0
 */
public class Alphabet extends HashSet<Symbol> implements Serializable, Comparable<Alphabet> {

	private static final long serialVersionUID = 1L;
	
	private String name;
    private SymbolCoder coder = new DefaultSymbolCoder();
    private Symbol gapSymbol = AlphabetTools.STANDARD_GAP;
    private final Map<Symbol, Collection<Symbol>> ambiguities = new HashMap<Symbol, Collection<Symbol>>();

    /**
     * Construct a new alphabet, empty, with the given name and the default
     * gap symbol ("-").
     * @param name the name.
     */
    public Alphabet(String name) {
        this.setName(name);
    }

    /**
     * Set the name of this alphabet. It may be {@code null}.
     * @param name the name of this alphabet.
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * Get the name of this alphabet. It may be {@code null}.
     * @return the name.
     */
    public String getName() {
        return this.name;
    }

    /**
     * Set the gap symbol for this alphabet. If {@code null}, then this alphabet 
     * will not support gaps.
     * @param gapSymbol the gap symbol.
     */
    public void setGapSymbol(Symbol gapSymbol) {
        this.gapSymbol = gapSymbol;
    }

    /**
     * Get the gap symbol for this alphabet.
     * @return the gap symbol. If {@code null}, it means this alphabet does not 
     * support gaps.
     */
    public Symbol getGapSymbol() {
        return this.gapSymbol;
    }

    /**
     * Get the coder for this alphabet that will convert strings to/from symbols.
     * @return the decoder.
     */
    public SymbolCoder getCoder() {
        return coder;
    }

    /**
     * Sets the coder for this alphabet. 
     * @param coder the coder to use. If {@code null}, the default one 
     * will be used instead, {@link DefaultSymbolCoder}.
     */
    public void setCoder(SymbolCoder coder) {
        if (coder == null) {
            coder = new DefaultSymbolCoder();
        }
        this.coder = coder;
    }

    /**
     * Makes sure the key exists in the ambiguity map for this symbol.
     * @param sym the symbol.
     */
    private void ensureAmbiguityKey(Symbol sym) {
        if (!this.ambiguities.containsKey(sym)) {
            this.ambiguities.put(sym, new HashSet<Symbol>());
        }
    }

    /**
     * If the ambiguity map for this symbol exists and is empty, it removes
     * the entry.
     * @param sym the symbol to check.
     */
    private void checkAmbiguityKeyNotEmpty(Symbol sym) {
        if (this.ambiguities.get(sym).isEmpty()) {
            this.ambiguities.remove(sym);
        }
    }

    /**
     * Adds an ambiguity entry.
     * @param sym the symbol.
     * @param ambiguousMatch the ambiguous partner.
     */
    private void _addAmbiguity(Symbol sym, Symbol ambiguousMatch) {
        this.ambiguities.get(sym).add(ambiguousMatch);
    }

    /**
     * Removes am ambiguity entry.
     * @param sym the symbol.
     * @param ambiguousMatch the ambiguous partner.
     */
    private void _removeAmbiguity(Symbol sym, Symbol ambiguousMatch) {
        this.ambiguities.get(sym).remove(ambiguousMatch);
    }

    /**
     * Adds ambiguities to a symbol.
     * @param sym the symbol.
     * @param ambiguousMatches all the symbols it could match ambiguously.
     */
    public void addAmbiguity(Symbol sym, Collection<Symbol> ambiguousMatches) {
        this.ensureAmbiguityKey(sym);
        for (Symbol ambiguousMatch : ambiguousMatches) {
            this._addAmbiguity(sym, ambiguousMatch);
        }
    }

    /**
     * Adds ambiguity to a symbol.
     * @param sym the symbol.
     * @param ambiguousMatch the symbol it could match ambiguously.
     */
    public void addAmbiguity(Symbol sym, Symbol ambiguousMatch) {
        this.ensureAmbiguityKey(sym);
        this._addAmbiguity(sym, ambiguousMatch);
    }

    /**
     * Removes ambiguities from a symbol.
     * @param sym the symbol.
     * @param ambiguousMatches all the symbols it no longer matches ambiguously.
     */
    public void removeAmbiguity(Symbol sym, Collection<Symbol> ambiguousMatches) {
        this.ensureAmbiguityKey(sym);
        for (Symbol ambiguousMatch : ambiguousMatches) {
            this._removeAmbiguity(sym, ambiguousMatch);
        }
        this.checkAmbiguityKeyNotEmpty(sym);
    }

    /**
     * Removes ambiguity from a symbol.
     * @param sym the symbol.
     * @param ambiguousMatch the symbol it no longer matches ambiguously.
     */
    public void removeAmbiguity(Symbol sym, Symbol ambiguousMatch) {
        this.ensureAmbiguityKey(sym);
        this._removeAmbiguity(sym, ambiguousMatch);
        this.checkAmbiguityKeyNotEmpty(sym);
    }

    /**
     * Check to see if the symbol has ambiguities.
     * @param sym the symbol.
     * @return {@code true} if it has ambiguities.
     */
    public boolean hasAmbiguities(Symbol sym) {
        return this.ambiguities.containsKey(sym);
    }

    /**
     * Gets all ambiguities for a symbol.
     * @param sym the symbol.
     * @return a collection of symbols that it matches ambiguously.
     */
    public Collection<Symbol> getAmbiguities(Symbol sym) {
        if (!this.hasAmbiguities(sym)) {
            return Collections.emptySet();
        }
        return this.ambiguities.get(sym);
    }

    /**
     * Within the context of this alphabet, check to see how two symbols
     * might match each other.
     * @param a the first symbol.
     * @param b the second symbol.
     * @return the way in which the two symbols match. See 
     * {@link SymbolMatchType} to see what the possible outcomes might be.
     */
    public SymbolMatchType getSymbolMatchType(Symbol a, Symbol b) {
        if (a.equals(b)) {
            return SymbolMatchType.EXACT_MATCH;
        }
        if ((a.equals(this.gapSymbol) || b.equals(this.gapSymbol)) && !a.equals(b)) {
            return SymbolMatchType.GAP_MATCH;
        }
        if (a.toString().equalsIgnoreCase(b.toString())) {
            return SymbolMatchType.EXACT_STRING_MATCH;
        }
        if (a.toString().equalsIgnoreCase(b.toString())) {
            return SymbolMatchType.MIXED_CASE_MATCH;
        }
        for (Map.Entry<Symbol, Collection<Symbol>> entry : ambiguities.entrySet()) {
            Symbol testFor = null;
            boolean mixedCase = false;
            boolean exactString = false;
            if (entry.getKey().equals(a)) {
                testFor = b;
            } else if (entry.getKey().equals(b)) {
                testFor = a;
            } else if (entry.getKey().toString().equals(a.toString())) {
                testFor = b;
                exactString = true;
            } else if (entry.getKey().toString().equals(b.toString())) {
                testFor = a;
                exactString = true;
            } else if (entry.getKey().toString().equalsIgnoreCase(a.toString())) {
                testFor = b;
                mixedCase = true;
            } else if (entry.getKey().toString().equalsIgnoreCase(b.toString())) {
                testFor = a;
                mixedCase = true;
            }
            if (testFor != null) {
                for (Symbol sym : entry.getValue()) {
                    if (sym.equals(testFor)) {
                        if (exactString) {
                            return SymbolMatchType.AMBIGUOUS_STRING_MATCH;
                        }
                        if (mixedCase) {
                            return SymbolMatchType.AMBIGUOUS_MIXED_CASE_MATCH;
                        }
                        return SymbolMatchType.AMBIGUOUS_MATCH;
                    }
                    if (sym.toString().equals(testFor.toString())) {
                        return SymbolMatchType.AMBIGUOUS_STRING_MATCH;
                    }
                    if (sym.toString().equalsIgnoreCase(testFor.toString())) {
                        return SymbolMatchType.AMBIGUOUS_MIXED_CASE_MATCH;
                    }
                }
            }
        }
        return SymbolMatchType.MISMATCH;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (this.getClass() != obj.getClass()) {
            return false;
        }
        final Alphabet other = (Alphabet) obj;
        return this.name.equals(other.name);
    }

    @Override
    public int hashCode() {
        return this.name.hashCode();
    }

    public int compareTo(Alphabet o) {
        return this.name.compareTo(o.name);
    }

    @Override
    public String toString() {
        return this.name;
    }
}
