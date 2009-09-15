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

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import org.biojava3.core.symbol.SymbolCoder.BooleanSymbolCoder;
import org.biojava3.core.symbol.SymbolCoder.CharacterSymbolCoder;
import org.biojava3.core.symbol.SymbolCoder.DefaultSymbolCoder;
import org.biojava3.core.symbol.SymbolCoder.DoubleSymbolCoder;
import org.biojava3.core.symbol.SymbolCoder.IntegerSymbolCoder;
import org.biojava3.core.symbol.SymbolCoder.StringSymbolCoder;

/**
 * Tools for working with alphabets.
 * @author Richard Holland
 * @since 3.0
 */
public final class AlphabetTools {

    /**
     * Given a collection of symbols, construct an alphabet using them as 
     * strings. If the alphabet already exists, the existing alphabet will be
     * returned with the additional symbols added to it. Any existing symbols
     * will be retained.
     * @param name the name of the alphabet.
     * @param syms the symbols to use.
     * @return the alphabet.
     */
    public static Alphabet createSimpleAlphabet(String name, Collection<Symbol> syms) {
        return createSimpleAlphabet(name, syms, new DefaultSymbolCoder());
    }

    /**
     * Given a collection of symbols, construct an alphabet using them as 
     * strings. If the alphabet already exists, the existing alphabet will be
     * returned with the additional symbols added to it. Any existing symbols
     * will be retained.
     * @param name the name of the alphabet.
     * @param syms the symbols to use.
     * @param coder the coder to use.
     * @return the alphabet.
     */
    public static Alphabet createSimpleAlphabet(String name, Collection<Symbol> syms,
            SymbolCoder coder) {
        Alphabet alphabet = new Alphabet(name);
        alphabet.addAll(syms);
        alphabet.setCoder(coder);
        return alphabet;
    }

    /**
     * Given a string of single-character symbol names, return a collection of
     * symbols.
     * @param syms the symbol string to use one character at a time.
     * @param useBothCases {@code true} if both upper- and lower-case versions
     * of each character should be added to the collection.
     */
    public static Collection<Symbol> quickSymbols(String str, boolean useBothCases) {
        Collection<Symbol> syms = new HashSet<Symbol>();
        if (useBothCases) {
            str = str.toLowerCase() + str.toUpperCase();
        }
        for (int i = 0; i < str.length(); i++) {
            syms.add(Symbol.get(str.charAt(i)));
        }
        return syms;
    }
    /**
     * The standard gap symbol which all alphabets use by default.
     */
    public static final Symbol STANDARD_GAP = Symbol.get('-');
    /**
     * Use this alphabet to represent sets of characters. It is empty by default.
     */
    public static final Alphabet CHARACTER_ALPHABET = createSimpleAlphabet("CHARACTER", new ArrayList<Symbol>(), new CharacterSymbolCoder());
    /**
     * Use this alphabet to represent sets of integers. It is empty by default.
     */
    public static final Alphabet INTEGER_ALPHABET = createSimpleAlphabet("INTEGER", new ArrayList<Symbol>(), new IntegerSymbolCoder());
    /**
     * Use this alphabet to represent sets of strings. It is empty by default.
     */
    public static final Alphabet STRING_ALPHABET = createSimpleAlphabet("STRING", new ArrayList<Symbol>(), new StringSymbolCoder());
    /**
     * Use this alphabet to represent sets of doubles. It is empty by default.
     */
    public static final Alphabet DOUBLE_ALPHABET = createSimpleAlphabet("DOUBLE", new ArrayList<Symbol>(), new DoubleSymbolCoder());
    /**
     * Use this alphabet to represent sets of booleans. It is empty by default.
     */
    public static final Alphabet BOOLEAN_ALPHABET = createSimpleAlphabet("BOOLEAN", new ArrayList<Symbol>(), new BooleanSymbolCoder());
}
