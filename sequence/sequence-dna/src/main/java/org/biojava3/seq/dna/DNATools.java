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
package org.biojava3.seq.dna;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.biojava3.core.symbol.Alphabet;
import org.biojava3.core.symbol.AlphabetTools;
import org.biojava3.core.symbol.Symbol;
import org.biojava3.core.symbol.SymbolCoder.CharacterSymbolCoder;

public final class DNATools {

    /**
     * This is the DNA alphabet, prepopulated with both upper- and lower-case
     * versions of the DNA symbols A, T, G and C.
     */
    public static final Alphabet DNA_ALPHABET = AlphabetTools.createSimpleAlphabet("DNA", AlphabetTools.quickSymbols("ACGTURYKMSWBDHVN", true),
            new CharacterSymbolCoder());
    private static final Map<Symbol, Symbol> complementMap = new HashMap<Symbol, Symbol>();

    static {
        addComplement("A", "T");
        addComplement("C", "G");
        addComplement("R", "Y");
        addComplement("K", "M");
        addComplement("B", "V");
        addComplement("D", "H");
        
        addAmbiguity("R", "GA");
        addAmbiguity("Y", "TC");
        addAmbiguity("K", "GT");
        addAmbiguity("S", "AC");
        addAmbiguity("W", "GC");
        addAmbiguity("B", "AT");
        addAmbiguity("D", "GTC");
        addAmbiguity("H", "GAT");
        addAmbiguity("V", "GCA");
        addAmbiguity("N", "AGCT");
    }
    
    /**
     * Shortcut private method for setting up the complement map 
     * in both directions and both cases.
     * @param a symbol a.
     * @param b symbol b.
     */
    private static void addComplement(String a, String b) {
        char lca = a.toLowerCase().charAt(0);
        char uca = a.toUpperCase().charAt(0);
        char lcb = b.toLowerCase().charAt(0);
        char ucb = b.toUpperCase().charAt(0);
        complementMap.put(Symbol.get(lca), Symbol.get(lcb));
        complementMap.put(Symbol.get(lcb), Symbol.get(lca));
        complementMap.put(Symbol.get(uca), Symbol.get(ucb));
        complementMap.put(Symbol.get(ucb), Symbol.get(uca));
    }
    
    /**
     * Shortcut private method for setting up the ambiguity map 
     * in both cases.
     * @param a symbol a.
     * @param b symbol b.
     */
    private static void addAmbiguity(String a, String b) {
        char lca = a.toLowerCase().charAt(0);
        char uca = a.toUpperCase().charAt(0);
        String lcbs = b.toLowerCase();
        String ucbs = b.toUpperCase();
        for (int i = 0; i < b.length(); i++) {
            char lcb = lcbs.charAt(i);
            char ucb = ucbs.charAt(i);
            DNA_ALPHABET.addAmbiguity(Symbol.get(lca), Symbol.get(lcb));
            DNA_ALPHABET.addAmbiguity(Symbol.get(uca), Symbol.get(ucb));
        }
    }

    /**
     * Reverse the DNA.
     * @param symList the DNA to be reversed in-place.
     */
    public static void reverse(List<Symbol> symList) {
        Collections.reverse(symList);
    }

    /**
     * Complement the DNA.
     * @param symList the DNA to be complemented in-place. Anything in the list 
     * that isn't DNA will remain untouched.
     */
    public static void complement(List<Symbol> symList) {
        for (int i = 0; i < symList.size(); i++) {
            Symbol sym = symList.get(i);
            Symbol comp = complementMap.get(sym);
            if (comp != null) {
                symList.set(i, comp);
            }
        }
    }

    /**
     * Reverse-complement the DNA.
     * @param symList the DNA to be reverse-complemented in-place. Anything in 
     * the list that isn't DNA will remain untouched.
     */
    public static void reverseComplement(List<Symbol> symList) {
        reverse(symList);
        complement(symList);
    }
    private static final Collection<Symbol> GC = AlphabetTools.quickSymbols("GC", true);

    /**
     * Calculate the % GC content of the DNA. This will count both upper-
     * and lower-case versions of G and C and return their number as a 
     * percentage of the length of the list.
     * @param symList the list to count % GC percentage for.
     * @return the % of the list that is G or C in either upper- or lower-case.
     */
    public static double percentGC(List<Symbol> symList) {
        List<Symbol> workingList = new ArrayList<Symbol>(symList);
        workingList.retainAll(GC);
        return (double) workingList.size() / (double) symList.size();
    }
}
