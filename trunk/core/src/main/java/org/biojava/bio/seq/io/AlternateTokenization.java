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


package org.biojava.bio.seq.io;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.Unchangeable;

/**
 * <p>Implementation of SymbolTokenization which binds symbols to
 * strings of characters. These tokenizations are intented to provide
 *  alternate way of writing sequences into Strings.  Therefore they cannot be
 * used for parsing files.</p>
 *
 * <p>As this release, alternate tokenizations are available for the built-in
 * DNA alphabet (write symbols as capital letter) and PROTEIN-TERM alphabet
 * (write symbol as triplets of characters with the first one being a capital
 * letter as in "Glu".</p>
 *
 * <p>By convention, instances of AlternateTokenization should have an associated
 * token starting by the word 'alternate'.
 *
 * @author George Waldon
 * @since 1.5
 */

public class AlternateTokenization extends Unchangeable
        implements SymbolTokenization, Serializable {
    
    private Alphabet alphabet;
    private Map symbolsToStrings = new HashMap();
    private boolean caseSensitive;
    private boolean initiated = false;
    private int width = 0;
    
    public AlternateTokenization(Alphabet alpha, boolean caseSensitive) {
        alphabet = alpha;
        this.caseSensitive = caseSensitive;
    }
    
    public Alphabet getAlphabet() {
        return alphabet;
    }
    
    /** Tokens have fixed size.
     */
    public TokenType getTokenType() {
        return FIXEDWIDTH;
    }
    
    public Annotation getAnnotation() {
        return Annotation.EMPTY_ANNOTATION;
    }
    
    private synchronized void init(String str) {
        if(initiated) return;
        width = str.length();
        initiated = true;
    }
    
    /** Get the width of the tokens.
     */
    public int getWidth() {
        if(initiated==false)
            throw new IllegalStateException("Tokenization not initialize yet");
        return width;
    }
    
    /** Bind a Symbol to a string.
     *
     * @param s  the Symbol to bind
     * @param str  the string to bind it to
     */
    public void bindSymbol(Symbol s, String str) {
        if(!initiated)
            init(str);
        if(str.length()!=width)
            throw new IllegalArgumentException("This tokenization must have all its tokens with the same size");
        if (!symbolsToStrings.containsKey(s)) {
            symbolsToStrings.put(s, str);
        }
    }
    
    /** Will throw an exception.
     */
    public Symbol parseToken(String token)
    throws IllegalSymbolException {
        throw new UnsupportedOperationException("AlternateTokenization are for writing only");
    }
    
    public String tokenizeSymbol(Symbol s) throws IllegalSymbolException {
        String str = (String) symbolsToStrings.get(s);
        if (str == null) {
            Alphabet alpha = getAlphabet();
            alphabet.validate(s);
            if (alpha instanceof FiniteAlphabet) {
                str = (String) symbolsToStrings.get(AlphabetManager.getAllAmbiguitySymbol((FiniteAlphabet) alpha));
            }
            if (str == null) {
                throw new IllegalSymbolException("No mapping for symbol " + s.getName());
            }
        }
        return str;
    }
    
    public String tokenizeSymbolList(SymbolList sl)
    throws IllegalAlphabetException {
        if (sl.getAlphabet() != getAlphabet()) {
            throw new IllegalAlphabetException("Alphabet " + sl.getAlphabet().getName() + " does not match " + getAlphabet().getName());
        }
        StringBuffer sb = new StringBuffer();
        for (Iterator i = sl.iterator(); i.hasNext(); ) {
            Symbol sym = (Symbol) i.next();
            try {
                String str = tokenizeSymbol(sym);
                sb.append(str);
            } catch (IllegalSymbolException ex) {
                throw new IllegalAlphabetException(ex, "Couldn't tokenize");
            }
        }
        
        return sb.substring(0);
    }
    
    /** Will throw an exception.
     */
    public StreamParser parseStream(SeqIOListener listener) {
        throw new UnsupportedOperationException("AlternateTokenization are for writing only");
    }
}
