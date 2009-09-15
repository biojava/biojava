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
 * single unicode characters.</p>
 *
 * <p>Many alphabets (and all simple built-in alphabets like DNA, RNA
 * and Protein) will have an instance of CharacterTokenization
 * registered under the name 'token', so that you could say
 * <code>CharacterTokenization ct = (CharacterTokenization)
 * alpha.getTokenization('token');</code> and expect it to work. When
 * you construct a new instance of this class for an alphabet, there
 * will be no initial associations of Symbols with characters. It is
 * your responsibility to populate the new tokenization appropriately.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Greg Cox
 * @author Keith James
 * @since 1.2
 */

public class CharacterTokenization
  extends
    Unchangeable
  implements
    SymbolTokenization, Serializable
{
    private Alphabet alphabet;
    private Map symbolsToCharacters = new HashMap();
    private Map charactersToSymbols = new HashMap();
    private transient Symbol[] tokenTable;
    private boolean caseSensitive;

    public CharacterTokenization(Alphabet alpha, boolean caseSensitive) {
        alphabet = alpha;
        this.caseSensitive = caseSensitive;
    }

    public Alphabet getAlphabet() {
        return alphabet;
    }

    public TokenType getTokenType() {
        return CHARACTER;
    }

    public Annotation getAnnotation() {
        return Annotation.EMPTY_ANNOTATION;
    }

    /**
     * <p>
     * Bind a Symbol to a character.
     * </p>
     *
     * <p>
     * This method will ensure that when this char is observed, it resolves to
     * this symbol. If it was previously associated with another symbol, the old
     * binding is removed.
     * If this is the first time the symbol has been bound to any character,
     * then this character is taken to be the default tokenization of the
     * Symbol. This means that when converting symbols into characters, this
     * char will be used. If the symbol has previously been bound to another
     * character, then this char will not be produced for the symbol when
     * stringifying the symbol, but this symbol will be produced when tokenizing
     * this character.
     * </p>
     *
     * @param s  the Symbol to bind
     * @param c  the char to bind it to
     */
    public void bindSymbol(Symbol s, char c) {
        Character chr = new Character(c);

        if (!symbolsToCharacters.containsKey(s)) {
            symbolsToCharacters.put(s, chr);
        }
        if (!charactersToSymbols.containsKey(chr)) {
            charactersToSymbols.put(chr, s);
        }
        tokenTable = null;
    }

    public Symbol parseToken(String token)
        throws IllegalSymbolException
    {
        if (token.length() != 1) {
            throw new IllegalSymbolException("This Tokenization only accepts single-character tokens");
        }
        return parseTokenChar(token.charAt(0));
    }

    protected Symbol[] getTokenTable() {
        if (tokenTable == null) {
            int maxChar = 0;
            for (Iterator i = charactersToSymbols.keySet().iterator(); i.hasNext(); ) {
                Character c = (Character) i.next();
                char cv = c.charValue();
                if (caseSensitive) {
                    maxChar = Math.max(maxChar, cv);
                } else {
                    maxChar = Math.max(maxChar, Character.toUpperCase(cv));
                    maxChar = Math.max(maxChar, Character.toLowerCase(cv));
                }
            }

            tokenTable = new Symbol[maxChar + 1];

            for (Iterator i = charactersToSymbols.entrySet().iterator(); i.hasNext(); ) {
                Map.Entry me = (Map.Entry) i.next();
                Symbol sym = (Symbol) me.getValue();
                Character c = (Character) me.getKey();
                char cv = c.charValue();
                if (caseSensitive) {
                    tokenTable[cv] = sym;
                } else {
                    tokenTable[Character.toUpperCase(cv)] = sym;
                    tokenTable[Character.toLowerCase(cv)] = sym;
                }
            }
        }

        return tokenTable;
    }

    protected Symbol parseTokenChar(char c)
        throws IllegalSymbolException
    {
        Symbol[] tokenTable = getTokenTable();
        Symbol sym = null;
        if (c < tokenTable.length) {
            sym = tokenTable[c];
        }
        if (sym == null) {
            throw new IllegalSymbolException("This tokenization doesn't contain character: '" + c + "'");
        }

        return sym;
    }

    private Character _tokenizeSymbol(Symbol s)
        throws IllegalSymbolException
    {
        Character c = (Character) symbolsToCharacters.get(s);
        if (c == null) {
            Alphabet alpha = getAlphabet();
            alphabet.validate(s);
            if (alpha instanceof FiniteAlphabet) {
                c = (Character) symbolsToCharacters.get(AlphabetManager.getAllAmbiguitySymbol((FiniteAlphabet) alpha));
            }
            if (c == null) {
                throw new IllegalSymbolException("No mapping for symbol " + s.getName());
            }
        }

        return c;
    }

    public String tokenizeSymbol(Symbol s) throws IllegalSymbolException {
        return String.valueOf(_tokenizeSymbol(s).charValue());
    }

    public String tokenizeSymbolList(SymbolList sl)
        throws IllegalAlphabetException
    {
        if (sl.getAlphabet() != getAlphabet()) {
            throw new IllegalAlphabetException("Alphabet " + sl.getAlphabet().getName() + " does not match " + getAlphabet().getName());
        }
        StringBuffer sb = new StringBuffer();
        for (Iterator i = sl.iterator(); i.hasNext(); ) {
            Symbol sym = (Symbol) i.next();
            try {
                Character c = _tokenizeSymbol(sym);
                sb.append(c.charValue());
            } catch (IllegalSymbolException ex) {
                throw new IllegalAlphabetException(ex, "Couldn't tokenize");
            }
        }

        return sb.substring(0);
    }

    public StreamParser parseStream(SeqIOListener listener) {
        return new TPStreamParser(listener);
    }

    private class TPStreamParser implements StreamParser {
        private SeqIOListener listener;
        private Symbol[] buffer;

        {
            buffer = new Symbol[256];
        }

        public TPStreamParser(SeqIOListener l) {
            this.listener = l;
        }

        public void characters(char[] data, int start, int len)
            throws IllegalSymbolException
        {
            int cnt = 0;
            while (cnt < len) {
                int bcnt = 0;
                while (cnt < len && bcnt < buffer.length) {
                    buffer[bcnt++] = parseTokenChar(data[start + (cnt++)]);
                }
                try {
                    listener.addSymbols(getAlphabet(),
                                        buffer,
                                        0,
                                        bcnt);
                } catch (IllegalAlphabetException ex) {
                    throw new BioError( "Assertion failed: can't add symbols.", ex);
                }
            }
        }

        public void close() {
        }
    }
}
