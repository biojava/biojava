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
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.Annotation;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.Unchangeable;

/**
 * Base class for tokenizations which accept whitespace-separated
 * `words'.  Splits at whitespace, except when it is quoted by
 * either double-quotes ("), brackets (), or square brackets [].
 *
 * @author Thomas Down
 * @author Greg Cox
 * @author Keith James
 * @since 1.2
 */

public abstract class WordTokenization
  extends
    Unchangeable
  implements
    SymbolTokenization, Serializable
{
    private Alphabet alphabet;

    public WordTokenization(Alphabet fab) {
	this.alphabet = fab;
    }

    public Alphabet getAlphabet() {
	return alphabet;
    }

    public TokenType getTokenType() {
	return SEPARATED;
    }

    public Annotation getAnnotation() {
	return Annotation.EMPTY_ANNOTATION;
    }

    public String tokenizeSymbolList(SymbolList sl)
        throws IllegalSymbolException, IllegalAlphabetException
    {
	if (sl.getAlphabet() != getAlphabet()) {
	    throw new IllegalAlphabetException("Alphabet " + sl.getAlphabet().getName() + " does not match " + getAlphabet().getName());
	}
	StringBuffer sb = new StringBuffer();
	Iterator i = sl.iterator();
	while (i.hasNext()) {
	    Symbol sym = (Symbol) i.next();
	    sb.append(tokenizeSymbol(sym));
	    if (i.hasNext()) {
		sb.append(' ');
	    }
	}
	return sb.substring(0);
    }

    public StreamParser parseStream(SeqIOListener siol) {
	return new WordStreamParser(siol);
    }

    protected List splitString(String str)
	throws IllegalSymbolException
    {
	int ptr = 0;
	List sl = new ArrayList();

	while (ptr < str.length()) {
	    char c = str.charAt(ptr);
	    if (Character.isWhitespace(c)) {
		++ptr;
	    } else if (c == '(') {
		int nextPtr = findMatch(str, ptr, '(', ')');
		sl.add(str.substring(ptr, nextPtr));
		ptr = nextPtr;
	    } else if (c == '[') {
		int nextPtr = findMatch(str, ptr, '[', ']');
		sl.add(str.substring(ptr, nextPtr));
		ptr = nextPtr;
	    } else {
		int nextPtr = ptr;
		char nc;
		boolean quoted = false;
		do {
		    nextPtr++;
		    if (nextPtr == str.length()) {
			nc = ' ';
		    } else {
			nc = str.charAt(nextPtr);
		    }
		    if (nc == '"') {
			quoted = !quoted;
		    }
		} while (!Character.isWhitespace(nc));

		sl.add(str.substring(ptr, nextPtr));
		ptr = nextPtr;
	    }
	}

	return sl;
    }

    protected Symbol[] parseString(String s)
        throws IllegalSymbolException
    {
	List split = splitString(s);
	Symbol[] syms = new Symbol[split.size()];
	for (int i = 0; i < split.size(); ++i) {
	    syms[i] = parseToken((String) split.get(i));
	}
	return syms;
    }

    private class WordStreamParser implements StreamParser {
	SeqIOListener listener;
	StringBuffer sb = new StringBuffer();

	WordStreamParser(SeqIOListener l) {
	    listener = l;
	}

	public void characters(char[] data, int start, int len) {
	    sb.append(data, start, len);
	}

	public void close()
	    throws IllegalSymbolException
	{
	    String str = sb.substring(0);
	    Symbol[] syms = parseString(str);
	    try {
		listener.addSymbols(alphabet, syms, 0, syms.length);
	    } catch (IllegalAlphabetException ex) {
		throw new IllegalSymbolException("Mismatched alphabets");
	    }
	}
    }

    private int findMatch(String str, int ptr, char openChar, char closeChar) {
	int level = 0;
	do {
	    char c = str.charAt(ptr++);
	    if (c == openChar) {
		++level;
	    } else if (c == closeChar) {
		--level;
	    }
	} while (level > 0);
	return ptr;
    }
}
