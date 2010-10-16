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

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;

/**
 * Simple implementation of SymbolTokenization which uses the `name'
 * field of the symbols.  This class works with any non-cross-product
 * FiniteAlphabet, and doesn't need any extra data to be provided.
 *
 * @author Thomas Down
 * @since 1.2 
 */

public class NameTokenization extends WordTokenization {
    private transient Map nameToSymbol = null;
    private boolean caseSensitive;

    public NameTokenization(FiniteAlphabet fab, boolean caseSensitive) {
	super(fab);
	fab.addChangeListener(ChangeListener.ALWAYS_VETO, ChangeType.UNKNOWN);
	this.caseSensitive = caseSensitive;
    }

    /**
     * Construct a new NameTokenization, defaulting to case-insensitive.
     */

    public NameTokenization(FiniteAlphabet fab) {
	this(fab, false);
    }

    protected void finalize() throws Throwable {
	super.finalize();
	getAlphabet().removeChangeListener(ChangeListener.ALWAYS_VETO, ChangeType.UNKNOWN);
    }

    protected Map getNameToSymbol() {
        if (nameToSymbol == null) {
	    nameToSymbol = new HashMap();
	    for (Iterator i = ((FiniteAlphabet) getAlphabet()).iterator(); i.hasNext(); ) {
		Symbol sym = (Symbol) i.next();
		if (caseSensitive) {
		    nameToSymbol.put(sym.getName(), sym);
		} else {
		    nameToSymbol.put(sym.getName().toLowerCase(), sym);
		}
	    }
	    nameToSymbol.put("gap", getAlphabet().getGapSymbol());
	}

	return nameToSymbol;
    }

    public Symbol parseToken(String token)
        throws IllegalSymbolException
    {
	Symbol sym;
	if (caseSensitive) {
	    sym = (Symbol) getNameToSymbol().get(token);
	} else {
	    sym = (Symbol) getNameToSymbol().get(token.toLowerCase());
	}

	if (sym == null) {
	    char c = token.charAt(0);
	    if (c == '[') {
		if (token.charAt(token.length() - 1) != ']') {
		    throw new IllegalSymbolException("Mismatched parentheses: " + token);
		} else {
		    Symbol[] syms = parseString(token.substring(1, token.length() - 1));
		    Set ambigSet = new HashSet();
		    for (int i = 0; i < syms.length; ++i) {
			ambigSet.add(syms[i]);
		    }
		    return getAlphabet().getAmbiguity(ambigSet);
		}
	    } else {
		throw new IllegalSymbolException("Token `" + token + "' does not appear as a named symbol in alphabet `" + getAlphabet().getName() + "'");
	    }
	}
	return sym;
    }

    public String tokenizeSymbol(Symbol s) throws IllegalSymbolException {
	getAlphabet().validate(s);
	return s.getName();
    }
}
