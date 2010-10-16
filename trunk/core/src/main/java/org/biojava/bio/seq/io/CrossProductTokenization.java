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

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.BasisSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * Tokenization for cross-product alphabets.  This class handles
 * the general case of tokens of the form (foo bar baz), where
 * each element is handled by a sub-tokenization.  By default,
 * these will be the "name" tokenizations of each of the sub-alphabets,
 * but any tokenization can be used.
 *
 * @author Thomas Down
 * @author Greg Cox
 * @since 1.2
 */

public class CrossProductTokenization extends WordTokenization {
    private List subTokenizations;  // List<SymbolTokenization>

    public CrossProductTokenization(Alphabet alpha)
        throws BioException
    {
	super(alpha);
	subTokenizations = new ArrayList();
	for (Iterator i = alpha.getAlphabets().iterator(); i.hasNext(); ) {
	    Alphabet subAlpha = (Alphabet) i.next();
	    subTokenizations.add(subAlpha.getTokenization("name"));
	}
    }

    public CrossProductTokenization(Alphabet alpha,
				    List tokenizers)
    {
	super(alpha);
	this.subTokenizations = tokenizers;
	// Ought to validate...
    }

    public Symbol parseToken(String token)
        throws IllegalSymbolException
    {
	char c = token.charAt(0);
	if (c == '(') {
	    if (token.charAt(token.length() - 1) != ')') {
		throw new IllegalSymbolException("Mismatched parentheses: " + token);
	    } else {
		List split = splitString(token.substring(1, token.length() - 1));
		List syms = new ArrayList();

		Iterator si = split.iterator();
		Iterator ti = subTokenizations.iterator();
		while (si.hasNext()) {
		    String subToken = (String) si.next();
		    SymbolTokenization subTokenization = (SymbolTokenization) ti.next();
		    syms.add(subTokenization.parseToken(subToken));
		}

		return getAlphabet().getSymbol(syms);
	    }
	} else if (c == '[') {
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
	    throw new IllegalSymbolException("Not in standard cross-product form: " + token);
	}
    }

    public String tokenizeSymbol(Symbol s) throws IllegalSymbolException {
	getAlphabet().validate(s);

	if (s instanceof BasisSymbol) {
	    StringBuffer sb = new StringBuffer();
	    sb.append('(');
	    Iterator si = ((BasisSymbol) s).getSymbols().iterator();
	    Iterator ti = subTokenizations.iterator();

	    while (si.hasNext()) {
		Symbol subSym = (Symbol) si.next();
		SymbolTokenization subToke = (SymbolTokenization) ti.next();

		sb.append(subToke.tokenizeSymbol(subSym));
		if (si.hasNext()) {
		    sb.append(' ');
		}
	    }
	    sb.append(')');
	    return sb.substring(0);
	} else {
	    StringBuffer sb = new StringBuffer();
	    sb.append('[');
	    Iterator si = ((FiniteAlphabet) s.getMatches()).iterator();
	    while (si.hasNext()) {
		Symbol aSym = (Symbol) si.next();
		sb.append(tokenizeSymbol(aSym));
		if (si.hasNext()) {
		    sb.append(' ');
		}
	    }
	    sb.append(']');
	    return sb.substring(0);
	}
    }
}
