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

package org.biojava.bio.symbol;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Stack;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.SymbolTokenization;

/**
 * <code>MotifTools</code> contains utility methods for sequence
 * motifs.
 *
 * @author Keith James
 */
public class MotifTools
{
    private static Symbol [] symProto = new Symbol [0];

    /**
     * <p><code>createRegex</code> creates a regular expression which
     * matches the <code>SymbolList</code>. Ambiguous
     * <code>Symbol</code>s are simply transformed into character
     * classes. For example the nucleotide sequence "AAGCTT" becomes
     * "A{2}GCT{2}" and "CTNNG" is expanded to
     * "CT[ABCDGHKMNRSTVWY]{2}G". The character class is generated
     * using the <code>getMatches</code> method of an ambiguity symbol
     * to obtain the alphabet of <code>AtomicSymbol</code>s it
     * matches, followed by calling <code>getAllSymbols</code> on this
     * alphabet, removal of any gap symbols and then tokenization of
     * the remainder. The ordering of the tokens in a character class
     * is by ascending numerical order of their tokens as determined
     * by <code>Arrays.sort(char [])</code>.</p>
     *
     * <p>The <code>Alphabet</code> of the <code>SymbolList</code>
     * must be finite and must have a character token type. Regular
     * expressions may be generated for any such
     * <code>SymbolList</code>, not just DNA, RNA and protein.</p>
     *
     * @param motif a <code>SymbolList</code>.
     *
     * @return a <code>String</code> regular expression.
     */
    public static String createRegex(SymbolList motif)
    {
        if (motif.length() == 0)
            throw new IllegalArgumentException("SymbolList was empty");

        StringBuffer regex = new StringBuffer();
        Stack stack = new Stack();

        try
        {
            SymbolTokenization sToke =
                motif.getAlphabet().getTokenization("token");
            if (sToke.getTokenType() != SymbolTokenization.CHARACTER)
                throw new IllegalArgumentException("SymbolList alphabet did not have a character token type");

            int motifLen = motif.length();

            for (int i = 1; i <= motifLen; i++)
            {
                StringBuffer sb = new StringBuffer();
                Symbol sym = motif.symbolAt(i);
                FiniteAlphabet ambiAlpha = (FiniteAlphabet) sym.getMatches();

                if (ambiAlpha.size() == 1)
                {
                    sb.append(sToke.tokenizeSymbol(sym));
                }
                else
                {
                    Set rawSyms = AlphabetManager.getAllSymbols(ambiAlpha);
                    List gapSyms = new ArrayList();

                    for (Iterator si = rawSyms.iterator(); si.hasNext();)
                    {
                        Symbol rawSym = (Symbol) si.next();
                        // Crude check for gap symbol
                        if (((FiniteAlphabet) rawSym.getMatches()).size() == 0)
                        {
                            gapSyms.add(rawSym);
                        }
                    }

                    rawSyms.removeAll(gapSyms);

                    // getAllSymbols returns a Set (i.e. unordered) so
                    // we convert to char array so we can sort tokens
                    Symbol [] ambiSyms = (Symbol []) rawSyms.toArray(symProto);
                    char [] ambiChars = new char [ambiSyms.length];

                    for (int j = 0; j < ambiSyms.length; j++)
                    {
                        ambiChars[j] =
                            sToke.tokenizeSymbol(ambiSyms[j]).charAt(0);
                    }

                    Arrays.sort(ambiChars);
                    sb.append(ambiChars);
                }

                String result = sb.substring(0);

                if (i == 1)
                {
                    stack.push(result);
                }
                else if (i < motifLen)
                {
                    if (! stack.isEmpty() && stack.peek().equals(result))
                    {
                        stack.push(result);
                    }
                    else
                    {
                        regex = extendRegex(stack, regex);
                        stack.push(result);
                    }
                }
                else
                {
                    if (! stack.isEmpty() && stack.peek().equals(result))
                    {
                        stack.push(result);
                        regex = extendRegex(stack, regex);
                    }
                    else
                    {
                        regex = extendRegex(stack, regex);
                        stack.push(result);
                        regex = extendRegex(stack, regex);
                    }
                }
            }
        }
        catch (IllegalSymbolException ise)
        {
            throw new BioError("Internal error: failed to tokenize a Symbol from an existing SymbolList", ise);
        }
        catch (BioException be)
        {
            throw new BioError("Internal error: failed to get SymbolTokenization for SymbolList alphabet", be);
        }

        return regex.substring(0);
    }

    private static StringBuffer extendRegex(Stack stack, StringBuffer regex)
    {
        String component = (String) stack.peek();

        if (component.length() == 1)
        {
            regex.append(component);

            if (stack.size() > 1)
            {
                regex.append("{");
                regex.append(stack.size());
                regex.append("}");
            }
        }
        else
        {
            regex.append("[");
            regex.append(component);
            regex.append("]");

            if (stack.size() > 1)
            {
                regex.append("{");
                regex.append(stack.size());
                regex.append("}");
            }
        }

        stack.clear();

        return regex;
    }
}
