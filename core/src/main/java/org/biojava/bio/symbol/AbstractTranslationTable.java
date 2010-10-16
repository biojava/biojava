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

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * an abstract class implementing basic functionality
 * of a translation table that translates Symbols from
 * one Alphabet to another.
 *
 * @author Matthew Pocock
 * @author Keith James (docs)
 * @author Thomas Down
 * @author Greg Cox
 * @author Mark Schreiber
 * @author David Huen (refactoring)
 */
abstract class AbstractTranslationTable implements TranslationTable
{
    public abstract Alphabet getSourceAlphabet();
    public abstract Alphabet getTargetAlphabet();
    /**
     * this method is expected to translate any symbol
     * in the source alphabet.  Failure can be indicated
     * by returning a null if, for example, your method
     * only handles AtomicSymbols and you want BasisSymbols
     * to be taken apart.  If you are sure the symbol is
     * illegal, you can throw the IllegalSymbolException
     * immediately to bypass further processing.
     * <p>
     * As an optimisation, if your method is capable of
     * immediately translating an ambiguity Symbol, just
     * return it and the alternate route of establishing
     * the translation through doing an ambiguity
     * lookup will be avoided.
     */
    protected abstract Symbol doTranslate(Symbol sym) throws IllegalSymbolException;

    public Symbol translate(final Symbol sym)
    throws IllegalSymbolException {
        // make an attempt to translate immediately
        Symbol s = doTranslate(sym);

        // translation failed, validate and try an ambiguity lookup
        if(s == null) {
            if(sym instanceof AtomicSymbol) { //changed this from s to sym, since we already checked and s is null
                getSourceAlphabet().validate(sym);

                // the symbol was valid and still we can't handle it, bail out!
                throw new IllegalSymbolException("Unable to map " + sym.getName());
            } else {
                if(sym == null) {
                    throw new NullPointerException("Can't translate null");
                }
                Set syms = new HashSet();
                for (Iterator i = ((FiniteAlphabet) sym.getMatches()).iterator(); i.hasNext(); ) {
                    Symbol is = (Symbol) i.next();
                    syms.add(this.translate(is));
                }
                s = getTargetAlphabet().getAmbiguity(syms);
            }
        }
        return s;
    }
}
