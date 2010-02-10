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


/**
 * This class makes SimpleSymbolLists.
 * It could be refactored into the SimpleSymbolList class eventually.
 *
 * @author David Huen
 */
public class SimpleSymbolListFactory

    implements SymbolListFactory
{
    /**
     * Create a factory for SimpleSymbolLists.
     */

    public SymbolList makeSymbolList(Symbol [] symbolArray, int size, Alphabet alfa)
        throws IllegalAlphabetException
    {
        return new SimpleSymbolList(symbolArray, size, alfa);
    }
}

