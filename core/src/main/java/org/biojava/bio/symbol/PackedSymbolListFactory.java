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
 * This class makes PackedSymbolLists.
 * It could be refactored into the PackedSymbolList class eventually.
 *
 * @author David Huen
 */
public class PackedSymbolListFactory

    implements SymbolListFactory
{
    private final static int AUTO_SELECT = 1;
    private final static int USER_SELECT = 2;

    private boolean ambiguity;
    private int opMode = AUTO_SELECT;

    /**
     * Create a factory for PackedSymbolLists.
     * The use of ambiguity packing is determined automatically
     * as required.
     */
    public PackedSymbolListFactory()
    {
        opMode = AUTO_SELECT;
        ambiguity = false; // set to avoid javac bug
    }

    /**
     * Create a factory for PackedSymbolLists with specified packing type.
     *
     * @param ambiguity is ambiguity to be supported by the encoding?
     * @deprecated the argumentless constructor creates a SymbolListFactory
     *   that will autoselect the packing appropriately.
     */
    public PackedSymbolListFactory(boolean ambiguity)
    {
        opMode = USER_SELECT;
        this.ambiguity = ambiguity;
    }
 /**
   * Makes a packed SymbolList out of a list of Symbols.
   *
   * @param symbolArray the symbols to be make in a packed SymbolList
   * @param size the length of the symbolArray array.
   * @param alfa the Alphabet over which the SymbolList shoudl be
   * @return a packed SymbolList with the Symbols in symbolArray and the
   * Alphabet in alfa
   * @exception IllegalAlphabetException if alfa and the Symbols in
   *symbolArray disagree.
   */
    public SymbolList makeSymbolList(Symbol [] symbolArray, int size, Alphabet alfa)
        throws IllegalAlphabetException
    {
        switch (opMode) {
            case AUTO_SELECT:
                // check for ambiguity symbols
                ambiguity = false;
                for (int i=0; i < size; i++) {
                    if (!(symbolArray[i] instanceof AtomicSymbol)) { ambiguity = true; break; }
                }
            case USER_SELECT:
            default:
                return new PackedSymbolList(PackingFactory.getPacking((FiniteAlphabet) alfa, ambiguity), symbolArray, size, alfa);
        }
    }
}

