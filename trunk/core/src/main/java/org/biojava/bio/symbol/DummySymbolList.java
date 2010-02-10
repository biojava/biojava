/*
 * BioJava development code
 * 
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 * 
 * http://www.gnu.org/copyleft/lesser.html
 * 
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 * 
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 * 
 * http://www.biojava.org
 *
 */

package org.biojava.bio.symbol;

import java.io.Serializable;

/**  
 * Symbol list which just consists of non-informative symbols.
 * A DummySymbolList can be constructed over any Alphabet, and may
 * be of any length.  Calls to the symbolAt method will always return
 * the non-informative symbol for the alphabet in question (i.e.
 * 'n' for DNA, 'X' for protein, etc.).
 *
 * If you wish to work with <code>Feature</code> objects, but don't
 * have the actual sequence data available, you can construct a
 * <code>SimpleSequence</code> from a <code>DummySequence</code>,
 * and create features. on that.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.2
 */

public class DummySymbolList extends AbstractSymbolList implements Serializable {
    private final Symbol sym;
    private final Alphabet alpha;
    private final int length;

    public DummySymbolList(FiniteAlphabet alpha, int length) {
        super();
        this.alpha = alpha;
        this.length = length;
        sym = AlphabetManager.getAllAmbiguitySymbol(alpha);
    }
    
    public DummySymbolList(Alphabet alpha, int length, Symbol sym)
    throws IllegalSymbolException {
        alpha.validate(sym);
        
        this.alpha = alpha;
        this.length = length;
        this.sym = sym;
    }

    public Alphabet getAlphabet() {
        return alpha;
    }

    public int length() {
        return length;
    }

    public Symbol symbolAt(int i) {
        return sym;
    }
}
    
