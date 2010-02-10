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
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.utils.Unchangeable;

/**
 * a simple indexer for a full DNA alphabet including
 * ambiguities.
 *
 * @author David Huen
 * @author Matthew Pocock
 */
class FullDnaAlphabetIndex
  extends
    Unchangeable
  implements
    AlphabetIndex, java.io.Serializable
{
    static Symbol [] symbolArray = null;
    FiniteAlphabet dna = null;
    AlphabetIndex alphaIndex = null;

    public FullDnaAlphabetIndex()
        throws BioException
    {
        if (symbolArray == null) makeSymbolArray();
    }

    public FiniteAlphabet getAlphabet()
    {
        return dna;
    }

    public int indexForSymbol(Symbol sym)
        throws IllegalSymbolException
    {
        // linear search to find entry in my table
        for (int i=0; i<16; i++) {
            if (sym == symbolArray[i]) return i;
        }

        throw new IllegalSymbolException("can't find symbol " + sym);
    }

    public Symbol symbolForIndex(int index)
    {
        if ((index < 0) || (index > 15)) throw new IndexOutOfBoundsException();

        return symbolArray[index];
    }

    /**
     * private utilities
     */
    private void makeSymbolArray()
        throws IllegalSymbolException
    {

        // create an array of symbol mappings
        symbolArray = new Symbol[16];

        // get DNA alphabet
        if (dna == null) dna = DNATools.getDNA();

        // setup gap symbol
        symbolArray[4] = AlphabetManager.getGapSymbol();

        // set up atomic symbols
        symbolArray[0] = DNATools.a();
        symbolArray[1] = DNATools.c();
        symbolArray[2] = DNATools.g();
        symbolArray[3] = DNATools.t();

        // setup ambiguity symbols
        Set ambiguitySet = new HashSet(3);

        ambiguitySet.add(DNATools.g());
        ambiguitySet.add(DNATools.t());
        symbolArray[5] = dna.getAmbiguity(ambiguitySet);

        ambiguitySet.clear();
        ambiguitySet.add(DNATools.c());
        ambiguitySet.add(DNATools.t());
        symbolArray[6]=dna.getAmbiguity(ambiguitySet);

        ambiguitySet.clear();
        ambiguitySet.add(DNATools.c());
        ambiguitySet.add(DNATools.g());
        symbolArray[7]=dna.getAmbiguity(ambiguitySet);

        ambiguitySet.clear();
        ambiguitySet.add(DNATools.c());
        ambiguitySet.add(DNATools.g());
        ambiguitySet.add(DNATools.t());
        symbolArray[8]=dna.getAmbiguity(ambiguitySet);

        ambiguitySet.clear();
        ambiguitySet.add(DNATools.a());
        ambiguitySet.add(DNATools.t());
        symbolArray[9]=dna.getAmbiguity(ambiguitySet);

        ambiguitySet.clear();
        ambiguitySet.add(DNATools.a());
        ambiguitySet.add(DNATools.g());
        symbolArray[10]=dna.getAmbiguity(ambiguitySet);

        ambiguitySet.clear();
        ambiguitySet.add(DNATools.a());
        ambiguitySet.add(DNATools.g());
        ambiguitySet.add(DNATools.t());
        symbolArray[11]=dna.getAmbiguity(ambiguitySet);

        ambiguitySet.clear();
        ambiguitySet.add(DNATools.a());
        ambiguitySet.add(DNATools.c());
        symbolArray[12]=dna.getAmbiguity(ambiguitySet);

        ambiguitySet.clear();
        ambiguitySet.add(DNATools.a());
        ambiguitySet.add(DNATools.c());
        ambiguitySet.add(DNATools.t());
        symbolArray[13]=dna.getAmbiguity(ambiguitySet);

        ambiguitySet.clear();
        ambiguitySet.add(DNATools.a());
        ambiguitySet.add(DNATools.c());
        ambiguitySet.add(DNATools.g());
        symbolArray[14]=dna.getAmbiguity(ambiguitySet);

        ambiguitySet.clear();
        ambiguitySet.add(DNATools.a());
        ambiguitySet.add(DNATools.c());
        ambiguitySet.add(DNATools.g());
        ambiguitySet.add(DNATools.t());
        symbolArray[15]=dna.getAmbiguity(ambiguitySet);
    }

}

