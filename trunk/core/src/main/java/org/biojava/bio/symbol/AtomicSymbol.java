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
 * <p>
 * A symbol that is not ambiguous.
 * </p>
 *
 * <p>
 * Atomic symbols are the real underlying elements that a SymbolList is meant
 * to be composed of. DNA nucleotides are atomic, as are amino-acids. The
 * getMatches() method should return an alphabet containing just the one Symbol.
 * </p>
 *
 * <p>
 * The Symbol instances for single codons would be instances of AtomicSymbol as
 * they can only be represented as a Set of symbols from their alphabet that
 * contains just that one symbol.
 * </p>
 *
 * <p>
 * AtomicSymbol instances guarantee that getMatches returns an Alphabet
 * containing just that Symbol and each element of the List returned by
 * getSymbols is also atomic.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.1
 */
public interface AtomicSymbol extends BasisSymbol {
}
