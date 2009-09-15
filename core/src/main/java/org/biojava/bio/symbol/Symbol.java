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

import org.biojava.bio.Annotatable;

/**
 * A single symbol.
 * <p>
 * This is the atomic unit of a SymbolList, or a sequence. It allows
 * for fine-grain fly-weighting, so that there can be one instance
 * of each symbol that is referenced multiple times.
 * <p>
 * Symbols from finite alphabets are identifiable using the == operator.
 * Symbols from infinite alphabets may have some specific API to test for
 * equality, but should realy over-ride the equals() method.
 * <p>
 * Some symbols represent a single token in the sequence. For example, there is
 * a Symbol instance for adenine in DNA, and another one for cytosine.
 * Symbols can potentialy represent sets of Symbols. For example, n represents
 * any DNA Symbol, and X any protein Symbol. Gap represents the knowledge that
 * there is no Symbol. In addition, some symbols represent ordered lists of
 * other Symbols. For example, the codon agt can be represented by a single
 * Symbol from the Alphabet DNAxDNAxDNA. Symbols can represent ambiguity over
 * these complex symbols. For example, you could construct a Symbol instance
 * that represents the codons atn. This matches the codons {ata, att, atg, atc}.
 * It is also possible to build a Symbol instance that represents all stop
 * codons {taa, tag, tga}, which can not be represented in terms of a
 * single ambiguous n'tuple.
 * <p>
 * There are three Symbol interfaces. Symbol is the most generic. It has the
 * methods getToken and getName so that the Symbol can be textually represented.
 * In addition, it defines getMatches that returns an Alphabet over all the
 * AtomicSymbol instances that match the Symbol (N would return an Alphabet
 * containing {A, G, C, T}, and Gap would return {}).
 * <p>
 * BasisSymbol instances can always be represented by an n'tuple of BasisSymbol
 * instances. It adds the method getSymbols so that you can retrieve this list.
 * For example, the tuple [ant] is a BasisSymbol, as it is uniquely specified
 * with those three BasisSymbol instances a, n and t. n is a BasisSymbol
 * instance as it is uniquely represented by itself.
 * <p>
 * AtomicSymbol instances specialize BasisSymbol by guaranteeing that getMatches
 * returns a set containing only that instance. That is, they are indivisable.
 * The DNA nucleotides are instances of AtomicSymbol, as are individual codons.
 * The stop codon {tag} will have a getMatches method that returns {tag},
 * a getBases method that also returns {tag} and a getSymbols method that returns
 * the List [t, a, g]. {tna} is a BasisSymbol but not an AtomicSymbol as it
 * matches four AtomicSymbol instances {taa, tga, tca, tta}. It follows that
 * each symbol in getSymbols for an AtomicSymbol instance will also be
 * AtomicSymbol instances.
 *
 * @author Matthew Pocock
 */
public interface Symbol extends Annotatable {
  /**
   * The long name for the symbol.
   *
   * @return  the long name
   */
  String getName();
  
  /**
   * The alphabet containing the symbols matched by this ambiguity symbol.
   * <p>
   * This alphabet contains all of, and only, the symbols matched by this
   * symbol. For example, the symbol representing the DNA
   * ambiguity code for W would contain the symbol for A and T from the DNA
   * alphabet.
   *
   * @return  the Alphabet of symbols matched by this
   *          symbol
   */
  Alphabet getMatches();
}
