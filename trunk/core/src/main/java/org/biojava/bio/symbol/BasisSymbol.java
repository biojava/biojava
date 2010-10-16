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

import java.util.List;

/**
 * <p>
 * A symbol that can be represented as a string of Symbols.
 * </p>
 *
 * <p>
 * BasisSymbol instances can always be represented uniquely as a single List of
 * BasisSymbol instances. The symbol N is a BasisSymbol - it can be uniquely
 * represented by N. It matches {a, g, c, t}.
 * Similarly, the symbol atn is a BasisSymbol, as it can be uniquely
 * represented with a single list of symbols [a, t, n]. Its getMatches will
 * return the set {ata, att, atg, atc}.
 * </p>
 *
 * <p>
 * The getSymbols method returns the unique list of BasisSymbol instances that
 * this is composed from. For example, the codon ambiguity symbol atn will have
 * a getSymbols method that returns the list [a, t, n]. The getMatches method
 * will return an Alphabet containing each AtomicSymbol that can be made by
 * expanding the list of BasisSymbol instances.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.1
 */
public interface BasisSymbol extends Symbol {
  /**
   * <p>
   * The list of symbols that this symbol is composed from.
   * </p>
   *
   * <p>
   * In the usual case, this list will contain just this single symbol. In the
   * case where a symbol represents an ordered combination of other symbols,
   * the list will contain each of these BasisSymbols.
   * </p>
   *
   * @return the List of Symbols that this Symbol is built from
   */
  List getSymbols();
}
