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
package org.biojava.bio.program.abi;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.BioError;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.io.CharacterTokenization;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AbstractAlphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IntegerAlphabet;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolListViews;
import org.biojava.utils.ListTools;

/**
 * Useful functionality for working with fasta files where the quality of the
 * DNA is encoded as upper and lower case DNA characters.
 *
 * @author Matthew Pocock
 */
public class ABITools {
  /**
   * The quality alphabet. This is equivalent to DNA x [0,1] where 0 represents
   * poorly supported (lower case) and 1 represents strongly supported (upper
   * case).
   */
  public static final FiniteAlphabet QUALITY;

  /**
   * The poorly supported symbol.
   */
  public static final AtomicSymbol _0;

  /**
   * The well supported symbol.
   */
  public static final AtomicSymbol _1;

  /**
   * Alignment label for the DNA sequence row.
   */
  public static final Object SEQUENCE = "SEQUENCE";

  /**
   * Alignment label for the support row.
   */
  public static final Object SUPPORT = "SUPPORT";

  static {
    try {
      IntegerAlphabet.SubIntegerAlphabet _01
      = IntegerAlphabet.getSubAlphabet(0, 1);
      _0 = _01.getSymbol(0);
      _1 = _01.getSymbol(1);

      List alphas = new ArrayList();
      alphas.add(DNATools.getDNA());
      alphas.add(_01);

      // naughty here - we know because we are insiders that the result of this
      // call will be an AbstractAlphabet impl
      AbstractAlphabet quality = (AbstractAlphabet) AlphabetManager.getCrossProductAlphabet(alphas);
      CharacterTokenization tok = new CharacterTokenization(quality, true);

      // all lower case characters go to sym,0
      // all upper case characters go to sym,1
      SymbolList sl = DNATools.createDNA("agctrymkswhbvdn");
      ListTools.Doublet pair = new ListTools.Doublet();
      SymbolTokenization dnaTok = DNATools.getDNA().getTokenization("token");
      for(Iterator i = sl.iterator(); i.hasNext(); ) {
        pair.setA((Symbol) i.next());
        String c = dnaTok.tokenizeSymbol((Symbol) pair.getA());

        pair.setB(_1);
        tok.bindSymbol(quality.getSymbol(pair), c.toUpperCase().charAt(0));

        pair.setB(_0);
        tok.bindSymbol(quality.getSymbol(pair), c.toLowerCase().charAt(0));
      }

      quality.putTokenization("token", tok);
      QUALITY = quality;
    } catch (Exception e) {
      throw new BioError("Could not initialize ABI quality alphabet",e);
    }
  }

  /**
   * <p>
   * View a symbol list over the QUALITY alphabet as an alignment.
   * </p>
   *
   * <p>
   * The alignment will have labels of SEQUENCE and SUPPORT that retrieve the
   * DNA sequence and the binary support values respectively.
   * </p>
   *
   * @param abiSeq  the SymbolList over the QUALITY alphabet to view
   * @return an Alignment view of abiSeq
   * @throws IllegalAlphabetException if abiSeq is not over QUALITY
   */
  public static Alignment getAlignment(SymbolList abiSeq)
  throws IllegalAlphabetException {
    return SymbolListViews.alignment(
      new ListTools.Doublet(SEQUENCE, SUPPORT),
      abiSeq
    );
  }
}
