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
 */


package org.biojava.bio.structure;

import java.util.Collections;
import java.util.List;

import org.biojava.bio.BioError;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.BasisSymbol;
import org.biojava.bio.symbol.DoubleAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ListTools;


/**
 * AlphaCTools is a collection of static convenience methods for dealing with
 * Alpha Carbon Backbone Phi / Psi angles. In BioJava Phi Psi angles are represented
 * as Symbols from the ALPHA CARBON ANGLES alphabet. A Backbone could be represented
 * as a SymbolList. A Ramachandran plot might be represented as a Distribution over
 * Phi Psi Symbols.
 *
 * @author Mark Schreiber
 * @version 1.0
 */
public final class AlphaCTools {
    /** MAX_ANGLE . */
  public static final double MAX_ANGLE = 180.0;
    /** MIN_ANGLE . */
  public static final double MIN_ANGLE = -180.0;

  private static String ALPHA = "ALPHA CARBON ANGLES";
  private static DoubleAlphabet daInstance = DoubleAlphabet.getInstance();


  /**
   * Returns a reference to the Alphabet that contains Symbols that represent PHI,
   * PSI angles.
   *
   * @return a reference to the ALPHA CARBON ANGLES alphabet
   */
  public static Alphabet getAlphaCarbonAngleAlphabet(){
    if (AlphabetManager.registered(ALPHA)) {
      return AlphabetManager.alphabetForName(ALPHA);
    }
    else {
      List l = Collections.nCopies(2, DoubleAlphabet.getInstance());
      try {
        Alphabet a = AlphabetManager.getCrossProductAlphabet(l, ALPHA);
        AlphabetManager.registerAlphabet(ALPHA, a);

        return a;
      }
      catch (IllegalAlphabetException ex) {
        throw new BioError( "Cannot construct "+ALPHA+" alphabet",ex);
      }
    }
  }

  /**
   * Makes a Phi - Psi Symbol from the ALPHA CARBON ANGLES alphabet.
   *
   * @param phiAngle the phi angle between -180.0 and +180.0
   * @param psiAngle the psi angle between -180.0 and +180.0
   * @return a reference to the 'fly weight' Symbol.
   * @throws IllegalSymbolException if the bond angles are outside the specified range
   */
  public static Symbol getPhiPsiSymbol(double phiAngle, double psiAngle)
    throws IllegalSymbolException{

    if(phiAngle > MAX_ANGLE || phiAngle < MIN_ANGLE){
      throw new IllegalSymbolException("Phi angle must be between -180.0 and +180.0");
    }

    if(psiAngle > MAX_ANGLE || psiAngle < MIN_ANGLE){
      throw new IllegalSymbolException("Psi angle must be between -180.0 and +180.0");
    }

    Symbol phi = daInstance.getSymbol(phiAngle);
    Symbol psi = daInstance.getSymbol(psiAngle);

    return
      getAlphaCarbonAngleAlphabet().getSymbol(new ListTools.Doublet(phi, psi));
  }

  /**
   * extracts the Phi angle from a <code>Symbol</code>.
   *
   * @param phiPsiSym a <code>Symbol</code> from the ALPHA CARBON ANGLES
   * <code>Alphabet</code>
   * @return a double between -180.0 and +180.0
   * @throws IllegalSymbolException if the <code>Symbol</code> is not from
   *  the ALPHA CARBON ANGLES <code>Alphabet</code>
   */
  public static double getPhiAngle(Symbol phiPsiSym) throws IllegalSymbolException{
    //validate the Symbol
    getAlphaCarbonAngleAlphabet().validate(phiPsiSym);

    //get the phi angle
    List l = ((BasisSymbol)phiPsiSym).getSymbols();
    return ((DoubleAlphabet.DoubleSymbol)l.get(0)).doubleValue();
  }

  /**
   * extracts the Psi angle from a <code>Symbol</code>.
   * @param phiPsiSym a <code>Symbol</code> from the ALPHA CARBON ANGLES
   * <code>Alphabet</code>
   * @return a double between -180.0 and +180.0
   * @throws IllegalSymbolException if the <code>Symbol</code> is not from
   *  the ALPHA CARBON ANGLES <code>Alphabet</code>
   */
  public static double getPsiAngle(Symbol phiPsiSym) throws IllegalSymbolException{
    //validate the Symbol
    getAlphaCarbonAngleAlphabet().validate(phiPsiSym);

    //get the phi angle
    List l = ((BasisSymbol)phiPsiSym).getSymbols();
    return ((DoubleAlphabet.DoubleSymbol)l.get(1)).doubleValue();
  }
}
