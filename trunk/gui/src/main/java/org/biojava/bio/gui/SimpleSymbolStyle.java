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


package org.biojava.bio.gui;

import java.awt.Color;
import java.awt.Paint;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.BioError;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * A no-frills implementation of SymbolStyle.
 *
 * @author Matthew Pocock
 */
public class SimpleSymbolStyle implements SymbolStyle {
  private final Map outlinePaint;
  private final Map fillPaint;
  private final FiniteAlphabet alphabet;

  {
    outlinePaint = new HashMap();
    fillPaint = new HashMap();
  }

  public Alphabet getAlphabet() {
    return alphabet;
  }

  public SimpleSymbolStyle(FiniteAlphabet alphabet) {
    this.alphabet = alphabet;
    Map outline = getStandardOutlinePaints(alphabet);
    Map fill = getStandardFillPaints(alphabet);
    try {
      if(fill == null || outline == null) {
        for(Iterator i = alphabet.iterator(); i.hasNext(); ) {
          Symbol s = (Symbol) i.next();
          if(outline == null) {
            setOutlinePaint(s, Color.black);
          } else {
            setOutlinePaint(s, (Paint) outline.get(s));
          }
          if(fill == null) {
            setFillPaint(s, Color.gray);
          } else {
            setOutlinePaint(s, (Paint) fill.get(s));
          }
        }
      }
    } catch (IllegalSymbolException ire) {
      throw new BioError("Symbol dissapeared from my alphabet", ire);
    }
  }

  public Paint outlinePaint(Symbol s) throws IllegalSymbolException {
    getAlphabet().validate(s);
    return (Paint) outlinePaint.get(s);
  }

  public Paint fillPaint(Symbol s) throws IllegalSymbolException {
    getAlphabet().validate(s);
    return (Paint) fillPaint.get(s);
  }

  public void setOutlinePaint(Symbol s, Paint paint)
  throws IllegalSymbolException {
    getAlphabet().validate(s);
    outlinePaint.put(s, paint);
  }

  public void setFillPaint(Symbol s, Paint paint)
  throws IllegalSymbolException {
    getAlphabet().validate(s);
    fillPaint.put(s, paint);
  }

  private static Map standardFillPaints;
  private static Map standardOutlinePaints;

  public static Map getStandardFillPaints(Alphabet alpha) {
    return (Map) standardFillPaints.get(alpha);
  }

  public static Map getStandardOutlinePaints(Alphabet alpha) {
    return (Map) standardOutlinePaints.get(alpha);
  }

  static {
    standardFillPaints = new HashMap();
    standardOutlinePaints = new HashMap();

    Map dnaFill = new HashMap();
    dnaFill.put(DNATools.t(), Color.red);
    dnaFill.put(DNATools.g(), Color.blue);
    dnaFill.put(DNATools.c(), Color.yellow);
    dnaFill.put(DNATools.a(), Color.green);
    standardFillPaints.put(DNATools.getDNA(), dnaFill);

    Map dnaOutline = new HashMap();
    dnaOutline.put(DNATools.t(), Color.black);
    dnaOutline.put(DNATools.a(), Color.black);
    dnaOutline.put(DNATools.g(), Color.black);
    dnaOutline.put(DNATools.c(), Color.black);
    standardOutlinePaints.put(DNATools.getDNA(), dnaOutline);
  }
}
