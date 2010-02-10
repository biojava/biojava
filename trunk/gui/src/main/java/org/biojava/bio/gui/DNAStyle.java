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
import java.util.Map;

import org.biojava.bio.BioError;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * A simple implementation of SymbolStyle optimized for DNA.
 *
 * @author Matthew Pocock
 */
public class DNAStyle implements SymbolStyle {
  private Map outlinePaint;
  private Map fillPaint;

  {
    outlinePaint = new HashMap();
    fillPaint = new HashMap();
  }

  public Paint outlinePaint(Symbol s) throws IllegalSymbolException {
    DNATools.getDNA().validate(s);
    return (Paint) outlinePaint.get(s);
  }

  public Paint fillPaint(Symbol s) throws IllegalSymbolException {
    DNATools.getDNA().validate(s);
    return (Paint) fillPaint.get(s);
  }

  public void setOutlinePaint(Symbol s, Paint paint)
  throws IllegalSymbolException {
    DNATools.getDNA().validate(s);
    outlinePaint.put(s, paint);
  }

  public void setFillPaint(Symbol s, Paint paint)
  throws IllegalSymbolException {
    DNATools.getDNA().validate(s);
    fillPaint.put(s, paint);
  }

  public DNAStyle() {
    try {
      setOutlinePaint(DNATools.t(), Color.black);
      setFillPaint(DNATools.t(), Color.red);
      setOutlinePaint(DNATools.a(), Color.black);
      setFillPaint(DNATools.a(), Color.green);
      setOutlinePaint(DNATools.g(), Color.black);
      setFillPaint(DNATools.g(), Color.blue);
      setOutlinePaint(DNATools.c(), Color.black);
      setFillPaint(DNATools.c(), Color.yellow);
    } catch (IllegalSymbolException ire) {
      throw new BioError("DNA symbols dissapeared from DNA alphabet", ire);
    }
  }
}
