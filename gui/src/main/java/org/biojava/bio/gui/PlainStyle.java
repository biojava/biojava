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

import org.biojava.bio.symbol.Symbol;

/**
 * A simple implementation of SymbolStyle that just uses a single paint for
 * outlines and a single paint for filling.
 *
 * @author Matthew Pocock
 */
public class PlainStyle implements SymbolStyle {
  private Paint outlinePaint;
  private Paint fillPaint;
  
  public Paint outlinePaint(Symbol s) {
    return outlinePaint;
  }
  
  public Paint fillPaint(Symbol s) {
    return fillPaint;
  }
  
  public PlainStyle() {
    this(Color.black, Color.gray);
  }
  
  public PlainStyle(Paint outlinePaint, Paint fillPaint) {
    this.outlinePaint = outlinePaint;
    this.fillPaint = fillPaint;
  }
}
