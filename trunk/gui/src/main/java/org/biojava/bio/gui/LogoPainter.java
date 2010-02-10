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


/**
 * The interface for something that will draw the sequence logo for a state.
 * <p>
 * A StateLogo object claims the screen realestate for rendering, and does the
 * calculations for sizes & information and the like. The LogoPainter renders
 * this information onto a graphics context. It is given the StateLog to render,
 * so that a single LogoPainter can be shared among many state logos.
 *
 * @author Matthew Pocock
 */
public interface LogoPainter {
  /**
   * Render the state from sl onto the graphics object g.
   *
   * @param lCtxt the LogoContext to render
   */
  public void paintLogo(LogoContext lCtxt);
}
