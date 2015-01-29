/*
 *                  BioJava development code
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
 * Created on Aug 3, 2007
 */
package org.biojava.bio.structure.gui.util.color;

import java.awt.Color;

/**
 * @author Spencer Bliven
 *
 */
public interface ColorInterpolator {
	/**
	 * Interpolates to a color between a and b
	 * @param a First color
	 * @param b Second color
	 * @param mixing Mixing coefficient; the fraction of a in the result.
	 * @return The color between a and b
	 * @throws IllegalArgumentException if mixing is not between 0 and 1
	 */
	public Color interpolate(Color a, Color b, float mixing);
}
