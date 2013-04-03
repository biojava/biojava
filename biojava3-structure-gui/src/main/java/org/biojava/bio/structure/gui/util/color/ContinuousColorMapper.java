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
 * Provides a mapping between real numbers and Colors.
 * For instance, this could generate a gradient.
 * 
 * @author Spencer Bliven
 *
 */
public interface ContinuousColorMapper {

	/**
	 * 
	 * @param value The real to be mapped
	 * @return The color corresponding to value
	 */
	public Color getColor(double value);
}
