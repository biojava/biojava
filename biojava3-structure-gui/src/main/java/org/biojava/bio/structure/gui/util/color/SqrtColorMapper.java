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

/**
 * Performs a sqrt transform on input before passing the values off to another
 * colorMapper.
 * 
 * For instance, to map [0^2, 10^2] to a rainbow gradient, use
 * 
 * new LogColorMapper(GradientMapper.getGradientMapper(GradientMapper.RAINBOW_GRADIENT, 0, 10))
 * 
 * @author Spencer Bliven
 *
 */
public class SqrtColorMapper extends ContinuousColorMapperTransform {
	
	/**
	 * Creates a new SqrtColorMapper.
	 * @param sqrtspaceMapper 
	 */
	public SqrtColorMapper(ContinuousColorMapper sqrtspaceMapper) {
		super(sqrtspaceMapper);
	}
	
	/**
	 * Return sqrt(value).
	 * If value is negative, return the color corresponding to negative infinity.
	 * 
	 * @param value Value to be mapped 
	 * @return sqrt(value), or NEGATIVE_INFINITY
	 * @see org.biojava.bio.structure.gui.util.color.ContinuousColorMapper#getColor(double)
	 */
	public double transform(double value) {
		double sqrtValue = Double.NEGATIVE_INFINITY;
		if(value >= 0)
			sqrtValue = Math.sqrt(value);
		
		return sqrtValue;
	}

}
