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
 * Performs a log10 transform on input before passing the values off to another
 * colorMapper.
 * 
 * For instance, to map [10^0, 10^10] to a rainbow gradient, use
 * 
 * new LogColorMapper(GradientMapper.getGradientMapper(GradientMapper.RAINBOW_GRADIENT, 0, 10))
 * 
 * @author Spencer Bliven
 *
 */
public class LogColorMapper extends ContinuousColorMapperTransform {
	
	private int base;

	/**
	 * Creates a new LogColorMapper with base 10.
	 * @param logspaceMapper 
	 */
	public LogColorMapper(ContinuousColorMapper logspaceMapper) {
		this(logspaceMapper, 10);
	}
	
	/**
	 * If logspaceMapper maps values x1 to x2, this creates a
	 * mapper for values base^x1 to base^x2
	 * 
	 * @param logspaceMapper logspace mapper 
	 * @param base The base of the logorithm
	 */
	public LogColorMapper(ContinuousColorMapper logspaceMapper, int base) {
		super(logspaceMapper);
		this.base = base;
	}

	/**
	 * Apply log transform.
	 * @param value
	 * @return log_b(value)
	 * @see org.biojava.bio.structure.gui.util.color.ContinuousColorMapperTransform#transform(double)
	 */
	@Override
	public double transform(double value) {
		double logValue = Math.log(value>0?value:0)/Math.log(base);
		return logValue;
	}

}
