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
 * Maps colors by performing a transform of the input data and then passing 
 * the transformed value to a ContinuousColorMapper for rendering.
 * 
 * For instance, to map [10^0, 10^10] to a rainbow gradient, use
 * 
 * new LogColorMapper(GradientMapper.getGradientMapper(GradientMapper.RAINBOW_GRADIENT, 0, 10))
 * 
 * @author Spencer Bliven
 *
 */
public abstract class ContinuousColorMapperTransform implements ContinuousColorMapper {
	
	protected ContinuousColorMapper mapper;

	/**
	 * Creates a transform.
	 * @param submapper A color mapper which acts on the transformed input value. 
	 */
	public ContinuousColorMapperTransform(ContinuousColorMapper submapper) {
		this.mapper = submapper;
	}

	/**
	 * Return the color corresponding to submapper.getColor(transform(value)).
	 * 
	 * @param value Value to be mapped
	 * @return color corresponding to transform(value)
	 * @see org.biojava.bio.structure.gui.util.color.ContinuousColorMapper#getColor(double)
	 */
	public Color getColor(double value) {
		return mapper.getColor(transform(value));
	}

	/**
	 * An arbitrary transform over reals
	 * @param the input value
	 * @return the transformed value
	 */
	public abstract double transform(double value);
}
