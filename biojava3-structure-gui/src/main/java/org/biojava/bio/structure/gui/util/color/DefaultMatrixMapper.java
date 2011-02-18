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
 * Color Mapper which mimics the default coloring of JMatrixPanel pixels.
 * 
 * Assumes inputs in the range [0,max]. These are mapped to HSB colors such
 * that the hue and brightness are from [1,0].
 * @author Spencer Bliven
 *
 */
public class DefaultMatrixMapper implements ContinuousColorMapper {

	private double scalevalue;
	private float saturation;
	
	public DefaultMatrixMapper(double scale, float saturation ) {
		this.scalevalue = scale;
		this.saturation = saturation;
	}
	/**
	 * @param value
	 * @return
	 * @see org.biojava.bio.structure.gui.util.color.ContinuousColorMapper#getColor(double)
	 */
	public Color getColor(double value) {
		float hue = 1.0f;
		hue = (float)(1-(value/scalevalue));
		if (hue < 0)
			hue = 0;
		
		return Color.getHSBColor(hue,saturation,hue);
	}

	/**
	 * @return the scalevalue
	 */
	public double getScalevalue() {
		return scalevalue;
	}
	/**
	 * @param scalevalue the scalevalue to set
	 */
	public void setScalevalue(double scalevalue) {
		this.scalevalue = scalevalue;
	}
	/**
	 * @return the saturation
	 */
	public float getSaturation() {
		return saturation;
	}
	/**
	 * @param saturation the saturation to set
	 */
	public void setSaturation(float saturation) {
		this.saturation = saturation;
	}
}
