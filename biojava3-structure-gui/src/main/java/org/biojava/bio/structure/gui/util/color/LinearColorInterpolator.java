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
import java.awt.color.ColorSpace;

/**
 * @author Spencer Bliven
 *
 */
public class LinearColorInterpolator implements ColorInterpolator {

	public enum InterpolationDirection {
		INNER, /* within [a,b] */
		OUTER, /* outside ]a,b[ */
		UPPER, /* above a, eg in [a, a+(b-a)%1] */
		LOWER  /* below a */
	}

	private ColorSpace colorSpace;
	private InterpolationDirection[] interpolationDirection;
	
	public LinearColorInterpolator() {
		this(ColorSpace.getInstance(ColorSpace.CS_sRGB));
	}
	public LinearColorInterpolator(ColorSpace colorSpace) {
		super();
		this.setColorSpace(colorSpace);
	}
	
	/**
	 * Interpolates to a color between a and b
	 * @param a First color
	 * @param b Second color
	 * @param mixing Mixing coefficient; the fraction of a in the result.
	 * @return The color between a and b
	 * @throws IllegalArgumentException if mixing is not between 0 and 1
	 * @see org.biojava.bio.structure.gui.util.color.ColorInterpolator#interpolate(java.awt.Color, java.awt.Color, float)
	 */
	public Color interpolate(Color a, Color b, float mixing) {
		float[] compA, compB;
		// Get components
		// Don't convert colorSpaces unless necessary
		if(a.getColorSpace().equals(colorSpace) ) {
			compA = a.getComponents(null);
		} else {
			compA = a.getComponents(colorSpace, null);
		}
		if(b.getColorSpace().equals(colorSpace)) {
			compB = b.getComponents(null);
		} else {
			compB = b.getComponents(colorSpace, null);
		}

		float[] compMixed = new float[compA.length];
		
		for(int i=0;i<compA.length;i++){
			//Normalizing to [0,1] after the interpolation,
			// INNER means between a and b
			// OUTER means between max(a,b) and min(a,b)+1
			// UPPER means between a and b' s.t. b'>a and b' in {b, b+1}
			// LOWER means between a and b' s.t. b'<a and b' in {b, b-1}
			float left, right;
			left = compA[i];
			//Alpha uses INNER direction
			InterpolationDirection dir = i<interpolationDirection.length ?
					interpolationDirection[i] : InterpolationDirection.INNER;
			switch(dir) {
			case INNER:
				right = compB[i];
				break;
			case OUTER:
				if(compA[i]<compB[i]) {
					right = compB[i]-1;
				} else {
					right = compB[i]+1;
				}
				break;
			case UPPER:
				if(compA[i]<compB[i]) {
					right = compB[i];
				} else {
					right = compB[i]+1;
				}
				break;
			case LOWER:
				if(compA[i]<compB[i]) {
					right = compB[i]-1;
				} else {
					right = compB[i];
				}
				break;
			default: throw new IllegalStateException("Unkown interpolation Direction "+interpolationDirection[i]);
			}
			
			//Perform mixing
			compMixed[i] = mixing*left + (1-mixing)*right;
			
			if(dir != InterpolationDirection.INNER) {
				//Normalize to [0,1]
				if(compMixed[i] < 0)
					compMixed[i] += 1f;
				if(compMixed[i] > 1)
					compMixed[i] -= 1f;
			}
		}
				
		return new Color(colorSpace,compMixed,compMixed[compMixed.length-1]);
	}


	/**
	 * Sets the ColorSpace to use for interpolation.
	 * 
	 * The most common scheme for color spaces is to use linear components
	 * between 0 and 1 (for instance red,green,blue). For such a component, a
	 * linear interpolation between two colors is used.
	 * 
	 * Sometimes a component may be in cylindrical coordinates. In this case,
	 * the component can be mapped in a number of ways. These are set by
	 * InterpolationDirections.
	 *
	 * @param colorSpace The color space for interpolation
	 * @param interpDirection An array of size colorSpace.getNumComponents()
	 * 		giving the interpolation direction for each component.
	 */
	public void setColorSpace(ColorSpace colorSpace, InterpolationDirection[] dir) {
		if(dir.length < colorSpace.getNumComponents()) {
			throw new IllegalArgumentException( "Must specify an interpolation " +
					"direction for each colorspace component ("+colorSpace.getNumComponents()+")");
		}
		this.colorSpace = colorSpace;
		this.interpolationDirection = dir;
	}

	public void setColorSpace(ColorSpace colorSpace) {
		InterpolationDirection[] dir = new InterpolationDirection[colorSpace.getNumComponents()];
		for(int i=0;i<dir.length;i++)
			dir[i] = InterpolationDirection.INNER;
		this.setColorSpace(colorSpace, dir);
	}
	
	public void setInterpolationDirection(int componentIndex, InterpolationDirection dir) {
		interpolationDirection[componentIndex] = dir;
	}
	
	
}
