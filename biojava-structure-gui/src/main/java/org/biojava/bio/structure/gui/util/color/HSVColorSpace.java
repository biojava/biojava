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
public class HSVColorSpace extends ColorSpace {
	
	private static final long serialVersionUID = 8324413992279510075L;

	public static void main(String args[]) {
		HSVColorSpace csHSV = new HSVColorSpace();
		ColorSpace csRGB = ColorSpace.getInstance(CS_sRGB);
		Color c;
		float[] rgbComp;
		float[] hsvComp;
	
		c = Color.RED;
		rgbComp = c.getColorComponents(csRGB,null);
		hsvComp = c.getColorComponents(csHSV, null);
		assert(rgbComp.length == 3);
		assert(hsvComp.length == 3);
		System.out.format("RED\tRGB[%f %f %f] HSV[%f %f %f]\n",
				rgbComp[0], rgbComp[1], rgbComp[2],
				hsvComp[0], hsvComp[1], hsvComp[2] );
		
		
		c = Color.WHITE;
		rgbComp = c.getColorComponents(csRGB,null);
		hsvComp = c.getColorComponents(csHSV, null);
		System.out.format("WHITE\tRGB[%f %f %f] HSV[%f %f %f]\n",
				rgbComp[0], rgbComp[1], rgbComp[2],
				hsvComp[0], hsvComp[1], hsvComp[2] );
		
		c = Color.BLACK;
		rgbComp = c.getColorComponents(csRGB,null);
		hsvComp = c.getColorComponents(csHSV, null);
		System.out.format("BLACK\tRGB[%f %f %f] HSV[%f %f %f]\n",
				rgbComp[0], rgbComp[1], rgbComp[2],
				hsvComp[0], hsvComp[1], hsvComp[2] );
		
		c = Color.GRAY;
		rgbComp = c.getColorComponents(csRGB,null);
		hsvComp = c.getColorComponents(csHSV, null);
		System.out.format("GRAY\tRGB[%f %f %f] HSV[%f %f %f]\n",
				rgbComp[0], rgbComp[1], rgbComp[2],
				hsvComp[0], hsvComp[1], hsvComp[2] );
		
		c = Color.CYAN;
		rgbComp = c.getColorComponents(csRGB,null);
		hsvComp = c.getColorComponents(csHSV, null);
		System.out.format("CYAN\tRGB[%f %f %f] HSV[%f %f %f]\n",
				rgbComp[0], rgbComp[1], rgbComp[2],
				hsvComp[0], hsvComp[1], hsvComp[2] );
		
		
	}
	/**
	 * @param type
	 * @param numcomponents
	 */
	public HSVColorSpace() {
		super(ColorSpace.TYPE_HSV, 3);
	}

	/**
	 * @param ciexyzvalue The color components in CIEXYZ colorspace
	 * @return the equivalent HSV components
	 * @see java.awt.color.ColorSpace#fromCIEXYZ(float[])
	 */
	@Override
	public float[] fromCIEXYZ(float[] ciexyzvalue) {
		ColorSpace CIEXYZcs = ColorSpace.getInstance(CS_CIEXYZ);
		float[] rgb = CIEXYZcs.toRGB(ciexyzvalue);
		return this.fromRGB(rgb);
	}

	/**
	 * @param rgbvalue
	 * @return
	 * @see java.awt.color.ColorSpace#fromRGB(float[])
	 */
	@Override
	public float[] fromRGB(float[] rgbvalue) {
		assert(rgbvalue.length==3);
		Color rgbColor = new Color(rgbvalue[0], rgbvalue[1], rgbvalue[2]);
		return Color.RGBtoHSB(rgbColor.getRed(),rgbColor.getGreen(),rgbColor.getBlue(), null);
	}

	/**
	 * @param hsv
	 * @return
	 * @see java.awt.color.ColorSpace#toCIEXYZ(float[])
	 */
	@Override
	public float[] toCIEXYZ(float[] hsv) {
		float[] rgb = this.toRGB(hsv);
		ColorSpace CIEXYZcs = ColorSpace.getInstance(CS_CIEXYZ);
		return CIEXYZcs.fromRGB(rgb);
	}

	/**
	 * @param hsv 3-component array specifying hue, saturation, value
	 * @return equivalent red,green,blue components
	 * @see java.awt.color.ColorSpace#toRGB(float[])
	 */
	@Override
	public float[] toRGB(float[] hsv) {
		int rgb = Color.HSBtoRGB(hsv[0], hsv[1], hsv[2]);
		Color rgbColor = new Color(rgb);
		return rgbColor.getColorComponents(null);
	}
	
	private static HSVColorSpace hsvSpace;
	/**
	 * The HSV color space
	 */
	public static final int CS_HSV = 1007;
	public static ColorSpace getInstance(int colorspace) {
		ColorSpace theColorSpace;
		switch( colorspace ) {
		case CS_HSV:
			synchronized(HSVColorSpace.class) { 
				if(hsvSpace == null) {
					hsvSpace = new HSVColorSpace();
				}
				theColorSpace = hsvSpace;
			}
			break;
		default:
			theColorSpace = ColorSpace.getInstance(colorspace);
		}
		
		return theColorSpace;
	}
	
	public static ColorSpace getHSVColorSpace() {
		return HSVColorSpace.getInstance(CS_HSV);
	}

}
