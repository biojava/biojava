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
package org.biojava.nbio.structure.symmetry.jmolScript;

import javax.vecmath.Color4f;
import java.awt.*;

public class ColorConverter {

	public static Color4f convertColor4f(Color color) {
		return new Color4f(color);
	}

	public static Color4f[] convertColor4f(Color[] colors) {
		Color4f[] colors4 = new Color4f[colors.length];
		for (int i = 0; i < colors.length; i++) {
			colors4[i] = convertColor4f(colors[i]);
		}
		return colors4;
	}
}
