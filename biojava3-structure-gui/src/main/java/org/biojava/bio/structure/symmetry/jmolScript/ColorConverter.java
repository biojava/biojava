package org.biojava.bio.structure.symmetry.jmolScript;

import java.awt.Color;

import javax.vecmath.Color4f;

public class ColorConverter {

	public static Color4f convertColor4f(Color color) {
		return new Color4f(color);
	}
	
	public static Color4f[] convertColor4f(Color colors[]) {
		Color4f[] colors4 = new Color4f[colors.length];
		for (int i = 0; i < colors.length; i++) {
			colors4[i] = convertColor4f(colors[i]);
		}
		return colors4;
	}
}
