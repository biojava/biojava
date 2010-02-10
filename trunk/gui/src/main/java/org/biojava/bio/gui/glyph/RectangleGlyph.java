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
package org.biojava.bio.gui.glyph;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.geom.Rectangle2D;

/**
 * A Glyph that paints a rectangle shape within the bounds.
 *
 * @author Mark Southern
 * @author <a href="mailto:andreas.draeger@uni-tuebingen.de>Andreas Dr&auml;ger</a>
 * @since 1.5
 */
public class RectangleGlyph implements Glyph {
	private Paint forePaint;
	private Rectangle2D.Float bounds = new Rectangle2D.Float(0, 0, 0, 0);

	public RectangleGlyph() {
		forePaint = Color.RED.brighter();
	}

	public RectangleGlyph(Rectangle2D.Float bounds) {
		this();
		setBounds(bounds);
	}

	public RectangleGlyph(Paint paint) {
		this.forePaint = paint;
	}

	public Rectangle2D.Float getBounds() {
		return bounds;
	}

	public void setBounds(Rectangle2D.Float r) {
		if (bounds.equals(r)) {
			return;
		}

		bounds = r;
	}

	public void render(Graphics2D g) {
		g.setPaint(forePaint);
		g.fill(bounds);
	}

	/**
	 *
	 * @return The currently set paint properties of this glyph.
	 */
	public Paint getPaint() {
		return forePaint;
	}

	/**
	 * Allows you to set the paint properties of this glyph, i.e., its color.
	 *
	 * @param forePaint
	 */
	public void setPaint(Paint forePaint) {
		this.forePaint = forePaint;
	}
}
