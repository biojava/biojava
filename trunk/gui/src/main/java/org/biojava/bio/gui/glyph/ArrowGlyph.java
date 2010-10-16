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
import java.awt.Shape;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;

import org.biojava.bio.seq.StrandedFeature;

/**
 * A Glyph that paints an arrow shape within the bounds. The
 * <code>setDirection</code> method allows the decision of the direction, to
 * which the arrow points.
 *
 * @author Mark Southern
 * @author <a href="mailto:andreas.draeger@uni-tuebingen.de">Andreas Dr&auml;ger</a>
 * @since 1.5
 */
public class ArrowGlyph implements Glyph {
	private Paint	            fillPaint;

	private Paint	            outerPaint;

	private Rectangle2D.Float	bounds;

	private Shape	            arrowShape;

	/**
	 * Creates a new <code>ArrowGlyph</code>, which is filled with the color
	 * blue by default.
	 */
	public ArrowGlyph() {
		this(Color.BLUE, Color.BLACK);
	}

	/**
	 * Creates a new <code>ArrowGlyph</code>, which is filled with the given
	 * color.
	 *
	 * @param fillPaint Paint properties to fill this arrow.
	 * @param outerPaint Paint properties of the outer border of this arrow.
	 */
	public ArrowGlyph(Paint fillPaint, Paint outerPaint) {
		this.fillPaint = fillPaint;
		this.outerPaint = outerPaint;
		this.bounds = new Rectangle2D.Float(0, 0, 0, 0);
	}

	/**
	 * This constructs an arrow in the given bounds, which is colored blue.
	 *
	 * @param bounds
	 */
	public ArrowGlyph(Rectangle2D.Float bounds) {
		this(bounds, Color.BLUE, Color.BLACK);
	}

	/**
	 * Constructor which sets both the size of this arrow and its color.
	 *
	 * @param bounds
	 * @param fillPaint
	 */
	public ArrowGlyph(Rectangle2D.Float bounds, Paint fillPaint, Paint outerPaint) {
		this(fillPaint, outerPaint);
		setBounds(bounds);
	}

	/*
	 * (non-Javadoc)
	 *
	 * see org.biojava.bio.gui.glyph.Glyph#getBounds()
	 */
	public Rectangle2D.Float getBounds() {
		return bounds;
	}

	/*
	 * (non-Javadoc)
	 *
	 * see org.biojava.bio.gui.glyph.Glyph#setBounds(java.awt.geom.Rectangle2D.Float)
	 */
	public void setBounds(Rectangle2D.Float r) {
		if (bounds.equals(r)) return;
		bounds = r;
	}

	/**
	 * This method allows you to decide in which direction the arrow has to point.
	 * The definition of directions is equal to the definition of
	 *  {@link StrandedFeature}.
	 *
	 * @param direction
	 *          A +1 means to the right, -1 to the left an 0 (and any other value)
	 *          produces a rectangular shape without arrows at its end.
	 */
	public void setDirection(int direction) {
		float q1 = bounds.height * 0.25F;
		float q2 = bounds.height * 0.5F;
		float q3 = bounds.height * 0.75F;
		float arrowHeadSize = bounds.height;
		GeneralPath p = new GeneralPath();

		switch (direction) {
		case +1: // to the right
			if ((bounds.width - arrowHeadSize) > 0) {
				p.moveTo(bounds.x, bounds.y + q1);
				p.lineTo(bounds.x + bounds.width - arrowHeadSize, bounds.y + q1);
				p.lineTo(bounds.x + bounds.width - arrowHeadSize, bounds.y);
				p.lineTo(bounds.x + bounds.width, bounds.y + q2);
				p.lineTo(bounds.x + bounds.width - arrowHeadSize, bounds.y
				    + bounds.height);
				p.lineTo(bounds.x + bounds.width - arrowHeadSize, bounds.y + q3);
				p.lineTo(bounds.x, bounds.y + q3);
			} else {
				p.moveTo(bounds.x, bounds.y);
				p.lineTo(bounds.x + bounds.width, bounds.y + q2);
				p.lineTo(bounds.x, bounds.y + bounds.height);
			}
			break;
		case -1: // to the left
			if ((bounds.width - arrowHeadSize) > 0) {
				p.moveTo(bounds.x + bounds.width, bounds.y + q1);
				p.lineTo(bounds.x + arrowHeadSize, bounds.y + q1);
				p.lineTo(bounds.x + arrowHeadSize, bounds.y);
				p.lineTo(bounds.x, bounds.y + q2);
				p.lineTo(bounds.x + arrowHeadSize, bounds.y + bounds.height);
				p.lineTo(bounds.x + arrowHeadSize, bounds.y + q3);
				p.lineTo(bounds.x + bounds.width, bounds.y + q3);
			} else {
				p.moveTo(bounds.x + bounds.width, bounds.y);
				p.lineTo(bounds.x + bounds.width, bounds.y + bounds.height);
				p.lineTo(bounds.x, bounds.y + q2);
			}
			break;
		default: // unknown
			// we cannot draw an arrow, we should draw a rectangle shape instead.
			p.moveTo(bounds.x, bounds.y + q1);
			p.lineTo(bounds.x + bounds.width, bounds.y + q1);
			p.lineTo(bounds.x + bounds.width, bounds.y + q3);
			p.lineTo(bounds.x, bounds.y + q3);
			break;
		}
		p.closePath();
		arrowShape = p;
	}

	public void render(Graphics2D g) {
		if ((bounds.height > 0) && (bounds.width > 0) && (arrowShape == null))
		  setDirection(0);
		if (arrowShape != null) {
			g.setPaint(fillPaint);
			g.fill(arrowShape);
			g.setPaint(outerPaint);
			g.draw(arrowShape);
		}
	}

	/**
	 * Returns the paint properties of this glyph.
	 *
	 * @return the fillPaint
	 */
	public Paint getFillPaint() {
		return fillPaint;
	}

	/**
	 * Allows you to set the paint properties of this glyph.
	 *
	 * @param fillPaint
	 */
	public void setFillPaint(Paint fillPaint) {
		this.fillPaint = fillPaint;
	}

	/**
	 * Returns the paint properties of the outer line of this glyph.
	 *
	 * @return the outerPaint
	 */
	public Paint getOuterPaint() {
		return outerPaint;
	}

	/**
	 * Allows setting the paint properties of the outer line of this glyph to the
	 * given value.
	 *
	 * @param outerPaint
	 */
	public void setOuterPaint(Paint outerPaint) {
		this.outerPaint = outerPaint;
	}

}
