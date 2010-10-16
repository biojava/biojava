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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;


/**
 * A Glyph that paints a wide 'H' line within the bounds
 *
 * @author Mark Southern
 * @since 1.5
 */
public class TurnGlyph implements Glyph {
    private Paint forePaint;
    private Stroke stroke;
    private Rectangle2D.Float bounds = new Rectangle2D.Float(0, 0, 0, 0);
    private Shape turnShape;

    public TurnGlyph() {
        forePaint = Color.YELLOW.darker();
        stroke = new BasicStroke(4f);
    }

    public TurnGlyph(Rectangle2D.Float bounds) {
        this();
        setBounds(bounds);
    }

    public TurnGlyph(Paint paint, Stroke stroke) {
        this.forePaint = paint;
        this.stroke = stroke;
    }

    public Rectangle2D.Float getBounds() {
        return bounds;
    }

    public void setBounds(Rectangle2D.Float r) {
        if (bounds.equals(r)) {
            return;
        }

        bounds = r;

        float q1 = r.height * 0.25F;
        float q2 = r.height * 0.5F;
        float q3 = r.height * 0.75F;

        GeneralPath p = new GeneralPath();
        p.moveTo(r.x, r.y + q3);
        p.lineTo(r.x, r.y + q1);
        p.lineTo(r.x, r.y + q2);
        p.lineTo(r.x + r.width, r.y + q2);
        p.lineTo(r.x + r.width, r.y + q1);
        p.lineTo(r.x + r.width, r.y + q3);
        turnShape = p;
    }

    public void render(Graphics2D g) {
        if (turnShape != null) {
            g.setStroke(stroke);
            g.setPaint(forePaint);
            g.draw(turnShape);
        }
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
	 * @param forePaint
	 */
	public void setPaint(Paint forePaint) {
		this.forePaint = forePaint;
	}
}
