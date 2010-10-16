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


/**
 * A Glyph that paints a two headed arrow within the bounds
 * 
 * @author Mark Southern
 * @since 1.5
 */
public class TwoHeadedArrowGlyph implements Glyph {
    private Paint fillPaint;
    private Rectangle2D.Float bounds = new Rectangle2D.Float(0, 0, 0, 0);
    private Shape arrowShape;

    public TwoHeadedArrowGlyph() {
        fillPaint = Color.BLUE;
    }

    public TwoHeadedArrowGlyph(Rectangle2D.Float bounds) {
        this();
        setBounds(bounds);
    }

    public Rectangle2D.Float getBounds() {
        return bounds;
    }

    public void setBounds(Rectangle2D.Float r) {
        if (bounds.equals(r)) {
            return;
        }

        float q1 = r.height * 0.25F;
        float q2 = r.height * 0.5F;
        float q3 = r.height * 0.75F;
        float arrowHeadSize = r.height;
        GeneralPath p = new GeneralPath();

        if ((r.width - (2F * arrowHeadSize)) > 0) {
            p.moveTo(r.x, r.y + q2);
            p.lineTo(r.x + arrowHeadSize, r.y);
            p.lineTo(r.x + arrowHeadSize, r.y + q1);
            p.lineTo((r.x + r.width) - arrowHeadSize, r.y + q1);
            p.lineTo((r.x + r.width) - arrowHeadSize, r.y);
            p.lineTo(r.x + r.width, r.y + q2);
            p.lineTo((r.x + r.width) - arrowHeadSize, r.y + r.height);
            p.lineTo((r.x + r.width) - arrowHeadSize, r.y + q3);
            p.lineTo(r.x + arrowHeadSize, r.y + q3);
            p.lineTo(r.x + arrowHeadSize, r.y + r.height);
        } else {
            p.moveTo(r.x, r.y + q1);
            p.lineTo(r.x + r.width, r.y + q1);
            p.lineTo(r.x + r.width, r.y + q3);
            p.lineTo(r.x, r.y + q3);
        }

        p.closePath();
        arrowShape = p;

        bounds = r;
    }

    public void render(Graphics2D g) {
        if (arrowShape != null) {
            g.setPaint(fillPaint);
            g.fill(arrowShape);
        }
    }
}
