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
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;


/**
 * A Glyph that paints a Helix within the bounds
 *
 * @author Mark Southern
 * @since 1.5
 */
public class HelixGlyph implements Glyph {
    private Paint forePaint;
    private Paint backPaint;
    private Rectangle2D.Float bounds = new Rectangle2D.Float(0, 0, 0, 0);
    private Shape foreShape;
    private Shape backShape;

    public HelixGlyph() {
        forePaint = Color.RED;
        backPaint = Color.RED.darker().darker();
    }

    public HelixGlyph(Rectangle2D.Float bounds) {
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

        int[] vals = null;
        float hh = 0;

        //if (r.height != bounds.height) {
        hh = r.height / 2F;
        vals = new int[ 2 * ( int ) hh ];
        //System.out.println("length\t"+vals.length);

        double dd = (2 * Math.PI) / vals.length;
        //System.out.println("dd\t"+dd);
        
        for (int i = 0; i < vals.length; i++) {
            double v = ( double ) hh * Math.cos(dd * ( double ) i);
            vals[ i ] = ( int ) v;
            //System.out.println(i + "\t" + vals[i] );
        }

        //}
        //if(r.width != bounds.width && vals != null){
        List p1 = new ArrayList(2 * vals.length);
        List p2 = new ArrayList(2 * vals.length);

        for (int i = 0; i < r.width; i++) {
            int j = i % vals.length;

            if (j < (vals.length / 2)) {
                p1.add(new Point2D.Float(r.x + i, (r.y + hh) - vals[ j ]));
            } else {
                p1.add(new Point2D.Float(r.x + i, r.y));
            }

            if (vals[ j ] < 0) {
                p2.add(new Point2D.Float(r.x + i, r.y + hh + vals[ j ]));
            }
        }

        for (int i = ( int ) r.width - 1; i >= 0; i--) {
            int j = i % vals.length;

            if (j < (vals.length / 2)) {
                p1.add(new Point2D.Float(r.x + i, r.y + r.height));
            } else {
                p1.add(new Point2D.Float(r.x + i, r.y + hh + vals[ j ]));
            }

            if (vals[ j ] < 0) {
                p2.add(new Point2D.Float(r.x + i, (r.y + hh) - vals[ j ]));
            }
        }

        foreShape = GlyphUtils.getShape(p1);
        backShape = GlyphUtils.getShape(p2);

        //}
        bounds = r;
    }

    public void render(Graphics2D g) {
        if (backShape != null) {
            g.setPaint(backPaint);
            g.fill(backShape);
        }

        if (foreShape != null) {
            g.setPaint(forePaint);
            g.fill(foreShape);
        }
    }
}
