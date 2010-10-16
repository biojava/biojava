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
package org.biojava.bio.gui.sequence;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.event.MouseEvent;
import java.awt.geom.Area;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;

import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.symbol.Location;


/**
 * A Feature Renderer that paints the Feature as a right facing arrow Based heavily on
 * BasicFeatureRenderer
 *
 * @author Mark Southern
 * @since 1.5
 */
public class ArrowedFeatureRenderer implements FeatureRenderer {
    private Paint fill;
    private Paint outline;
    private double arrowSize = 15.0;
    private double arrowScoop = 4.0;
    private double arrowHeadWidth = 5;

    public ArrowedFeatureRenderer() {
        fill = Color.red;
        outline = Color.black;
    }

    public void setFill(Paint p) {
        fill = p;
    }

    public Paint getFill() {
        return fill;
    }

    public void setOutline(Paint p) {
        outline = p;
    }

    public Paint getOutline() {
        return outline;
    }

    public void setArrowSize(double arrowSize) {
        this.arrowSize = arrowSize;
    }

    public double getArrowSize() {
        return arrowSize;
    }

    public void setArrowScoop(double arrowScoop) {
        this.arrowScoop = arrowScoop;
    }

    public double getArrowScoop() {
        return arrowScoop;
    }

    public void setArrowHeadSize(double d) {
        this.arrowHeadWidth = d;
    }

    public double getArrowHeadSize() {
        return arrowHeadWidth;
    }

    public void renderFeature(Graphics2D g, Feature f, SequenceRenderContext src) {
        Shape s = null;
        Location loc = f.getLocation();
        float depth = ( float ) (arrowSize + (2.0 * arrowScoop));

        double minD;
        double maxD;

        if (src.getScale() > 1.0) {
            minD = src.sequenceToGraphics(loc.getMin());
            maxD = src.sequenceToGraphics(loc.getMax() + 1) - 1.0;
        } else {
            minD = src.sequenceToGraphics(loc.getMin());
            maxD = src.sequenceToGraphics(loc.getMax());
        }

        float min = ( float ) minD;
        float max = ( float ) maxD;

        float minBounds = ( float ) src.sequenceToGraphics(src.getRange().getMin() - 1);
        float maxBounds = ( float ) src.sequenceToGraphics(src.getRange().getMax() + 1);
        Shape bounds;

        if (src.getDirection() == SequenceRenderContext.HORIZONTAL) {
            bounds = new Rectangle2D.Double(minBounds, 0, maxBounds - minBounds, depth);
        } else {
            bounds = new Rectangle2D.Double(0, minBounds, depth, maxBounds - minBounds);
        }

        if ((max - min) >= arrowSize) {
            if (src.getDirection() == SequenceRenderContext.HORIZONTAL) {
                float minY = 0.0f;
                float maxY = depth;
                float minYS = minY + ( float ) arrowScoop;
                float maxYS = maxY - ( float ) arrowScoop;
                float midY = (minY + maxY) * 0.5f;
                float minX = min;
                float maxX = max;

                float midX1 = maxX - ( float ) getArrowHeadSize();
                float midX2 = minX + ( float ) getArrowHeadSize();

                GeneralPath path = new GeneralPath();
                path.moveTo(minX, midY);
                path.lineTo(midX2, minY);
                path.lineTo(midX2, minYS);
                path.lineTo(midX1, minYS);
                path.lineTo(midX1, minY);
                path.lineTo(maxX, midY);
                path.lineTo(midX1, maxY);
                path.lineTo(midX1, maxYS);
                path.lineTo(midX2, maxYS);
                path.lineTo(midX2, maxY);
                path.closePath();
                s = path;
            }
        }

        if (! bounds.contains(s.getBounds())) {
            //	System.err.println("Edgeclipping");
            s = new Area(s);
            (( Area ) s).intersect(new Area(bounds));
        }

        if (fill != null) {
            g.setPaint(fill);
            g.fill(s);
        }

        if ((outline != null) && ((maxD - minD) > 4.0)) {
            g.setPaint(outline);
            g.draw(s);
        } else {
            //	System.err.println("Not drawing outline...");
        }
    }

    public double getDepth(SequenceRenderContext src) {
        return arrowSize + (2.0 * arrowScoop) + 1.0;
    }

    public FeatureHolder processMouseEvent(FeatureHolder hits, SequenceRenderContext src,
        MouseEvent me
    ) {
        return hits;
    }
}
