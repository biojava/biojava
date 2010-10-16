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

import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Stroke;
import java.awt.geom.RoundRectangle2D;

import org.biojava.bio.seq.Feature;
import org.biojava.bio.symbol.Location;

/**
 * <code>RoundRectangularBeadRenderer</code> renders features
 * as rectangles with rounded corners. Their outline and fill
 * <code>Paint</code>, <code>Stroke</code>, feature depth, Y-axis
 * displacement are configurable.
 *
 * @author Keith James
 * @since 1.2
 */
public class RoundRectangularBeadRenderer extends AbstractBeadRenderer
{
    protected RoundRectangle2D rect;
    protected double arcWidth;
    protected double arcHeight;

    /**
     * Creates a new <code>RoundRectangularBeadRenderer</code>
     * object with the default settings.
     */
    public RoundRectangularBeadRenderer()
    {
        super();
        rect = new RoundRectangle2D.Double();
        arcWidth  = 5.0;
        arcHeight = 5.0;
    }

    /**
     * Creates a new <code>RoundRectangularBeadRenderer</code>.
     *
     * @param beadDepth a <code>double</code>.
     * @param beadDisplacement a <code>double</code>.
     * @param beadOutline a <code>Paint</code>.
     * @param beadFill a <code>Paint</code>.
     * @param beadStroke a <code>Stroke</code>.
     * @param arcWidth a <code>double</code> value which sets the arc
     * width of the corners.
     * @param arcHeight a <code>double</code> value which sets the arc
     * height of the corners.
     */
    public RoundRectangularBeadRenderer(double beadDepth,
                                        double beadDisplacement,
                                        Paint  beadOutline,
                                        Paint  beadFill,
                                        Stroke beadStroke,
                                        double arcWidth,
                                        double arcHeight)
    {
        super(beadDepth, beadDisplacement, beadOutline, beadFill, beadStroke);
        rect = new RoundRectangle2D.Double();
        this.arcWidth  = arcWidth;
        this.arcHeight = arcHeight;
    }

    /**
     * <code>renderBead</code> renders features as a rectangle with
     * rounded corners.
     *
     * @param g2 a <code>Graphics2D</code>.
     * @param f a <code>Feature</code> to render.
     * @param context a <code>SequenceRenderContext</code> context.
     */
    public void renderBead(Graphics2D            g2,
                           Feature               f,
                           SequenceRenderContext context)
    {
        Location loc = f.getLocation();

        int min = loc.getMin();
        int max = loc.getMax();
        int dif = max - min;

        if (context.getDirection() == SequenceRenderContext.HORIZONTAL)
        {
            double  posXW = context.sequenceToGraphics(min);
            double  posYN = beadDisplacement;
            double  width = Math.max(((double) (dif + 1)) * context.getScale(), 1.0);
            double height = Math.min(beadDepth, width / 2.0);

            // If the bead height occupies less than the full height
            // of the renderer, move it down so that it is central
            if (height < beadDepth)
                posYN += ((beadDepth - height) / 2.0);

            rect.setRoundRect(posXW, posYN,
                              Math.floor(width),
                              Math.floor(height),
                              arcWidth, arcHeight);
        }
        else
        {
            double  posXW = beadDisplacement;
            double  posYN = context.sequenceToGraphics(min);
            double height = Math.max(((double) dif + 1) * context.getScale(), 1.0);
            double  width = Math.min(beadDepth, height / 2.0);

            if (width < beadDepth)
                posXW += ((beadDepth - width) /  2.0);

            rect.setRoundRect(posXW, posYN,
                              Math.floor(width),
                              Math.floor(height),
                              arcWidth, arcHeight);
        }

        g2.setPaint(beadFill);
        g2.fill(rect);

        g2.setStroke(beadStroke);
        g2.setPaint(beadOutline);
        g2.draw(rect);
    }

    /**
     * <code>getDepth</code> calculates the depth required by this
     * renderer to display its beads.
     *
     * @param context a <code>SequenceRenderContext</code>.
     *
     * @return a <code>double</code>.
     */
    public double getDepth(SequenceRenderContext context)
    {
        // Get max depth of delegates using base class method
        double maxDepth = super.getDepth(context);
        return Math.max(maxDepth, (beadDepth + beadDisplacement));
    }
}
