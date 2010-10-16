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
import java.awt.geom.Rectangle2D;

import org.biojava.bio.seq.Feature;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * <p><code>RectangularBeadRenderer</code> renders features as simple
 * rectangles. Their outline and fill <code>Paint</code>,
 * <code>Stroke</code>, feature depth, Y-axis displacement are
 * configurable. The height of the rectangle will be equal to half its
 * width, but not greater than the <code>beadDepth</code> set in the
 * constructor.</p>
 *
 * <p>An alternative bead height behaviour is available where the
 * rectangle height does not scale with its current width. The
 * <code>setHeightScaling</code> method should be passed a boolean
 * value to change this. The default is to use height scaling.</p>
 *
 * @author Keith James
 * @since 1.2
 */
public class RectangularBeadRenderer extends AbstractBeadRenderer
{
    /**
     * Constant <code>HEIGHTSCALING</code> indicating a change to the
     * feature height scaling policy.
     */
    public static final ChangeType HEIGHTSCALING =
	new ChangeType("The height scaling policy of the features has changed",
		       "org.biojava.bio.gui.sequence.RectangularBeadRenderer",
		       "HEIGHTSCALING", SequenceRenderContext.LAYOUT);

    protected Rectangle2D rect;
    protected boolean scaleHeight;

    /**
     * Creates a new <code>RectangularBeadRenderer</code> with the
     * default settings.
     */
    public RectangularBeadRenderer()
    {
        super();
        rect = new Rectangle2D.Double();
        scaleHeight = true;
    }

    /**
     * Creates a new <code>RectangularBeadRenderer</code>.
     *
     * @param beadDepth a <code>double</code>.
     * @param beadDisplacement a <code>double</code>.
     * @param beadOutline a <code>Paint</code>.
     * @param beadFill a <code>Paint</code>.
     * @param beadStroke a <code>Stroke</code>.
     */
    public RectangularBeadRenderer(double beadDepth,
                                   double beadDisplacement,
                                   Paint  beadOutline,
                                   Paint  beadFill,
                                   Stroke beadStroke)
    {
        super(beadDepth, beadDisplacement, beadOutline, beadFill, beadStroke);
        rect = new Rectangle2D.Double();
        scaleHeight = true;
    }

    /**
     * <code>renderBead</code> renders features as simple rectangle.
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

        double posXW, posYN, width, height;

        if (context.getDirection() == SequenceRenderContext.HORIZONTAL)
        {
            posXW  = context.sequenceToGraphics(min);
            posYN  = beadDisplacement;
            width  = Math.max(((double) (dif + 1)) * context.getScale(), 1.0);

            if (scaleHeight)
            {
                height = Math.min(beadDepth, width / 2.0);

                // If the bead height occupies less than the full height
                // of the renderer, move it down so that it is central
                if (height < beadDepth)
                    posYN += ((beadDepth - height) / 2.0);
            }
            else
            {
                height = beadDepth;
            }

            rect.setRect(posXW, posYN,
                         Math.floor(width),
                         Math.floor(height));
        }
        else
        {
            posXW  = beadDisplacement;
            posYN  = context.sequenceToGraphics(min);
            height = Math.max(((double) dif + 1) * context.getScale(), 1.0);

            if (scaleHeight)
            {
                width = Math.min(beadDepth, height / 2.0);

                if (width < beadDepth)
                    posXW += ((beadDepth - width) /  2.0);
            }
            else
            {
                width = beadDepth;
            }

            rect.setRect(posXW, posYN,
                         Math.floor(width),
                         Math.floor(height));
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

    /**
     * <code>getHeightScaling</code> returns the state of the height
     * scaling policy.
     *
     * @return a <code>boolean</code> true if height scaling is
     * enabled.
     */
    public boolean getHeightScaling()
    {
        return scaleHeight;
    }

    /**
     * <code>setHeightScaling</code> sets the height scaling
     * policy. Default behaviour is for this to be enabled leading to
     * features being drawn with a height equal to half their width,
     * subject to a maximum height restriction equal to the
     * <code>beadDepth</code> set in the constructor. If disabled,
     * features will always be drawn at the maximum height allowed by
     * the <code>beadDepth</code> parameter.
     *
     * @param isEnabled a <code>boolean</code>.
     *
     * @exception ChangeVetoException if an error occurs.
     */
    public void setHeightScaling(boolean isEnabled) throws ChangeVetoException
    {
        if (hasListeners())
	{
	    ChangeSupport cs = getChangeSupport(SequenceRenderContext.LAYOUT);
	    synchronized(cs)
	    {
		ChangeEvent ce = new ChangeEvent(this, SequenceRenderContext.LAYOUT,
						 null, null,
						 new ChangeEvent(this, HEIGHTSCALING,
								 new Boolean(scaleHeight),
								 new Boolean(isEnabled)));
		cs.firePreChangeEvent(ce);
                scaleHeight = isEnabled;
		cs.firePostChangeEvent(ce);
	    }
	}
	else
	{
            scaleHeight = isEnabled;
	}
    }
}
