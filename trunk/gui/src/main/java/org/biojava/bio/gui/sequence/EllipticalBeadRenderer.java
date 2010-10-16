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
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;

import org.biojava.bio.seq.Feature;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * <p><code>EllipticalBeadRenderer</code> renders features as simple
 * ellipses. Their outline and fill <code>Paint</code>,
 * <code>Stroke</code>, feature depth, Y-axis displacement are
 * configurable. Also configurable is the minimum ratio of long axis
 * to short axis of the ellipse - this prevents long features also
 * becoming ever wider and obscuring neighbours.</p>
 *
 * @author Keith James
 * @since 1.2
 */
public class EllipticalBeadRenderer extends AbstractBeadRenderer
{
    /**
     * Constant <code>RATIO</code> indicating a change to the minimum
     * allowed ratio of long axis to short axis of the features.
     */
    public static final ChangeType RATIO =
	new ChangeType("The shape of the features has changed",
		       "org.biojava.bio.gui.sequence.EllipticalBeadRenderer",
		       "RATIO", SequenceRenderContext.LAYOUT);

    protected double dimensionRatio;

    /**
     * Creates a new <code>EllipticalBeadRenderer</code> object
     * with the default settings.
     */
    public EllipticalBeadRenderer()
    {
	super();
	dimensionRatio = 2.0F;
    }

    /**
     * Creates a new <code>EllipticalBeadRenderer</code>.
     *
     * @param beadDepth a <code>double</code>.
     * @param beadDisplacement a <code>double</code>.
     * @param beadOutline a <code>Paint</code>.
     * @param beadFill a <code>Paint</code>.
     * @param beadStroke a <code>Stroke</code>.
     * @param dimensionRatio a <code>double</code>.
     */
    public EllipticalBeadRenderer(double beadDepth,
				  double beadDisplacement,
				  Paint  beadOutline,
				  Paint  beadFill,
				  Stroke beadStroke,
				  double dimensionRatio)
    {
	super(beadDepth, beadDisplacement, beadOutline, beadFill, beadStroke);
	dimensionRatio = 2.0F;
    }

    /**
     * <code>renderBead</code> renders features as simple ellipse.
     *
     * @param g2 a <code>Graphics2D</code> context.
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

	Shape shape;

	if (context.getDirection() == SequenceRenderContext.HORIZONTAL)
	{
	    double  posXW = context.sequenceToGraphics(min);
	    double  posYN = beadDisplacement;
	    double  width = Math.max(((double) (dif + 1)) * context.getScale(), 1.0f);
	    double height = Math.min(beadDepth, width / dimensionRatio);

	    // If the bead height occupies less than the full height
	    // of the renderer, move it down so that it is central
	    if (height < beadDepth)
		posYN += ((beadDepth - height) / dimensionRatio);

	    shape = new Ellipse2D.Double(posXW, posYN, width, height);
	}
	else
	{
	    double  posXW = beadDisplacement;
	    double  posYN = context.sequenceToGraphics(min);
	    double height = Math.max(((double) dif + 1) * context.getScale(), 1.0f);
	    double  width = Math.min(beadDepth, height / dimensionRatio);

	    if (width < beadDepth)
		posXW += ((beadDepth - width) /  dimensionRatio);

	    shape = new Ellipse2D.Double(posXW, posYN, width, height);
	}

	g2.setPaint(beadFill);
	g2.fill(shape);

	g2.setStroke(beadStroke);
	g2.setPaint(beadOutline);
	g2.draw(shape);
    }

    /**
     * <code>getDepth</code> calculates the depth required by this
     * renderer to display its beads.
     *
     * @param context a <code>SequenceRenderContext</code> object.
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
     * <code>getDimensionRatio</code> returns the maximum ratio of
     * long dimension to short dimension of the bead. This should be
     * equal, or greater than 1.
     *
     * @return a <code>double</code>.
     */
    public double getDimensionRatio()
    {
	return dimensionRatio;
    }

    /**
     * <code>setDimensionRatio</code> sets the minimum ratio of
     * long dimension to short dimension of the bead. This should be
     * equal, or greater than 1.
     *
     * @param ratio a <code>double</code> ratio of depth.
     *
     * @exception ChangeVetoException if an error occurs.
     */
    public void setDimensionRatio(double ratio) throws ChangeVetoException
    {
	if (ratio < 1.0F)
	    throw new ChangeVetoException("The long dimension may not be less than the short dimension (ratio >= 1.0)");

	if (hasListeners())
	{
	    ChangeSupport cs = getChangeSupport(SequenceRenderContext.LAYOUT);
	    synchronized(cs)
	    {
		ChangeEvent ce = new ChangeEvent(this, SequenceRenderContext.LAYOUT,
						 null, null,
						 new ChangeEvent(this, RATIO,
								 new Double(dimensionRatio),
								 new Double(ratio)));
		cs.firePreChangeEvent(ce);
		dimensionRatio= ratio;
		cs.firePostChangeEvent(ce);
	    }
	}
	else
	{
	    dimensionRatio = ratio;
	}
    }
}
