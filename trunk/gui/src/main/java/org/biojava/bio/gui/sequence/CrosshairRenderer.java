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
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.geom.Line2D;
import java.util.List;

import javax.swing.SwingUtilities;

import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * <p>
 * <code>CrosshairRenderer</code> draws a crosshair, optionally
 * with coordinates. The crosshair is set to a sequence position by a
 * click and then stays there through scrolls/rescales until the next
 * click. See the <code>processMouseEvent</code> documentation for
 * details of responses to various mouse actions.
 * </p>
 *
 * @author Keith James
 * @since 1.2
 */
public class CrosshairRenderer extends AbstractChangeable
    implements PairwiseSequenceRenderer 
{
    /**
     * Constant <code>OUTLINE</code> indicating a change to the
     * crosshair paint.
     */
    public static final ChangeType OUTLINE =
        new ChangeType("The outline paint has changed",
                       "org.biojava.bio.gui.sequence.CrosshairRenderer",
                       "OUTLINE", SequenceRenderContext.REPAINT);
    
    /**
     * <code>xHair</code> is the vertical line positioned along the
     * X-axis.
     */
    protected Line2D xHair;

    /**
     * <code>yHair</code> is the horizontal line positioned along the
     * Y-axis.
     */
    protected Line2D yHair;
    /**
     * <code>point</code> is the current location (in sequence
     * coordinates) of the crosshair in the X and Y sequences.
     */
    protected Point point;

    // Crosshair location in sequence coordinates
    private int sPosX, sPosY;
    // Crosshair colour
    private Paint outline;
    // Display coordinates?
    private boolean display;

    /**
     * Creates a new <code>CrosshairRenderer</code> in light grey with
     * coordinates displayed.
     */
    public CrosshairRenderer()
    {
        this(Color.lightGray);
    }

    /**
     * Creates a new <code>CrosshairRenderer</code> of the specified
     * colour, with coordinates displayed.
     *
     * @param outline a <code>Paint</code>.
     */
    public CrosshairRenderer(Paint outline)
    {
        xHair = new Line2D.Double();
        yHair = new Line2D.Double();
        point = new Point();

        sPosX = 1;
        sPosY = 1;

        display = true;

        this.outline = outline;
    }

    public void paint(Graphics2D g2, PairwiseRenderContext context)
    {
        Rectangle clip = g2.getClipBounds();

        double xMin = clip.getMinX();
        double xMax = clip.getMaxX();
        double yMin = clip.getMinY();
        double yMax = clip.getMaxY();

        // Offset to get the hair to line up with the centre of a
        // residue
        double residueCentre = context.getScale() * 0.5;

        double gPosX = context.sequenceToGraphics(sPosX);
        gPosX += residueCentre;

        double gPosY = context.secondarySequenceToGraphics(sPosY);
        gPosY += residueCentre;

        if (context.getDirection() == PairwiseRenderContext.HORIZONTAL)
        {
            xHair.setLine(gPosX, yMin, gPosX, yMax);
            yHair.setLine(xMin, gPosY, xMax, gPosY);
        }
        else
        {
            xHair.setLine(xMin, gPosY, xMax, gPosY);
            yHair.setLine(gPosX, yMin, gPosX, yMax);
        }

        g2.setPaint(outline);
        g2.draw(xHair);
        g2.draw(yHair);

        if (display)
        {
            g2.setFont(context.getFont());
            g2.drawString(sPosX + ", " + sPosY,
                          (float) (gPosX + 5.0),
                          (float) (gPosY - 5.0));
        }
    }

    /**
     * <code>coordinateDisplayOn</code> toggles the display of
     * sequence coordinates.
     *
     * @param display a <code>boolean</code>.
     */
    public void coordinateDisplayOn(boolean display)
    {
        this.display = display;
    }

    /**
     * <code>getOutline</code> returns the colour used to draw the
     * lines.
     *
     * @return a <code>Paint</code>.
     */
    public Paint getOutline()
    {
        return outline;
    }

    /**
     * <code>setOutline</code> sets the the colour used to draw the
     * lines.
     *
     * @param outline a <code>Paint</code>.
     *
     * @exception ChangeVetoException if an error occurs.
     */
    public void setOutline(Paint outline) throws ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(SequenceRenderContext.REPAINT);
            synchronized(cs)
            {
                ChangeEvent ce = new ChangeEvent(this, SequenceRenderContext.REPAINT,
                                                 null, null,
                                                 new ChangeEvent(this, OUTLINE,
                                                                 outline,
                                                                 this.outline));
                cs.firePreChangeEvent(ce);
                this.outline = outline;
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            this.outline = outline;
        }
    }

    /**
     * <p><code>processMouseEvent</code> processes any
     * <code>MouseEvent</code>s directed to the renderer.</p>
     *
     * <p>
     * Mouse actions are as follows (all are button-1 only):
     * <ul>
     *  <li>Click sets the crosshair position and returns an event
     *      whose target object is the <code>Point</code> in sequence
     *      coordinates. The X coordinate is in the primary sequence,
     *      the Y coordinate is in the secondary sequence.</li>
     *  <li>Press same as Click</li>
     *  <li>Drag same as Click, except that the <code>Point</code> is
     *      not set</li>
     *  <li>Release same as Click, except that the <code>Point</code>
     *      is not set</li>
     *  <li>Move same as Click, except that the <code>Point</code> is
     *      not set and the target is null</li>
     * </ul>
     * </p>
     *
     * @param context a <code>PairwiseRenderContext</code>.
     * @param me a <code>MouseEvent</code>.
     * @param path a <code>List</code>.
     *
     * @return a <code>SequenceViewerEvent</code>.
     */
    public SequenceViewerEvent processMouseEvent(PairwiseRenderContext context,
                                                 MouseEvent            me,
                                                 List                  path)
    {
        path.add(this);

        // Only hit the point with left button
        if (! SwingUtilities.isLeftMouseButton(me))
            return new SequenceViewerEvent(this, null, sPosX, me, path);

        // Only hit the point with click/release/drag
        int id = me.getID();
        if (! (id == MouseEvent.MOUSE_CLICKED ||
               id == MouseEvent.MOUSE_PRESSED ||
               id == MouseEvent.MOUSE_DRAGGED ||
               id == MouseEvent.MOUSE_RELEASED))
            return new SequenceViewerEvent(this, null, sPosX, me, path);

        // Only move the point on clicks or presses
        if (id == MouseEvent.MOUSE_CLICKED ||
            id == MouseEvent.MOUSE_PRESSED)
        {
            double gPosX, gPosY;

            if (context.getDirection() == PairwiseRenderContext.HORIZONTAL)
            {
                gPosX = me.getPoint().getX();
                gPosY = me.getPoint().getY();
            }
            else
            {
                gPosX = me.getPoint().getY();
                gPosY = me.getPoint().getX();
            }

            sPosX = context.graphicsToSequence(gPosX);
            sPosY = context.graphicsToSecondarySequence(gPosY);

            // We were clicked, so set the point here
            point.setLocation(sPosX, sPosY);

            // Inform any REPAINT listeners that they need to repaint this
            if (hasListeners())
            {
                ChangeSupport cs = getChangeSupport(SequenceRenderContext.REPAINT);
                synchronized(cs)
                {
                    ChangeEvent ce = new ChangeEvent(this, SequenceRenderContext.REPAINT);
                    cs.firePostChangeEvent(ce);
                }
            }
        }

        return new SequenceViewerEvent(this, point, sPosX, me, path);
    }
}
