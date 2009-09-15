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
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.util.List;

import org.biojava.bio.symbol.Location;

/**
 * <p><code>RulerRenderer</code> renders numerical scales in sequence
 * coordinates. The tick direction may be set to point upwards (or
 * left when the scale is vertical) or downwards (right when the scale
 * is vertical).</p>
 *
 * <p>Note: The Compaq Java VMs 1.3.1 - 1.4.0 on Tru64 appear to have
 * a bug in font transformation which prevents a vertically oriented
 * ruler displaying correctly rotated text.</p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author David Huen
 * @author Keith James
 * @author Kalle Näslund
 */
public class RulerRenderer implements SequenceRenderer
{
    /**
     * <code>TICKS_UP</code> indicates that the ticks will point
     * upwards from a baseline.
     */
    public static final int TICKS_UP = 0;
    /**
     * <code>TICKS_DOWN</code> indicates that the ticks will point
     * downwards from a baseline.
     */
    public static final int TICKS_DOWN = 1;

    private Line2D            line;
    private double            depth;
    private AffineTransform   antiQuarter;
    private int               tickDirection;
    private float             tickHeight;
    private float             horizLabelOffset;
    private float             vertLabelOffset;

    /**
     * Creates a new <code>RulerRenderer</code> with the default
     * setting of ticks pointing downwards.
     */
    public RulerRenderer() throws IllegalArgumentException
    {
        this(TICKS_DOWN);
    }

    /**
     * Creates a new <code>RulerRenderer</code> with the specified
     * tick direction.
     *
     * @param tickDirection an <code>int</code>.
     * @exception IllegalArgumentException if an error occurs.
     */
    public RulerRenderer(int tickDirection) throws IllegalArgumentException
    {
        line   = new Line2D.Double();
        antiQuarter = AffineTransform.getRotateInstance(Math.toRadians(-90));

        if (tickDirection == TICKS_UP || tickDirection == TICKS_DOWN)
            this.tickDirection = tickDirection;
        else
            throw new IllegalArgumentException("Tick direction may only be set to RulerRenderer.TICKS_UP or RulerRenderer.TICKS_DOWN");

        depth      = 20.0;
        tickHeight = 4.0f;

        horizLabelOffset = ((float) depth) - tickHeight - 2.0f;
        vertLabelOffset  = ((float) depth) - ((tickHeight + 2.0f) * 2.0f);
    }

    public double getMinimumLeader(SequenceRenderContext context)
    {
        return 0.0;
    }

    public double getMinimumTrailer(SequenceRenderContext context)
    {
        return 0.0;
    }

    public double getDepth(SequenceRenderContext src)
    {
        return depth + 1.0;
    }

    public void paint(Graphics2D g2, SequenceRenderContext context)
    {
        AffineTransform prevTransform = g2.getTransform();

        g2.setPaint(Color.black);

        Location visible = GUITools.getVisibleRange(context, g2);
        if( visible == Location.empty ) {
            return;
        }
        
        int min = visible.getMin();
        int max = visible.getMax();
        double minX = context.sequenceToGraphics(min);
        double maxX = context.sequenceToGraphics(max);
        double scale = context.getScale();

        double halfScale = scale * 0.5;

        if (context.getDirection() == SequenceRenderContext.HORIZONTAL)
        {
            if (tickDirection == TICKS_UP)
            {
                line.setLine(minX - halfScale, depth,
                             maxX + halfScale, depth);
            }
            else
            {
                line.setLine(minX - halfScale, 0.0,
                             maxX + halfScale, 0.0);
            }
        }
        else
        {
            if (tickDirection == TICKS_UP)
            {
                line.setLine(depth, minX - halfScale,
                             depth, maxX + halfScale);
            }
            else
            {
                line.setLine(0.0, minX - halfScale,
                             0.0, maxX + halfScale);
            }
        }

        g2.draw(line);

        FontMetrics fMetrics = g2.getFontMetrics();

        // The widest (== maxiumum) coordinate to draw
        int coordWidth = fMetrics.stringWidth(Integer.toString(max));

        // Minimum gap getween ticks
        double minGap = (double) Math.max(coordWidth, 40);

        // How many symbols does a gap represent?
        int realSymsPerGap = (int) Math.ceil(((minGap + 5.0) / context.getScale()));

        // We need to snap to a value beginning 1, 2 or 5.
        double exponent = Math.floor(Math.log(realSymsPerGap) / Math.log(10));
        double characteristic = realSymsPerGap / Math.pow(10.0, exponent);

        int snapSymsPerGap;
        if (characteristic > 5.0)
        {
            // Use unit ticks
            snapSymsPerGap = (int) Math.pow(10.0, exponent + 1.0);
        }
        else if (characteristic > 2.0)
        {
            // Use ticks of 5
            snapSymsPerGap = (int) (5.0 * Math.pow(10.0, exponent));
        }
        else
        {
            snapSymsPerGap = (int) (2.0 * Math.pow(10.0, exponent));
        }

        int minP = min + (snapSymsPerGap - min) % snapSymsPerGap;

        for (int index = minP; index <= max; index += snapSymsPerGap)
        {
            double offset = context.sequenceToGraphics(index);
            String labelString = String.valueOf(index);
            float halfLabelWidth = fMetrics.stringWidth(labelString) / 2;

            if (context.getDirection() == SequenceRenderContext.HORIZONTAL)
            {
                if (tickDirection == TICKS_UP)
                {
                    line.setLine(offset + halfScale, depth - tickHeight,
                                 offset + halfScale, depth);
                    g2.drawString(String.valueOf(index),
                                  (float) (offset + halfScale - halfLabelWidth),
                                  horizLabelOffset);
                }
                else
                {
                    line.setLine(offset + halfScale, 0.0,
                                 offset + halfScale, tickHeight);
                    g2.drawString(String.valueOf(index),
                                  (float) (offset + halfScale - halfLabelWidth),
                                  horizLabelOffset);
                }
            }
            else
            {
                if (tickDirection == TICKS_UP)
                {
                    line.setLine(depth, offset + halfScale,
                                 depth - tickHeight, offset + halfScale);
                    g2.translate(vertLabelOffset,
                                 offset + halfScale + halfLabelWidth);
                    g2.transform(antiQuarter);
                    g2.drawString(String.valueOf(index), 0.0f, 0.0f);
                    g2.setTransform(prevTransform);
                }
                else
                {
                    line.setLine(0.0f, offset + halfScale,
                                 tickHeight, offset + halfScale);
                    g2.translate(vertLabelOffset,
                                 offset + halfScale + halfLabelWidth);
                    g2.transform(antiQuarter);
                    g2.drawString(String.valueOf(index), 0.0f, 0.0f);
                    g2.setTransform(prevTransform);
                }
            }
            g2.draw(line);
        }
    }

    public SequenceViewerEvent processMouseEvent(SequenceRenderContext context,
                                                 MouseEvent            me,
                                                 List                  path)
    {
        path.add(this);
        int sPos = context.graphicsToSequence(me.getPoint());
        return new SequenceViewerEvent(this, null, sPos, me, path);
    }
}
