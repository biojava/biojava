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

import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;


/**
 * OffsetRulerRenderer can render the ruler starting from an arbitrary offset from the sequence.
 * For example if the Protein contained an N-Terminal His tag then coordinate 1 should correspond
 * to the start of the protein, not the tag.  This implementation borrows heavily from
 * RulerRenderer
 *
 * @author Mark Southern
 *
 * @since 1.5
 */
public class OffsetRulerRenderer extends AbstractChangeable implements SequenceRenderer {
    public static final ChangeType OFFSET = new ChangeType("The ruler offset has changed",
            "org.biojava.bio.gui.sequence.OffsetRulerRenderer", "OFFSET",
            SequenceRenderContext.REPAINT
        );
    public static final ChangeType TICKS = new ChangeType("The ruler tick direction has changed",
            "org.biojava.bio.gui.sequence.OffsetRulerRenderer", "TICKS",
            SequenceRenderContext.REPAINT
        );
    public static final int TICKS_UP = 0;
    public static final int TICKS_DOWN = 1;
    private Line2D line;
    private double depth;
    private AffineTransform antiQuarter;
    private int tickDirection;
    private float tickHeight;
    private float horizLabelOffset;
    private float vertLabelOffset;
    private int sequenceOffset = 0;

    public OffsetRulerRenderer() throws IllegalArgumentException {
        this(TICKS_DOWN, 0);
    }

    public OffsetRulerRenderer(int tickDirection, int sequenceOffset)
        throws IllegalArgumentException {
        this.sequenceOffset = sequenceOffset;

        line = new Line2D.Double();
        antiQuarter = AffineTransform.getRotateInstance(Math.toRadians(-90));

        if ((tickDirection == TICKS_UP) || (tickDirection == TICKS_DOWN)) {
            this.tickDirection = tickDirection;
        } else {
            throw new IllegalArgumentException(
                "Tick direction may only be set to RulerRenderer.TICKS_UP or RulerRenderer.TICKS_DOWN"
            );
        }

        depth = 20.0;
        tickHeight = 4.0f;

        horizLabelOffset = (( float ) depth) - tickHeight - 2.0f;
        vertLabelOffset = (( float ) depth) - ((tickHeight + 2.0f) * 2.0f);
    }

    public void setSequenceOffset(int offset) throws ChangeVetoException {
        if (hasListeners()) {
            ChangeSupport cs = getChangeSupport(SequenceRenderContext.REPAINT);

            synchronized (cs) {
                ChangeEvent ce = new ChangeEvent(this, OFFSET);
                cs.firePreChangeEvent(ce);
                this.sequenceOffset = offset;
                cs.firePostChangeEvent(ce);
            }
        } else {
            this.sequenceOffset = offset;
        }
    }

    public int getSequenceOffset() {
        return this.sequenceOffset;
    }

    public void setTickDirection(int dir) throws ChangeVetoException {
        if (hasListeners()) {
            ChangeSupport cs = getChangeSupport(SequenceRenderContext.REPAINT);

            synchronized (cs) {
                ChangeEvent ce = new ChangeEvent(this, TICKS);
                cs.firePreChangeEvent(ce);
                tickDirection = dir;
                cs.firePostChangeEvent(ce);
            }
        } else {
            tickDirection = dir;
        }
    }

    public int getTickDirection() {
        return tickDirection;
    }

    public double getMinimumLeader(SequenceRenderContext context) {
        return 0.0;
    }

    public double getMinimumTrailer(SequenceRenderContext context) {
        return 0.0;
    }

    public double getDepth(SequenceRenderContext src) {
        return depth + 1.0;
    }

    public void paint(Graphics2D g2, SequenceRenderContext context) {
        g2.setStroke(new java.awt.BasicStroke(1F));

        AffineTransform prevTransform = g2.getTransform();

        g2.setPaint(Color.black);

        int min = context.getRange().getMin();
        int max = context.getRange().getMax();
        double minX = context.sequenceToGraphics(min);
        double maxX = context.sequenceToGraphics(max);
        double scale = context.getScale();

        double halfScale = scale * 0.5;

        if (context.getDirection() == SequenceRenderContext.HORIZONTAL) {
            if (tickDirection == TICKS_UP) {
                line.setLine(minX - halfScale, depth, maxX + halfScale, depth);
            } else {
                line.setLine(minX - halfScale, 0.0, maxX + halfScale, 0.0);
            }
        } else {
            if (tickDirection == TICKS_UP) {
                line.setLine(depth, minX - halfScale, depth, maxX + halfScale);
            } else {
                line.setLine(0.0, minX - halfScale, 0.0, maxX + halfScale);
            }
        }

        g2.draw(line);

        FontMetrics fMetrics = g2.getFontMetrics();

        // The widest (== maxiumum) coordinate to draw
        int coordWidth = fMetrics.stringWidth(Integer.toString(max));

        // Minimum gap getween ticks
        double minGap = ( double ) Math.max(coordWidth, 40);

        // How many symbols does a gap represent?
        int realSymsPerGap = ( int ) Math.ceil(((minGap + 5.0) / context.getScale()));

        // We need to snap to a value beginning 1, 2 or 5.
        double exponent = Math.floor(Math.log(realSymsPerGap) / Math.log(10));
        double characteristic = realSymsPerGap / Math.pow(10.0, exponent);

        int snapSymsPerGap;

        if (characteristic > 5.0) {
            // Use unit ticks
            snapSymsPerGap = ( int ) Math.pow(10.0, exponent + 1.0);
        } else if (characteristic > 2.0) {
            // Use ticks of 5
            snapSymsPerGap = ( int ) (5.0 * Math.pow(10.0, exponent));
        } else {
            snapSymsPerGap = ( int ) (2.0 * Math.pow(10.0, exponent));
        }

        min -= Math.abs(sequenceOffset);
        max += Math.abs(sequenceOffset);

        int minP = min + ((snapSymsPerGap - min) % snapSymsPerGap);

        for (int index = minP; index <= max; index += snapSymsPerGap) {
            double offset = context.sequenceToGraphics(index + sequenceOffset);
            String labelString = String.valueOf(index); // + sequenceOffset);
            float halfLabelWidth = fMetrics.stringWidth(labelString) / 2;

            if (context.getDirection() == SequenceRenderContext.HORIZONTAL) {
                if (tickDirection == TICKS_UP) {
                    line.setLine(offset + halfScale, depth - tickHeight, offset + halfScale, depth);
                    g2.drawString(labelString, ( float ) ((offset + halfScale) - halfLabelWidth),
                        horizLabelOffset
                    );
                } else {
                    line.setLine(offset + halfScale, 0.0, offset + halfScale, tickHeight);
                    g2.drawString(labelString, ( float ) ((offset + halfScale) - halfLabelWidth),
                        horizLabelOffset
                    );
                }
            } else { // vertical

                if (tickDirection == TICKS_UP) {
                    line.setLine(depth, offset + halfScale, depth - tickHeight, offset + halfScale);
                    g2.translate(vertLabelOffset, offset + halfScale + halfLabelWidth);
                    g2.transform(antiQuarter);
                    g2.drawString(labelString, 0.0f, 0.0f);
                    g2.setTransform(prevTransform);
                } else {
                    line.setLine(0.0f, offset + halfScale, tickHeight, offset + halfScale);
                    g2.translate(vertLabelOffset, offset + halfScale + halfLabelWidth);
                    g2.transform(antiQuarter);
                    g2.drawString(labelString, 0.0f, 0.0f);
                    g2.setTransform(prevTransform);
                }
            }

            g2.draw(line);
        }
    }

    public SequenceViewerEvent processMouseEvent(SequenceRenderContext context, MouseEvent me,
        List path
    ) {
        path.add(this);

        int sPos = context.graphicsToSequence(me.getPoint());

        return new SequenceViewerEvent(this, null, sPos, me, path);
    }
}
