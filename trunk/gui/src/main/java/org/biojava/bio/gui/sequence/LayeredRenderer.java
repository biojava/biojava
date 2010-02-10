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
import java.awt.Shape;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.seq.FeatureHolder;

/**
 * <code>LayeredRenderer</code> handles the lane offsets for
 * <code>MultiLineRender</code>s. For each successive lane it
 * translates the <code>Graphics2D</code> perpendicular to the
 * sequence rendering direction by an amount equal to the value
 * returned by the <code>getDepth()</code> method of that lane's
 * renderer.
 *
 * @author Matthew Pocock
 * @author Keith James
 * @since 1.1
 */
public class LayeredRenderer {

    /**
     * Static <code>LayeredRenderer</code> <code>INSTANCE</code> used
     * by <code>MultiLineRenderer</code>s.
     */
    public static final LayeredRenderer INSTANCE = new LayeredRenderer();

    /**
     * <code>getDepth</code> returns the total depth of a list of
     * <code>SequenceRenderer</code>s.
     *
     * @param srcL a <code>List</code> of
     * <code>SequenceRenderContext</code>s.
     * @param renderers a <code>List</code> of
     * <code>SequenceRenderer</code>s.
     *
     * @return a <code>double</code>.
     */
    public double getDepth(List srcL, List renderers) {
        if (srcL.size() != renderers.size()) {
            throw new IllegalArgumentException("srcL and renderers must be the same size: " +
                                               srcL.size() + ":" + renderers.size());
        }

        double depth = 0.0;
        Iterator srcI = srcL.iterator();
        Iterator i = renderers.iterator();

        while (srcI.hasNext() && i.hasNext()) {
            SequenceRenderContext src = (SequenceRenderContext) srcI.next();
            SequenceRenderer sRend = (SequenceRenderer) i.next();            
            if(sRend instanceof OverlayMarker){
                depth += 0.0; // maybe just do nothing here
            }else {
                depth += sRend.getDepth(src);                   
            }
        }
        return depth;
    }

    /**
     * <code>getMinimumLeader</code> returns the maximum value of
     * getMinimumLeader() for a list of <code>SequenceRenderer</code>s.
     *
     * @param srcL a <code>List</code> of
     * <code>SequenceRenderContext</code>s.
     * @param renderers a <code>List</code> of
     * <code>SequenceRenderer</code>s.
     *
     * @return a <code>double</code>.
     */
    public double getMinimumLeader(List srcL, List renderers) {
        if (srcL.size() != renderers.size()) {
            throw new IllegalArgumentException("srcL and renderers must be the same size: " +
                                               srcL.size() + ":" + renderers.size());
        }

        double max = 0.0;
        Iterator srcI = srcL.iterator();
        Iterator i = renderers.iterator();

        while (srcI.hasNext() && i.hasNext()) {
            SequenceRenderContext src = (SequenceRenderContext) srcI.next();
            SequenceRenderer sRend = (SequenceRenderer) i.next();
            max = Math.max(max, sRend.getMinimumLeader(src));
        }
        return max;
    }

    /**
     * <code>getMinimumTrailer</code> returns the maximum value of
     * getMinimumTrailer() for a list of <code>SequenceRenderer</code>s.
     *
     * @param srcL a <code>List</code> of
     * <code>SequenceRenderContext</code>s.
     * @param renderers a <code>List</code> of
     * <code>SequenceRenderer</code>s.
     *
     * @return a <code>double</code>.
     */
    public double getMinimumTrailer(List srcL, List renderers) {
        if (srcL.size() != renderers.size()) {
            throw new IllegalArgumentException("srcL and renderers must be the same size: " +
                                               srcL.size() + ":" + renderers.size());
        }

        double max = 0.0;
        Iterator srcI = srcL.iterator();
        Iterator i = renderers.iterator();

        while (srcI.hasNext() && i.hasNext()) {
            SequenceRenderContext src = (SequenceRenderContext) srcI.next();
            SequenceRenderer sRend = (SequenceRenderer) i.next();
            max = Math.max(max, sRend.getMinimumTrailer(src));
        }
        return max;
    }

    public void paint(Graphics2D g, List srcL, List renderers) {
        if (srcL.size() != renderers.size()) {
            throw new IllegalArgumentException("srcL and renderers must be the same size: " +
                                               srcL.size() + ":" + renderers.size());
        }

        // Offset perpendicular to sequence rendering direction
        double offset = 0.0;
        // Don't know what this is
        double allocatedOffset = 0.0;

        Iterator srcI = srcL.iterator();
        Iterator i = renderers.iterator();

        // New clipping rectangle
        Rectangle2D clip = new Rectangle2D.Double();

        while (srcI.hasNext() && i.hasNext()) {
            SequenceRenderContext src = (SequenceRenderContext) srcI.next();
            SequenceRenderer sRend = (SequenceRenderer) i.next();
            int dir = src.getDirection();
            double depth = sRend.getDepth(src);

            // Sequence range should be inclusive of the min/max
            // coordinates for sequenceToGraphics() so we use
            // src.getRange().getMax() + 1
            double minP = src.sequenceToGraphics(src.getRange().getMin()) -
                sRend.getMinimumLeader(src);
            double maxP = src.sequenceToGraphics(src.getRange().getMax() + 1) +
                sRend.getMinimumTrailer(src);

            // Added +1 to these as the outer edge of features was
            // being clipped off
            if (dir == SequenceRenderContext.HORIZONTAL) {
                clip.setFrame(minP, 0.0, maxP - minP + 1, depth + 1);
                g.translate(0.0, offset);
            } else {
                clip.setFrame(0.0, minP, depth + 1, maxP - minP + 1);
                g.translate(offset, 0.0);
            }

            Shape oldClip = g.getClip();
            g.clip(clip);
            sRend.paint(g, src);
            g.setClip(oldClip);

            if (dir == SequenceRenderContext.HORIZONTAL) {
                g.translate(0.0, -offset);
            } else {
                g.translate(-offset, 0.0);
            }

            if (sRend instanceof OverlayMarker)  {
                // overlay, just record maximum allocation
                allocatedOffset = Math.max(allocatedOffset, sRend.getDepth(src));
            } else {
                // non-overlaid: apply all relevant offsets
                offset += Math.max(sRend.getDepth(src), allocatedOffset);
                allocatedOffset = 0.0;  // clear it as it is now applied.
            }
        }
    }

    public SequenceViewerEvent processMouseEvent(List srcL, MouseEvent me,
                                                 List path, List renderers) {
        if (srcL.size() != renderers.size()) {
            throw new IllegalArgumentException("srcL and renderers must be the same size: " +
                                               srcL.size() + ":" + renderers.size());
        }

        // Offset perpendicular to sequence rendering direction
        double offset = 0.0;

        Iterator srcI = srcL.iterator();
        Iterator i = renderers.iterator();

        SequenceViewerEvent sve = null;
        
        while (srcI.hasNext() && i.hasNext()) {
            SequenceRenderContext src = (SequenceRenderContext) srcI.next();
            SequenceRenderer sRend = (SequenceRenderer) i.next();
            double depth = sRend.getDepth(src);

            SequenceViewerEvent thisSve = null;
            if (src.getDirection() == SequenceRenderContext.HORIZONTAL) {
                if ((me.getY() >= offset) && (me.getY() < offset + depth)) {
                    me.translatePoint(0, (int) -offset);
                    thisSve = sRend.processMouseEvent(src, me, path);
                    me.translatePoint(0, (int) +offset);
                }
            } else {
                if ((me.getX() >= offset) && (me.getX() < offset + depth)) {
                    me.translatePoint((int) -offset, 0);
                    thisSve = sRend.processMouseEvent(src, me, path);
                    me.translatePoint((int) +offset, 0);
                }
            }

            if (thisSve != null) {
                if (sve == null) {
                    sve = thisSve;
                } else if (thisSve.getTarget() instanceof FeatureHolder &&
                           ((FeatureHolder) thisSve.getTarget()).countFeatures() > 0) 
                {
                    // features trump anything else
                    sve = thisSve;
                }
            }

            if (! (sRend instanceof OverlayMarker)) offset += depth;
        }
        return sve;
    }
}

