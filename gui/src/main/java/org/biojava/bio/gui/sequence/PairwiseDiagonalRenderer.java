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
import java.awt.event.MouseEvent;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.io.Serializable;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.BioError;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.biojava.bio.seq.homol.SimilarityPairFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * <p><code>PairwiseDiagonalRenderer</code> renders a region of
 * similarity between two sequences as a straight line. The effect
 * produced is similar to a dotplot. This implementation requires that
 * these regions be represented by
 * <code>SimilarityPairFeature</code>s.<p>
 *
 * <p>Drawing outside the visible area using a range of valid
 * <code>double</code>s may cause Java to hang (Sun JDK 1.3.1 on
 * Linux, Compaq JDK 1.3.1 on Tru64, but not Sun JDK 1.4.0-beta2-b77
 * on Linux). I got round this by manual clipping of the lines to the
 * clip area. The code uses an implementation of the Cohen-Sutherland
 * line-clipping algorithm which clips lines to within a
 * rectangle.</p>
 *
 * <p>The clipping code is taken from Computer Graphics for Java
 * Programmers by Leen Ammeraal (1998, ISBN 0-471-98142-7) and
 * cosmetically altered to support Java2D objects. Any bugs introduced
 * are my responsibility.</p>
 *
 * @author Keith James
 * @author Leen Ammeraal
 * @since 1.2
 */
public class PairwiseDiagonalRenderer extends AbstractChangeable
    implements PairwiseSequenceRenderer, Serializable
{
    /**
     * Constant <code>OUTLINE</code> indicating a change to the fill of
     * the features.
     */
    public static final ChangeType OUTLINE =
        new ChangeType("The outline paint has changed",
                       "org.biojava.bio.gui.sequence.PairwiseDiagonalRenderer",
                       "OUTLINE", SequenceRenderContext.REPAINT);

    /**
     * <code>spf</code> is a filter which excludes all features except
     * <code>SimilarityPairFeature</code>s.
     */
    private static FeatureFilter spf;

    static
    {
        String className =
            "org.biojava.bio.seq.homol.SimilarityPairFeature";

        try
        {
            spf = new FeatureFilter.ByClass(Class.forName(className));
        }
        catch (Exception e)
        {
            throw new BioError("Failed to load Class for " + className, e);
        }
    }

    /**
     * <code>line</code> is the line to be drawn for each feature.
     */
    protected Line2D.Float line;

    /**
     * <code>outline</code> is the line colour.
     */
    protected Paint outline;

    /**
     * Creates a new <code>PairwiseDiagonalRenderer</code> which will
     * draw black lines.
     */
    public PairwiseDiagonalRenderer()
    {
        this(Color.black);
    }

    /**
     * Creates a new <code>PairwiseDiagonalRenderer</code> which will
     * draw lines using the specified <code>Paint</code>.
     *
     * @param outline a <code>Paint</code>.
     */
    public PairwiseDiagonalRenderer(Paint outline)
    {
        line = new Line2D.Float();
        this.outline = outline;
    }

    /**
     * <code>paint</code> renders the feature as a simple line.
     *
     * @param g2 a <code>Graphics2D</code>.
     * @param context a <code>PairwiseRenderContext</code>.
     */
    public void paint(Graphics2D g2, PairwiseRenderContext context)
    {
        FeatureHolder fh;

        if (context.getDirection() == PairwiseRenderContext.HORIZONTAL)
        {
            fh = context.getFeatures().filter(new
                FeatureFilter.And(new FeatureFilter.OverlapsLocation(context.getRange()),
                                  spf), false);
        }
        else
        {
            fh = context.getFeatures().filter(new
                FeatureFilter.And(new FeatureFilter.OverlapsLocation(context.getSecondaryRange()),
                                  spf), false);
        }

        for (Iterator fi = fh.features(); fi.hasNext();)
        {
            SimilarityPairFeature f1 = (SimilarityPairFeature) fi.next();
            SimilarityPairFeature f2 = f1.getSibling();

            Strand s1 = f1.getStrand();
            Strand s2 = f2.getStrand();

            Location loc1 = f1.getLocation();
            Location loc2 = f2.getLocation();

            int min1 = loc1.getMin();
            int max1 = loc1.getMax();

            int min2 = loc2.getMin();
            int max2 = loc2.getMax();

            if (context.getDirection() == PairwiseRenderContext.HORIZONTAL)
            {
                float posX1 = (float) context.sequenceToGraphics(min1);
                float posY1 = (float) context.secondarySequenceToGraphics(min2);
                float posX2 = (float) context.sequenceToGraphics(max1);
                float posY2 = (float) context.secondarySequenceToGraphics(max2);

                if (s1 == s2)
                    line.setLine(posX1, posY1, posX2, posY2);
                else
                    line.setLine(posX2, posY1, posX1, posY2);
            }
            else
            {
                float posY1 = (float) context.sequenceToGraphics(min1);
                float posX1 = (float) context.secondarySequenceToGraphics(min2);
                float posY2 = (float) context.sequenceToGraphics(max1);
                float posX2 = (float) context.secondarySequenceToGraphics(max2);

                if (s1 == s2)
                    line.setLine(posX1, posY1, posX2, posY2);
                else
                    line.setLine(posX2, posY1, posX1, posY2);
            }

            Rectangle2D clip = g2.getClipBounds();

            clipLine((float) clip.getMinX(), (float) clip.getMaxX(),
                     (float) clip.getMinY(), (float) clip.getMaxY(),line);

            g2.setPaint(outline);
            g2.draw(line);
        }
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
     * <code>processMouseEvent</code> acts on a mouse gesture. The
     * target object is a <code>FeatureHolder</code> containing the
     * features on the primary sequence which contain the mouse
     * pointer.
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

        double gPos;

        if (context.getDirection() == PairwiseRenderContext.HORIZONTAL)
            gPos = me.getPoint().getX();
        else
            gPos = me.getPoint().getY();

        int priMin = context.graphicsToSequence(gPos);
        int priMax = context.graphicsToSequence(gPos + 1);

        FeatureHolder fh = context.getFeatures().filter(new
            FeatureFilter.And(new FeatureFilter.OverlapsLocation(new
                RangeLocation(priMin, priMax)), spf), false);

        return new SequenceViewerEvent(this, fh, priMin, me, path);
    }

    /**
     * <code>clipLine</code> clips the line to within the rectangle.
     * Cohen-Sutherland clipping implementation by Leen Ammeraal.
     *
     * @param xMin a <code>float</code>.
     * @param xMax a <code>float</code>.
     * @param yMin a <code>float</code>.
     * @param yMax a <code>float</code>.
     * @param line a <code>Line2D.Float</code>.
     */
    private void clipLine(float xMin, float xMax,
                          float yMin, float yMax, Line2D.Float line)
    {
        int clipTypeX = calcClipType(xMin, xMax,
                                     yMin, yMax,
                                     line.x1, line.y1);

        int clipTypeY = calcClipType(xMin, xMax,
                                     yMin, yMax,
                                     line.x2, line.y2);

        float dx, dy;

        if ((clipTypeX | clipTypeY) != 0)
        {
            if ((clipTypeX & clipTypeY) != 0)
                return;

            dx = line.x2 - line.x1;
            dy = line.y2 - line.y1;

            if (clipTypeX != 0)
            {
                if ((clipTypeX & 8) == 8)
                {
                    line.y1 += (xMin - line.x1) * dy / dx;
                    line.x1 = xMin;
                }
                else if ((clipTypeX & 4) == 4)
                {
                    line.y1 += (xMax - line.x1) * dy / dx;
                    line.x1 = xMax;
                }
                else if ((clipTypeX & 2) == 2)
                {
                    line.x1 += (yMin - line.y1) * dx / dy;
                    line.y1 = yMin;
                }
                else if ((clipTypeX & 1) == 1)
                {
                    line.x1 += (yMax - line.y1) * dx / dy;
                    line.y1 = yMax;
                }
            }
            else if (clipTypeY != 0)
            {
                if ((clipTypeY & 8) == 8)
                {
                    line.y2 += (xMin - line.x2) * dy / dx;
                    line.x2 = xMin;
                }
                else if ((clipTypeY & 4) == 4)
                {
                    line.y2 += (xMax - line.x2) * dy / dx;
                    line.x2 = xMax;
                }
                else if ((clipTypeY & 2) == 2)
                {
                    line.x2 += (yMin - line.y2) * dx / dy;
                    line.y2 = yMin;
                }
                else if ((clipTypeY & 1) == 1)
                {
                    line.x2 += (yMax - line.y2) * dx / dy;
                    line.y2 = yMax;
                }
            }
        }
    }

    /**
     * <code>calcClipType</code> calculates which sort of clipping to
     * carry out and returns the relevant code.
     *
     * @param xMin a <code>float</code>.
     * @param xMax a <code>float</code>.
     * @param yMin a <code>float</code>.
     * @param yMax a <code>float</code>.
     * @param x a <code>float</code>.
     * @param y a <code>float</code>.
     *
     * @return an <code>int</code> code.
     */
    private int calcClipType(float xMin, float xMax,
                             float yMin, float yMax, float x, float y)
    {
        return
            (x < xMin ? 8 : 0) |
            (x > xMax ? 4 : 0) |
            (y < yMin ? 2 : 0) |
            (y > yMax ? 1 : 0);
    }
}
