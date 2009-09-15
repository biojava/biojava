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
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.geom.AffineTransform;
import java.io.Serializable;
import java.net.URL;

import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.OptimizableFilter;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.net.URLFactory;

/**
 * <code>RectangularImapRenderer</code> is a decorator for
 * <code>RectangularBeadRenderer</code> which adds the ability to
 * create HTML image map coordinates which correspond to the feature
 * rendering produced by the <code>RectangularBeadRenderer</code>.
 *
 * @author Keith James
 * @since 1.3
 */
public class RectangularImapRenderer
    implements BeadFeatureRenderer, ImageMapRenderer, Serializable
{
    private RectangularBeadRenderer renderer;
    private ImageMap imageMap;
    private URLFactory urlFactory;

    /**
     * Creates a new <code>RectangularImapRenderer</code>.
     *
     * @param renderer a <code>RectangularBeadRenderer</code>.
     * @param imageMap an <code>ImageMap</code>.
     * @param urlFactory a <code>URLFactory</code> which should be
     * capable of creating a suitable URL from each
     * <code>Feature</code> on the <code>Sequence</code> to be
     * rendered.
     */
    public RectangularImapRenderer(RectangularBeadRenderer renderer,
                                   ImageMap                imageMap,
                                   URLFactory              urlFactory)
    {
        this.renderer   = renderer;
        this.imageMap   = imageMap;
        this.urlFactory = urlFactory;
    }

    /**
     * <code>getImageMap</code> returns the current image map.
     *
     * @return an <code>ImageMap</code>.
     */
    public ImageMap getImageMap()
    {
        return imageMap;
    }

    /**
     * <code>setImageMap</code> sets the current image map.
     *
     * @param imageMap an <code>ImageMap</code>.
     */
    public void setImageMap(ImageMap imageMap)
    {
        this.imageMap = imageMap;
    }

    /**
     * <code>setDelegateRenderer</code> for the specified filter.
     *
     * @param filter an <code>OptimizableFilter</code>.
     * @param renderer a <code>BeadFeatureRenderer</code>.
     */
    public void setDelegateRenderer(OptimizableFilter   filter,
                                    BeadFeatureRenderer renderer)
    {
        this.renderer.setDelegateRenderer(filter, renderer);
    }

    /**
     * <p><code>renderImageMap</code> writes a set of image map
     * coordinates corresponding to the rectangle drawn by the
     * renderer. The hotspots created by this method have the rendered
     * <code>Feature</code> set as their user object.</p>
     *
     * <p>This method is called by <code>renderFeature</code> when a
     * raster image is rendered.</p>
     *
     * @param g2 a <code>Graphics2D</code>.
     * @param f a <code>Feature</code>.
     * @param context a <code>SequenceRenderContext</code>.
     */
    public void renderImageMap(Graphics2D                 g2,
                               Feature                     f,
                               SequenceRenderContext context)
    {
        Rectangle bounds = g2.getDeviceConfiguration().getBounds();

        // Safe to cast as bounds come from raster
        int  mapWidth = (int) bounds.getWidth();
        int mapHeight = (int) bounds.getHeight();

        URL url = urlFactory.createURL(f);

        double        beadDepth = getBeadDepth();
        double beadDisplacement = getBeadDisplacement();
        boolean scaleHeight     = getHeightScaling();

        AffineTransform t = g2.getTransform();
        double transX = t.getTranslateX();
        double transY = t.getTranslateY();

        int min, max, dif, x1, y1, x2, y2;
        double posXW, posYN, width, height;

        Location loc = f.getLocation();
        
        min = loc.getMin();
        max = loc.getMax();
        dif = max - min;

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
        }

        // Apply translation and round
        x1 = (int) Math.floor(posXW + transX);
        y1 = (int) Math.floor(posYN + transY);
        x2 = (int) Math.floor(posXW + width + transX);
        y2 = (int) Math.floor(posYN + height + transY);

        // If the whole rectangle is outside the image then ignore
        // it
        if (! (x1 > mapWidth || y1 > mapHeight || x2 < 0 || y2 < 0))
        {
            x1 = Math.max(x1, 0);
            y1 = Math.max(y1, 0);
            x2 = Math.min(x2, mapWidth);
            y2 = Math.min(y2, mapHeight);

            Integer [] coordinates = new Integer[4];
            coordinates[0] = new Integer(x1);
            coordinates[1] = new Integer(y1);
            coordinates[2] = new Integer(x2);
            coordinates[3] = new Integer(y2);

            imageMap.addHotSpot(new ImageMap.HotSpot(ImageMap.RECT,
                                                     url, coordinates, f));
        }
    }

    public void renderFeature(Graphics2D                 g2,
                              Feature                     f,
                              SequenceRenderContext context)
    {
        renderImageMap(g2, f, context);
        renderer.renderFeature(g2, f, context);
    }

    public void renderBead(Graphics2D                 g2,
                           Feature                     f,
                           SequenceRenderContext context)
    {
        renderer.renderBead(g2, f, context);
    }

    public double getDepth(SequenceRenderContext context)
    {
        return renderer.getDepth(context);
    }

    public double getBeadDepth()
    {
        return renderer.getBeadDepth();
    }

    public double getBeadDisplacement()
    {
        return renderer.getBeadDisplacement();
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
        return renderer.getHeightScaling();
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
        renderer.setHeightScaling(isEnabled);
    }

    public FeatureHolder processMouseEvent(FeatureHolder         holder,
                                           SequenceRenderContext context,
                                           MouseEvent            mEvent)
    {
        return renderer.processMouseEvent(holder, context, mEvent);
    }
}
