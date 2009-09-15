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

import org.biojava.bio.seq.Feature;

/**
 * <code>ImageMapRenderer</code>s create strings representing
 * <code>Feature</code>s suitable for use in HTML image
 * maps. Typically an <code>ImageMapRenderer</code> will be used as a
 * decorator on a <code>FeatureRenderer</code> which will draw the
 * corresponding image area(s).
 *
 * @author Keith James
 * @since 1.3
 */
public interface ImageMapRenderer extends FeatureRenderer
{
    /**
     * <code>renderImageMap</code> renders the <code>Feature</code> as
     * set of image map hotspots.
     *
     * @param g2 a <code>Graphics2D</code>.
     * @param f a <code>Feature</code>.
     * @param context a <code>SequenceRenderContext</code>.
     */
    public void renderImageMap(Graphics2D                 g2,
                               Feature                     f,
                               SequenceRenderContext context);
}
