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
import org.biojava.bio.seq.OptimizableFilter;

/**
 * <p><code>BeadFeatureRenderer</code>s use a 'string of beads'
 * metaphor for displaying features.</p>
 *
 * <p>A concrete <code>BeadFeatureRenderer</code> may render a series
 * of features in more than one style by delegating to other
 * <code>BeadFeatureRenderer</code>s for the additional style(s). This
 * is achieved using the <code>setDelegateRenderer()</code> method
 * which associates an <code>OptimizableFilter</code> with another
 * <code>BeadFeatureRenderer</code>. Any feature accepted by the
 * filter is rendered with that renderer, while the remainder are
 * rendered by the current renderer.</p>
 *
 * @author Keith James
 * @since 1.2
 */
public interface BeadFeatureRenderer extends FeatureRenderer
{
    /**
     * <code>getBeadDepth</code> returns the depth of a single bead
     * produced by the renderer.
     *
     * @return a <code>double</code>.
     */
    public double getBeadDepth();

    /**
     * <code>getBeadDisplacement</code> returns the displacement of
     * beads from the centre line of the renderer. A positive value
     * indicates displacment downwards (for horizontal renderers) or
     * to the right (for vertical renderers).
     *
     * @return a <code>double</code>.
     */
    public double getBeadDisplacement();

    /**
     * <code>setDelegateRenderer</code> associates an
     * <code>OptimizableFilter</code> with a
     * <code>BeadFeatureRenderer</code>. Any feature accepted by the
     * filter will be passed to the associated renderer for
     * drawing. The <code>OptimizableFilter</code>s should be disjoint
     * with respect to each other (a feature may not be rendered more
     * than once).
     *
     * @param filter an <code>OptimizableFilter</code>.
     * @param renderer a <code>BeadFeatureRenderer</code>.
     */
    public void setDelegateRenderer(OptimizableFilter   filter,
				    BeadFeatureRenderer renderer);

    /**
     * <code>renderBead</code> should implement rendering for this
     * bead type only. The <code>renderFeature</code> method is
     * expected to handle the calls to delegate renderers.
     *
     * @param g2 a <code>Graphics2D</code>.
     * @param f a <code>Feature</code> to render.
     * @param context a <code>SequenceRenderContext</code> context.
     */
    public void renderBead(Graphics2D            g2,
                           Feature               f,
                           SequenceRenderContext context);
}
