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
import java.awt.event.MouseEvent;
import java.util.List;

import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;

/**
 * <code>PairwiseSequenceRenderer</code>s render information about the
 * relationship between two sequences. Its function is analagous to
 * <code>SequenceRenderer</code> for single sequences and is
 * extensively based on that code.
 *
 * @author Keith James
 * @since 1.2
 */
public interface PairwiseSequenceRenderer
{
    /**
     * <code>paint</code>s some or all of the information about the
     * sequence pair.
     *
     * @param g2 a <code>Graphics2D</code>.
     * @param prc a <code>PairwiseRenderContext</code> encapsulating
     * the information to be displayed.
     */
    public void paint(Graphics2D g2, PairwiseRenderContext prc);

    /**
     * <code>processMouseEvent</code> produces a
     * <code>SequenceViewerEvent</code> in response to a mouse
     * gesture.
     *
     * @param prc a <code>PairwiseRenderContext</code>.
     * @param me a <code>MouseEvent</code> that caused the request.
     * @param path a <code>List</code> of
     * <code>PairwiseSequenceRenderer</code> instances passed through
     * so far.
     *
     * @return a <code>SequenceViewerEvent</code> encapsulating the
     * mouse gesture.
     */
    public SequenceViewerEvent processMouseEvent(PairwiseRenderContext prc,
                                                 MouseEvent            me,
                                                 List                  path);

    /**
     * <code>PairwiseRendererForwarder</code> forward events to other
     * renderers. This is closely based on the regular
     * <code>RendererForwarder</code>.
     */
    public static class PairwiseRendererForwarder extends ChangeForwarder
    {
        /**
         * Creates a new <code>PairwiseRendererForwarder</code>.
         *
         * @param source a <code>PairwiseSequenceRenderer</code>.
         * @param cs a <code>ChangeSupport</code>.
         */
        public PairwiseRendererForwarder(PairwiseSequenceRenderer source,
                                         ChangeSupport            cs)
        {
            super(source, cs);
        }
    
        /**
         * <code>generateEvent</code> generates events in response to
         * layout change and repaint requests.
         *
         * @param ce a <code>ChangeEvent</code>.
         *
         * @return a <code>ChangeEvent</code>. 
         */
        public ChangeEvent generateEvent(ChangeEvent ce)
        {
            ChangeType type = ce.getType();

            ChangeType newType;

            if (type.isMatchingType(SequenceRenderContext.LAYOUT))
                newType = SequenceRenderContext.LAYOUT;
            else if (type.isMatchingType(SequenceRenderContext.REPAINT))
                newType = SequenceRenderContext.REPAINT;
            else
                return null;

            return new ChangeEvent(getSource(), newType, null, null, ce);
        }
    }
}
