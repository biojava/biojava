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

package org.biojava.bio.program.ssbind;

import org.biojava.utils.stax.StAXContentHandler;

/**
 * <code>StAXHandlerFactory</code> is an interface for factories
 * producing <code>StAXContentHandler</code>s which are used by the
 * <code>SeqSimilarityStAXAdapter</code>.
 *
 * @author Thomas Down
 * @author Keith James
 * @since 1.3
 */
public interface StAXHandlerFactory
{
    /**
     * <code>getHandler</code> returns an appropriate
     * <code>StAXContentHandler</code> implementation containing a
     * reference to a parent context.
     *
     * @param ssContext a <code>SeqSimilarityStAXAdapter</code>.
     *
     * @return a <code>StAXContentHandler</code>.
     */
    public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext);
}
