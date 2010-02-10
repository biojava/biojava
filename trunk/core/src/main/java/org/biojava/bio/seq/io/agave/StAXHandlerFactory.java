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
package org.biojava.bio.seq.io.agave;

/**
 * Factory for StAX content handlers.
 * Roughly all xml elements are seperated into two parts:
 * those hold features and those hold annotations. (Hanningn)
 *
 * @author copied from Thomas Down
 */
public interface StAXHandlerFactory {
    /**
     * Return a suitable StAX content handler.
     *
     * @param staxenv An StAXFeatureHandler which defines the current
     *        XFF parsing `environment'.
     */

    public StAXContentHandler getHandler(StAXFeatureHandler staxenv);
}
