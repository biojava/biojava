/*
 *  BioJava development code
 *
 *  This code may be freely distributed and modified under the
 *  terms of the GNU Lesser General Public Licence.  This should
 *  be distributed with the code.  If you do not have a copy,
 *  see:
 *
 *  http://www.gnu.org/copyleft/lesser.html
 *
 *  Copyright for this code is held jointly by the individual
 *  authors.  These should be listed in @author doc comments.
 *
 *  For more information on the BioJava project and its aims,
 *  or to join the biojava-l mailing list, visit the home page
 *  at:
 *
 *  http://www.biojava.org/
 *
 */

package org.biojava.bio.seq.io.game12;

import org.biojava.utils.stax.StAXContentHandler;

/**
 *  Factory for StAX content handlers.
 *
 *@author     Thomas Down
 *@since      1.2
 */

public interface StAXHandlerFactory {
    /**
     *  Return a suitable StAX content handler.
     *
     * @param  staxenv   Description of the Parameter
     * @return           The handler value
     */

    public StAXContentHandler getHandler(StAXFeatureHandler staxenv);
}
