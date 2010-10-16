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

package org.biojava.bio.program.sax.blastxml;

import org.biojava.bio.seq.io.game.ElementRecognizer;
import org.biojava.utils.stax.StAXContentHandler;

/**
 * this class parses a preprocessed aggregation of &lt;BlastOutput&gt;
 * elements.
 * <p>
 * See BlastOputputHandler for more details.
 *
 * @author David Huen
 */
class BlastAggregator
    extends StAXFeatureHandler
{
    // constructor
    public BlastAggregator(StAXFeatureHandler staxenv)
    {
        // execute superclass
        super(staxenv);

        // delegate handling of <BlastOutput>
        super.addHandler(new ElementRecognizer.ByLocalName("BlastOutput"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new BlastOutputHandler(staxenv);
                }
            }
        );
    }
}
