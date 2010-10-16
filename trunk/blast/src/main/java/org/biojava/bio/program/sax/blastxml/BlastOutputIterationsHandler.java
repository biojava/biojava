/**
 *  BioJava development code This code may be freely distributed and modified
 *  under the terms of the GNU Lesser General Public Licence. This should be
 *  distributed with the code. If you do not have a copy, see:
 *  http://www.gnu.org/copyleft/lesser.html Copyright for this code is held
 *  jointly by the individual authors. These should be listed in
 *
 *@author    doc comments. For more information on the BioJava project and its
 *      aims, or to join the biojava-l mailing list, visit the home page at:
 *      http://www.biojava.org/
 */
package org.biojava.bio.program.sax.blastxml;

import org.biojava.bio.seq.io.game.ElementRecognizer;
import org.biojava.utils.stax.StAXContentHandler;

class BlastOutputIterationsHandler
    extends StAXFeatureHandler
{
    // create static factory class that makes an instance
    // of this class.
    public final static StAXHandlerFactory BLASTOUTPUT_ITERATIONS_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new BlastOutputIterationsHandler(staxenv);
            }
        };

    // constructor
    public BlastOutputIterationsHandler(StAXFeatureHandler staxenv)
    {
        super(staxenv);
//        System.out.println("BlastOutputIterationsHandler staxenv " + staxenv);

        // delegate handling of <Iteration>
        super.addHandler(new ElementRecognizer.ByLocalName("Iteration"),
            IterationHandler.ITERATION_HANDLER_FACTORY);
    }
}
