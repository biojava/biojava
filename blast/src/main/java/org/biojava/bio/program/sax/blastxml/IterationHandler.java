/*
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
import org.xml.sax.SAXException;

/**
 * @author David Huen
 */
class IterationHandler
    extends StAXFeatureHandler
{
    // create static factory class that makes an instance
    // of this class.
    public final static StAXHandlerFactory ITERATION_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new IterationHandler(staxenv);
            }
        };

    // constructor
    public IterationHandler(StAXFeatureHandler staxenv)
    {
        super(staxenv);
//        System.out.println("IterationHandler staxenv " + staxenv);
        // delegate handling of <Iteration_hits>
        super.addHandler(new ElementRecognizer.ByLocalName("Iteration_hits"),
            IterationHitsHandler.ITERATION_HITS_HANDLER_FACTORY);

//        // handle <Iteration_iter-num> internally.
//        super.addHandler(new ElementRecognizer.ByLocalName("Iteration_iter-num"),
//            HitPropertyHandler.HIT_PROPERTY_HANDLER_FACTORY);
    }

    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler handler)
             throws SAXException
    {
    }

}
