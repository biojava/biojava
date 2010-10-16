 /**
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
 *
 * @author David Huen
 * @author Hanning Ni     Doubletwist Inc
 */
public class AGAVEAssemblyHandler
               extends StAXFeatureHandler implements SequenceHandler

{
  public static final StAXHandlerFactory AGAVE_ASSEMBLY_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEAssemblyHandler(staxenv);
    }
  };

   AGAVEAssemblyHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    super( staxenv );
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("assembly", true);

    // setup handlers
        //
       super.addHandler(new ElementRecognizer.ByLocalName("bio_sequence"),
         AGAVEBioSeqHandler.AGAVE_BIO_SEQ_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("fragment_order"),
         AGAVEFragmentOrderHandler.AGAVE_FRAGMENT_ORDER_HANDLER_FACTORY);

  }

}
