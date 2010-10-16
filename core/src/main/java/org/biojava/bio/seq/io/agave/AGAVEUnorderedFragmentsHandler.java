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
 * unordered_fragments
 *
 *
 * @author Hanning Ni    Doubletwist Inc
 * 
 */
public class AGAVEUnorderedFragmentsHandler
               extends StAXFeatureHandler

{
  public static final StAXHandlerFactory AGAVE_UNORDERED_FRAGMENTS_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEUnorderedFragmentsHandler(staxenv);
    }
  };


  AGAVEUnorderedFragmentsHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    super( staxenv );
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("unordered_fragment", true);

    // setup handlers
        //
       super.addHandler(new ElementRecognizer.ByLocalName("bio_sequence"),
         AGAVEBioSeqHandler.AGAVE_BIO_SEQ_HANDLER_FACTORY);

  }


}
