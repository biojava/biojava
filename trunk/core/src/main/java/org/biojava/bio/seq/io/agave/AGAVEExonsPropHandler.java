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
 *
 * @author Hanning Ni    Doubletwist Inc
 */
public class AGAVEExonsPropHandler  extends StAXPropertyHandler{


   public static final StAXHandlerFactory AGAVE_EXONS_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEExonsPropHandler(staxenv);
    }
   };

   AGAVEExonsPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("exons", true);

    super.addHandler(new ElementRecognizer.ByLocalName("element_id"),
         AGAVEElementIdPropHandler.AGAVE_ELEMENT_ID_PROP_HANDLER_FACTORY);

  }


}
