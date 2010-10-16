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
 * @author Hanning Ni     Doubletwist Inc
 */
public class AGAVEAnnotationsHandler
               extends StAXPropertyHandler

{
  public static final StAXHandlerFactory AGAVE_ANNOTATIONS_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEAnnotationsHandler(staxenv);
    }
  };

 AGAVEAnnotationsHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    super( staxenv );
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("annotations", true);

    // setup handlers
       //
       super.addHandler(new ElementRecognizer.ByLocalName("seq_feature"),
         AGAVESeqFeatureHandler.AGAVE_SEQ_FEATURE_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("gene"),
         AGAVEGeneHandler.AGAVE_GENE_HANDLER_FACTORY);

      super.addHandler(new ElementRecognizer.ByLocalName("comp_result"),
         AGAVECompResultHandler.AGAVE_COMP_RESULT_HANDLER_FACTORY);

  }
  /**
  protected Feature.Template createTemplate() {
    // create Gene Template for this
    StrandedFeature.Template st = new StrandedFeature.Template();

    // assume feature set to describe a transcript
    st.type = "annotations";
    st.strand = StrandedFeature.UNKNOWN;
    // set up annotation bundle
    st.annotation = new SmallAnnotation();
    st.location = new  Location.EmptyLocation();


    if( staxenv != null )
        staxenv. subFeatures .add( this ) ;

    return st;
  }**/



}

