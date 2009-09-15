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

package org.biojava.bio.seq.io.game12;

import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.game.ElementRecognizer;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 *  Handles the GAME &lt;feature_set&gt; element
 *
 * this element is used to represent transcripts.
 *
 * @author     David Huen
 * @since      1.2
 */
public class GAMEFeatureSetHandler
         extends StAXFeatureHandler {
    // this creates a bundle of data on a specific
    // aspect of an annotation, e.g. a transcript.

    // database columns
//    private String id;
//    private String produces_seq;
    private String featureSetType;

    StrandedFeature.Template template;

    Location transcript;
    Location translateOffset;
    StrandedFeature.Strand transcriptStrand;

    // subclass GAMEFeatureSpanHandler to retrieve data
    private class FeatureSpanHandler extends GAMEFeatureSpanHandler
    {
        private FeatureSpanHandler(StAXFeatureHandler staxenv)
        {
            super(staxenv);
//            System.out.println("entering FeatureSpanHandler");
        }

        public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler contentHandler)
            throws SAXException
        {
            // validate
            super.endElementHandler(nsURI, localName, qName, contentHandler);

            // decide what to do on basis of type
//            System.out.println("SpanHandler end seeing " + featureSpanType);
            if (featureSpanType.equals("exon")) {
                // got an exon
                if (transcript == null) {
                    // store transcript data
                    transcript = featureSpanLoc;
                    transcriptStrand = featureSpanStrand;
                }
                else {
                    // validate this exon, is it on same strand?
                    if (transcriptStrand != featureSpanStrand) {
                        System.err.println("exons on differing strands!");
                    }

                    // update the transcript
//                    System.out.println("setting exon in transcript " + featureSpanLoc);
                    transcript = LocationTools.union(transcript, featureSpanLoc);
                }
            }
            else if (featureSpanType.equals("translate offset")) {
                // stash the location away
                if (translateOffset != null)
                    System.err.println("translate offset multiply defined.");
                else {
                    if (transcriptStrand == null) 
                        transcriptStrand = featureSpanStrand;
                    else if (transcriptStrand != featureSpanStrand)
                        System.err.println("exons on differing strands!");

                    translateOffset = featureSpanLoc;
                }
            }
        }
    }

    // set up factory method
    /**
     *  Description of the Field
     */
    public final static StAXHandlerFactory GAME_FEATURE_SET_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new GAMEFeatureSetHandler(staxenv);
            }
        };


    GAMEFeatureSetHandler(StAXFeatureHandler staxenv) {
        // setup environment
        super(staxenv);

        // setup handlers
        // <name>
        super.addHandler(new ElementRecognizer.ByLocalName("name"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new NameHandler();
                }
            }
                );
        // <type> local handler type
        super.addHandler(new ElementRecognizer.ByLocalName("type"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new TypeHandler();
                }
            }
                );
        // <author> local handler type
        super.addHandler(new ElementRecognizer.ByLocalName("author"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new AuthorHandler();
                }
            }
                );

        // <feature_span>
//        super.addHandler(new ElementRecognizer.ByLocalName("feature_span"),
//                GAMEFeatureSpanHandler.GAME_FEATURE_SPAN_HANDLER_FACTORY);
        super.addHandler(new ElementRecognizer.ByLocalName("feature_span"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new FeatureSpanHandler(staxenv);
                }
            }
        );
        // <seq>
//        super.addHandler(new ElementRecognizer.ByLocalName("seq"),
//                GAMESeqHandler.GAME_SEQ_HANDLER_FACTORY);
        // <property>
        super.addHandler(new ElementRecognizer.ByLocalName("property"),
                GAMEPropertyHandler.GAME_PROPERTY_HANDLER_FACTORY);
    }

    private class NameHandler extends StringElementHandlerBase {
        /**
         *  Sets the stringValue attribute of the NameHandler object
         *
         *@param  s  The new stringValue value
         */
        protected void setStringValue(String s) {
        }
    }

    private class TypeHandler extends StringElementHandlerBase {
        /**
         *  Sets the stringValue attribute of the TypeHandler object
         *
         *@param  s  The new stringValue value
         */
        protected void setStringValue(String s) {
            featureSetType = s.trim();
        }
    }

    private class AuthorHandler extends StringElementHandlerBase {
        /**
         *  Sets the stringValue attribute of the TypeHandler object
         *
         *@param  s  The new stringValue value
         */
        protected void setStringValue(String s) {
        }
    }

    public void startElementHandler(
            String nsURI,
            String localName,
            String qName,
            Attributes attrs) 
        throws SAXException
    {
        // acquire attributes here
        String id = attrs.getValue("id");
        String produces_seq = attrs.getValue("produces_seq");

        // create template
        // this element is expected to encode a transcript
        template = new StrandedFeature.Template();
        template.annotation = new SimpleAnnotation();

        // fill in annotation bundle
        try {
            try {
                if (id != null) template.annotation.setProperty("id", id);
                if (produces_seq != null) template.annotation.setProperty("produces_seq", produces_seq);
            }
            catch (ChangeVetoException cve) {
                cve.printStackTrace();
                throw new SAXException("unexpected ChangeVetoException when setting id!");
            }

            listener.startFeature(template);
        }
        catch (ParseException pe) {
            pe.printStackTrace();
            throw new SAXException("error in GAMEFeatureSetHandler.");
        }
    }


    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler contentHandler) 
        throws SAXException
    {
        // issue warning if it is not a transcript
        if (!featureSetType.equals("transcript")) {
             System.err.println("<feature_set> of type " + featureSetType + " encountered when transcript expected");
        }

        // fill in template as best can
        template.type = featureSetType;
        template.source = "";
        template.location = transcript;
        template.strand = transcriptStrand;

        // indicate end of sequence
        try {

            // set the translate offset if available
            try {
                if (translateOffset != null) template.annotation.setProperty("translate_offset", translateOffset);
            }
            catch (ChangeVetoException cve) {
                cve.printStackTrace();
                throw new SAXException("unexpected ChangeVetoException when setting id!");
            }

            listener.endFeature();
        }
        catch (ParseException pe) {
            pe.printStackTrace();
            throw new SAXException("error in GAMEFeatureSetHandler.");
        }
    }
}

