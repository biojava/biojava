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
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 *  Handles the GAME &lt;feature_span&gt; element
 *
 *@author     David Huen
 *@since      1.2
 */
public class GAMEFeatureSpanHandler
         extends StAXFeatureHandler {
    // <annotation> is a container for all features of a "gene".
    // the only important property of this container is its id
    // which I need to capture and supply nested classes.

    // database columns
    private String featureSpanId;
    String featureSpanType;
    Location featureSpanLoc;
    StrandedFeature.Strand featureSpanStrand;

    private StrandedFeature.Template template;

//    private GAMESeqRelHandler.SeqRelTemplate seqRelTmpl = null;

    public class SeqRelHandler extends GAMESeqRelHandler
    {
        private SeqRelHandler(StAXFeatureHandler staxenv)
        {
            super(staxenv);
//            System.out.println("entering SeqRelHandler");
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

            featureSpanLoc = seqRelLoc;
            featureSpanStrand = seqRelStrand;
        }
    }

    // set up factory method
    /**
     *  Description of the Field
     */
    public final static StAXHandlerFactory GAME_FEATURE_SPAN_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new GAMEFeatureSpanHandler(staxenv);
            }
        };

    GAMEFeatureSpanHandler(StAXFeatureHandler staxenv) {
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

        // <type>
        super.addHandler(new ElementRecognizer.ByLocalName("type"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new TypeHandler();
                }
            }
                );
        // <seq_relationship> external handler type
        super.addHandler(new ElementRecognizer.ByLocalName("seq_relationship"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new SeqRelHandler(staxenv);
                }
            }
        );
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
            featureSpanType = s.trim();
        }
    }

    public void startElementHandler(
            String nsURI,
            String localName,
            String qName,
            Attributes attrs) 
        throws SAXException
    {
        featureSpanId = attrs.getValue("id");

        try {
            template = new StrandedFeature.Template();
            template.annotation = new SimpleAnnotation();
            template.source = null;

            try {
                if (featureSpanId != null) template.annotation.setProperty("id", featureSpanId);
            }
            catch (ChangeVetoException cve) {
                throw new SAXException("unexpected ChangeVetoException when setting id!");
            }

            listener.startFeature(template);
        }
        catch (ParseException pe) {
        }
    }

    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler contentHandler)
        throws SAXException 
    {
        // issue warning if it is not a translate offset or exon
        if ( !( (featureSpanType.equals("exon")) || (featureSpanType.equals("translate offset")) ) ) {
             System.err.println("<feature_span> of unexpected type " + featureSpanType);
        }

        // fill in template
        // i would like to avoid returning something for translate offsets
        // but I can't because the start feature has already allocated a template
        template.type = featureSpanType;
        template.source = "";
        template.location = featureSpanLoc;
        template.strand = featureSpanStrand;

        // indicate end of sequence
        try {
            listener.endFeature();
        }
        catch (ParseException pe) {
            pe.printStackTrace();
            throw new SAXException("error in GAMEFeatureSetHandler.");
        }
    }
}

