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

import java.util.HashSet;
import java.util.Set;

import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.game.ElementRecognizer;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 *  Handles the GAME &lt;annotation&gt; element
 *
 * @author     David Huen
 * @since      1.2
 */
public class GAMEAnnotationHandler
         extends StAXFeatureHandler {
    // <annotation> is a container for all features of a "gene".
    // the only important property of this container is its id
    // which I need to capture and supply nested classes.

    private Set knownTypes;

    // database columns
    String annotationName;
    String annotationType;

    // there can be multiple transcripts (aka feature_sets) in a
    // <annotation>.  We must get the full extent of all transcripts
    // to use as the limits.
    int annotationLocMin = Integer.MAX_VALUE;
    int annotationLocMax = Integer.MIN_VALUE;

    StrandedFeature.Template annotationTemplate;

    // subclass GAMEFeatureSetHandler to retrieve transcript info
    private class FeatureSetHandler extends GAMEFeatureSetHandler
    {
        private FeatureSetHandler(StAXFeatureHandler staxenv)
        {
            super(staxenv);
//            System.out.println("entering FeatureSetHandler");
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

            // fill in the template with location and strand info
            annotationTemplate.strand = transcriptStrand;

            annotationLocMin = Math.min(annotationLocMin, transcript.getMin());
            annotationLocMax = Math.max(annotationLocMax, transcript.getMax());

            if (annotationTemplate.strand == null) 
                annotationTemplate.strand = transcriptStrand;
            else if (annotationTemplate.strand != transcriptStrand) {
                // conflicting info
                System.err.println("inconsistent strand info from transcripts.");
            }
        }
    }

    // subclass <dbxref> to write a feature property here
    private class DbxrefHandler extends GAMEDbxrefHandler
    {
        private DbxrefHandler(StAXFeatureHandler staxenv)
        {
            super(staxenv);
        }

        public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler contentHandler)
            throws SAXException
        {
            // validate before going further
            super.endElementHandler(nsURI, localName, qName, contentHandler);

            try {
                listener.addFeatureProperty("dbxref", "dbxref:" + db_xref_db + "//" + db_xref_id);
            }
            catch (ParseException pe) {
                pe.printStackTrace();
                throw new SAXException("unexpected exception while add <dbxref> as a feature property.");
            }
        }
    }

    // set up factory method
    /**
     *  Description of the Field
     */
    public final static StAXHandlerFactory GAME_ANNOTATION_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new GAMEAnnotationHandler(staxenv);
            }
        };

    /**
     *  Constructor for the GAMEAnnotationHandler object
     *
     *@param  staxenv   Description of the Parameter
     *@param  parentID  Description of the Parameter
     */
    GAMEAnnotationHandler(StAXFeatureHandler staxenv) {
        // setup environment
        super(staxenv);

        // initialise known types
        knownTypesInitialiser();

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

        // <seq>: never seen it used yet.
//       super.addHandler(new ElementRecognizer.ByLocalName("seq"),
//         GAMESeqPropHandler.GAME_SEQ_PROP_HANDLER_FACTORY);
        // <gene>
        super.addHandler(new ElementRecognizer.ByLocalName("gene"),
                GAMEGeneHandler.GAME_GENE_HANDLER_FACTORY);
        // <feature_set>
        super.addHandler(new ElementRecognizer.ByLocalName("feature_set"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new FeatureSetHandler(staxenv);
                }
            }
        );
        // <dbxref>
        super.addHandler(new ElementRecognizer.ByLocalName("dbxref"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new DbxrefHandler(staxenv);
                }
            }
        );
        // <Aspect>
        super.addHandler(new ElementRecognizer.ByLocalName("aspect"),
                GAMEAspectHandler.GAME_ASPECT_HANDLER_FACTORY);
        // <property>
        super.addHandler(new ElementRecognizer.ByLocalName("property"),
                GAMEPropertyHandler.GAME_PROPERTY_HANDLER_FACTORY);
    }


    /**
     *  Description of the Class
     *
     *@author     david
     *@created    19 January 2002
     */
    private class NameHandler extends StringElementHandlerBase {
        /**
         *  Sets the stringValue attribute of the NameHandler object
         *
         *@param  s  The new stringValue value
         */
        protected void setStringValue(String s) {
            annotationName = s.trim();
        }
    }

    /**
     *  Description of the Class
     *
     *@author     david
     *@created    19 January 2002
     */
    private class TypeHandler extends StringElementHandlerBase {
        /**
         *  Sets the stringValue attribute of the TypeHandler object
         *
         *@param  s  The new stringValue value
         */
        protected void setStringValue(String s) {
            annotationType = s.trim();
        }
    }

    private void knownTypesInitialiser()
    {
        // initialise a String array
        String [] types = {"gene", "tRNA", "snRNA", "pseudogene", "transposon", "snoRNA", "misc. non-coding RNA", 
            "transposable_element", "miscellaneous curator's observation"};

        // now initialise the knownTypes Set
        knownTypes = new HashSet();

        for (int i=0; i < types.length; i++) {
            knownTypes.add(types[i]);
        }
    }

    public void startElementHandler(
            String nsURI,
            String localName,
            String qName,
            Attributes attrs) 
        throws SAXException 
    {
        // indicate start of sequence
        try {
//            System.out.println("local name is " + localName);
            annotationTemplate = new StrandedFeature.Template();
            annotationTemplate.annotation = new SimpleAnnotation();

            // there should be an id in the attributes, get it
            String id = attrs.getValue("id");

            try {
                if (id != null) annotationTemplate.annotation.setProperty("id", id);
            }
            catch (ChangeVetoException cve) {
                cve.printStackTrace();
                throw new SAXException("unexpected ChangeVetoException when setting id!");
            }

            listener.startFeature(annotationTemplate);
        }
        catch (ParseException pe) {
            pe.printStackTrace();
            throw new SAXException("error in GAMEAnnotationHandler.");
        }
    }


    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler contentHandler) 
        throws SAXException
    {
        // fill in the template as best can
        // confirm that it is something I expect
        if (!knownTypes.contains(annotationType)) {
             System.err.println("<annotation> of type " + annotationType + " encountered when gene expected");
        }

        // ganbatte!
        annotationTemplate.type = annotationType;
        annotationTemplate.source = "";
        annotationTemplate.location = new RangeLocation(annotationLocMin, annotationLocMax);

        // bear in mind we are completely dependent on the SeqIOListener
        // to set the location and strand!!!!

        // indicate end of sequence
        try {
            listener.endFeature();
        }
        catch (ParseException pe) {
            pe.printStackTrace();
            throw new SAXException("error in GAMEAnnotationHandler.");
        }
    }
}

