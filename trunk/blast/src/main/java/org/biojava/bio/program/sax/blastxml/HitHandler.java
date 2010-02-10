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
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.ContentHandler;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;

/**
 * @author David Huen
 */
class HitHandler
    extends StAXFeatureHandler
{
    // create static factory class that makes an instance
    // of this class.
    public final static StAXHandlerFactory HIT_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new HitHandler(staxenv);
            }
        };

    // class variables
    AttributesImpl hitAttrs = null;
    AttributesImpl hitIdAttrs = null;
    String sHit_def = null;

    // constructor
    public HitHandler(StAXFeatureHandler staxenv)
    {
        super(staxenv);
//        System.out.println("HitHandler staxenv " + staxenv);
//        // acquire value of <Hit_num> with inner class.
//        super.addHandler(new ElementRecognizer.ByLocalName("Hit_num"),
//            HitPropertyHandler.HIT_PROPERTY_HANDLER_FACTORY);

        // acquire value of <Hit_id> with inner class.
        super.addHandler(new ElementRecognizer.ByLocalName("Hit_id"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {

                        public void startElement(
                            String nsURI,
                            String localName,
                            String qName,
                            Attributes attrs,
                            DelegationManager dm)
                            throws SAXException
                        {
                            // generate start of <biojava:HitDescription>
                            hitIdAttrs = new AttributesImpl();

                            // must call superclass to keep track of levels
                            super.startElement(nsURI, localName, qName, attrs, dm);
                        }

                        public void setStringValue(String s) {
                            hitIdAttrs.addAttribute(biojavaUri, "id", "id", CDATA, s.trim());
                        }
                    };
                }
            }
        );

        // acquire value of <Hit_def> with inner class.
        super.addHandler(new ElementRecognizer.ByLocalName("Hit_def"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {

                    return new StringElementHandlerBase() {

                        public void setStringValue(String s)  throws SAXException {
                            sHit_def = s.trim();
                        }
                    };
                }
            }
        );

        // acquire value of <Hit_len> with inner class.  This needs to become and attribute of Hit.
        super.addHandler(new ElementRecognizer.ByLocalName("Hit_len"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void startElement(
                            String nsURI,
                            String localName,
                            String qName,
                            Attributes attrs,
                            DelegationManager dm)
                            throws SAXException
                        {
                            if (hitAttrs == null) hitAttrs = new AttributesImpl();

                            // must call superclass to keep track of levels
                            super.startElement(nsURI, localName, qName, attrs, dm);                            
                        }

                        public void setStringValue(String s)  throws SAXException {
                            hitAttrs.addAttribute(biojavaUri, "sequenceLength", "sequenceLength", CDATA, s.trim());
                        }

                        public void endElement(
                            String nsURI,
                            String localName,
                            String qName,
                            StAXContentHandler handler)
                            throws SAXException
                        {
                            // necessary as staxenv cannot be final and therefore
                            // staxenv.listener cannot be accessed from inner class
                            ContentHandler listener = getListener();

                            // get superclass to process the PCDATA for this element
                            super.endElement(nsURI, localName, qName, handler);

                            // we now generate the Hit element
                            listener.startElement(biojavaUri, "Hit", biojavaUri + ":Hit", hitAttrs);

                            // create <biojava:HitId> element
                            if (hitIdAttrs != null) {
                                hitIdAttrs.addAttribute(biojavaUri, "metaData", "metaData", CDATA, "none");
                                listener.startElement(biojavaUri, "HitId", biojavaUri + ":HitId", hitIdAttrs);
                                listener.endElement(biojavaUri, "HitId", biojavaUri + ":HitId");
                            }

                            // generate start of <biojava:HitDescription>
                            if (sHit_def != null) {
                                listener.startElement(biojavaUri, "HitDescription", biojavaUri + ":HitDescription", new AttributesImpl());
                                listener.characters(sHit_def.toCharArray(), 0, sHit_def.length());
                                listener.endElement(biojavaUri, "HitDescription", biojavaUri + ":HitDescription");
                            }
                        }
                    };
                }
            }
        );

        // handle <Hit_hsps> with its own handler.
        // the handling here is a tad perverse in that the sequence length has 
        // to be saved as an attribute of <Hit> although it is present as
        // an element.  This would mean that it cannot be created in the
        // startElementHandler.  Creating it in the endElementHandler would
        // cause it to fail to contain the child elements correctly.
        super.addHandler(new ElementRecognizer.ByLocalName("Hit_hsps"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new HitHspsHandler(staxenv) {
                        public void startElementHandler(
                            String nsURI,
                            String localName,
                            String qName,
                            Attributes attrs)
                            throws SAXException
                        {
                            // now I generate my own start element
                            super.startElementHandler(nsURI, localName, qName, attrs);
                        }
                    };
                }            
            }
        );
    }

    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler handler)
             throws SAXException
    {
        staxenv.listener.endElement(biojavaUri, "Hit", biojavaUri + ":Hit");
    }    
}
