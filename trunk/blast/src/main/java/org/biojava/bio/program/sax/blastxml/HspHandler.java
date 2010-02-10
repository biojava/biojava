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
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;

/**
 * This class handles the <Hsp> element in NCBI Blast XML.
 * It generates events of the type handled by SearchContentHandler.
 * Most events will be generated just using the element local name
 * as the key and the CDATA as the value.
 *
 * The events do not conform to the BlastLikeDataSetCollection DTD
 * conformity is achieved with an adaptor class that intercepts
 * certain event types and translates them to that used in the
 * above DTD.  This was done to focus all the changes necessary
 * to achieve conformity in the adaptor class.
 */
class HspHandler
    extends StAXFeatureHandler
{
    // create static factory class that makes an instance
    // of this class.
    public final static StAXHandlerFactory HSP_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new HspHandler(staxenv);
            }
        };

    // local constants
    private static final String bitScore = "bitScore";
    private static final String score = "score";
    private static final String expectValue = "expectValue";
    private static final String numberOfIdentities = "numberOfIdentities";
    private static final String numberOfPositives = "numberOfPositives";
    private static final String alignmentSize = "alignmentSize";
    private static final String queryFrame = "queryFrame";
    private static final String hitFrame = "hitFrame";
    private static final String queryStrand = "queryStrand";
    private static final String hitStrand = "hitStrand";
    private static final String percentageIdentity = "percentageIdentity";

    // class variables
    AttributesImpl hspAttrs;
    AttributesImpl alignAttrs;

    // variables for temp storage
    int iNumberOfIdentities = Integer.MIN_VALUE;
    int iAlignmentSize = Integer.MIN_VALUE;
    String sHsp_qseq = null;
    String sHsp_hseq = null;
    String sHsp_midline = null;
    String sHsp_hit_from = null;
    String sHsp_hit_to = null;
    String sHsp_query_from = null;
    String sHsp_query_to = null;


    // constructor
    public HspHandler(StAXFeatureHandler staxenv)
    {
        super(staxenv);
//        System.out.println("HspHandler staxenv " + staxenv);
        // delegate handling of <Hsp_num>
//        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_num"),
//            SubHitPropertyHandler.SUBHIT_PROPERTY_HANDLER_FACTORY);

        // delegate handling of <Hsp_bit-score>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_bit-score"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) {
                            hspAttrs.addAttribute(biojavaUri, bitScore, bitScore, CDATA, s);
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_score>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_score"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) {
                            hspAttrs.addAttribute(biojavaUri, score, score, CDATA, s);
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_evalue>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_evalue"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) {
                            hspAttrs.addAttribute(biojavaUri, expectValue, expectValue, CDATA, s);
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_query-from>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_query-from"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) {
                            sHsp_query_from = s.trim();
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_query-to>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_query-to"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) {
                            sHsp_query_to = s.trim();
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_hit-from>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_hit-from"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) {
                            sHsp_hit_from = s.trim();
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_hit-to>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_hit-to"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) {
                            sHsp_hit_to = s.trim();
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_query-frame>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_query-frame"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {

                        public void setStringValue(String s) throws SAXException {
                            // save this to compute the percentage identity later.
                            int frameNo = Integer.parseInt(s.trim());

                            // convert the frame to the required format and return it
                            if (hitSequenceType.equals("protein")) {
                                hspAttrs.addAttribute(biojavaUri, queryFrame, queryFrame, CDATA, 
                                    stringifyFrame(frameNo));
                            }
                            else if (hitSequenceType.equals("dna")) {
                                // for some peculiar reason, when Hsp_hit-frame is reversed, it is
                                // the query frame sequence that is depicted inverted!!
                                // I assume it works the other way too although that never happens.
                                hspAttrs.addAttribute(biojavaUri, hitStrand, hitStrand, CDATA,
                                    stringifyStrand(frameNo));
                            }
                            else throw new SAXException("illegal sequence type");
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_hit-frame>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_hit-frame"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {

                        public void setStringValue(String s) throws SAXException {
                            // save this to compute the percentage identity later.
                            int frameNo = Integer.parseInt(s.trim());

                            // convert the frame to the required format and return it
                            if (hitSequenceType.equals("protein")) {
                                hspAttrs.addAttribute(biojavaUri, hitFrame, hitFrame, CDATA, 
                                    stringifyFrame(frameNo));
                            }
                            else if (hitSequenceType.equals("dna")) {
                                // for some peculiar reason, when Hsp_hit-frame is reversed, it is
                                // the query frame sequence that is depicted inverted!!
                                hspAttrs.addAttribute(biojavaUri, queryStrand, queryStrand, CDATA,
                                    stringifyStrand(frameNo));
                            }
                            else throw new SAXException("illegal sequence type");
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_identity>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_identity"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) {
                            iNumberOfIdentities = Integer.parseInt(s.trim());
                            hspAttrs.addAttribute(biojavaUri, numberOfIdentities, numberOfIdentities, CDATA, s);
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_positive>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_positive"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) {
                            hspAttrs.addAttribute(biojavaUri, numberOfPositives, numberOfPositives, CDATA, s);
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_align-len>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_align-len"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) {
                            iAlignmentSize = Integer.parseInt(s.trim());
                            hspAttrs.addAttribute(biojavaUri, alignmentSize, alignmentSize, CDATA, s);
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_qseq>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_qseq"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) {
                            sHsp_qseq = s;
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_hseq>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_hseq"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) {
                            sHsp_hseq = s;
                        }
                    };
                }
            }
        );

        // delegate handling of <Hsp_midline>
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp_midline"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) {
                            sHsp_midline = s;
                        }
                    };
                }
            }
        );


    }

    private String stringifyFrame(int frame) throws SAXException {
        switch (frame) {
            case -3: return "minus3";
            case -2: return "minus2";
            case -1: return "minus1";
            case 1: return "plus1";
            case 2: return "plus2";
            case 3: return "plus3";           
            default: throw new SAXException("illegal frame number encountered. ("+frame+")");
        }
    }

    private String stringifyStrand(int strand) throws SAXException {
        if (strand > 0) return "plus";
        else if (strand < 0) return "minus";
        else throw new SAXException("illegal strand number encountered.");
    }

    public void startElementHandler(
            String nsURI,
            String localName,
            String qName,
            Attributes attrs)
             throws SAXException 
    { 
        // create an AttributesImpl to save the attributes to
        hspAttrs = new AttributesImpl();
    }


    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler handler)
             throws SAXException
    {
        // only generate the element if key parts are present
        if ((iNumberOfIdentities == Integer.MIN_VALUE)
            || (iAlignmentSize == Integer.MIN_VALUE)
            || (sHsp_qseq == null)
            || (sHsp_hseq == null)
            || (sHsp_midline == null)
            || (sHsp_hit_from == null)
            || (sHsp_hit_to == null)
            || (sHsp_query_from == null)
            || (sHsp_query_to == null) ) {
            throw new SAXException("<Hsp> is non-compliant.");
        }

        // compute percentage identity and report it
        hspAttrs.addAttribute(biojavaUri, 
            percentageIdentity, 
            percentageIdentity, 
            CDATA, 
            Float.toString( ((float) (100 * iNumberOfIdentities)) / ((float) iAlignmentSize))
            );

        // generate start of <biojava:Hsp>
        staxenv.listener.startElement(biojavaUri, "HSP", biojavaUri + ":HSP", new AttributesImpl());

            // generate <biojava:HSPSummary>
            staxenv.listener.startElement(biojavaUri, "HSPSummary", biojavaUri + ":" + "HSPSummary", hspAttrs);
            staxenv.listener.endElement(biojavaUri, "HSPSummary", biojavaUri + ":" + "HSPSummary");

            // generate the <biojava:BlastLikeAlignment>
            staxenv.listener.startElement(biojavaUri, "BlastLikeAlignment", biojavaUri + ":BlastLikeAlignment", new AttributesImpl());

                // generate start of <biojava:QuerySequence>
                AttributesImpl queryAttrs = new AttributesImpl();
                queryAttrs.addAttribute(biojavaUri, "startPosition", "startPosition", CDATA, sHsp_query_from);
                queryAttrs.addAttribute(biojavaUri, "stopPosition", "stopPosition", CDATA, sHsp_query_to);
                staxenv.listener.startElement(biojavaUri, "QuerySequence", biojavaUri + ":QuerySequence", queryAttrs);

                // pass the sequence symbol tokens over
                staxenv.listener.characters(sHsp_qseq.toCharArray(), 0, sHsp_qseq.length());

                // generate end of <biojava:QuerySequence>
                staxenv.listener.endElement(biojavaUri, "QuerySequence", biojavaUri + ":QuerySequence");

                // generate start of <biojava:MatchConsensus>
                AttributesImpl matchAttrs = new AttributesImpl();
                matchAttrs.addAttribute("xml", "space", "xml:space", CDATA, "preserve");
                staxenv.listener.startElement(biojavaUri, "MatchConsensus", biojavaUri + ":MatchConsensus", matchAttrs);

                // pass the sequence symbol tokens over
                staxenv.listener.characters(sHsp_midline.toCharArray(), 0, sHsp_midline.length());

                // generate end of <biojava:MatchConsensus>
                staxenv.listener.endElement(biojavaUri, "MatchConsensus", biojavaUri + ":MatchConsensus");

                // generate start of <biojava:HitSequence>
                AttributesImpl hitAttrs = new AttributesImpl();
                hitAttrs.addAttribute(biojavaUri, "startPosition", "startPosition", CDATA, sHsp_hit_from);
                hitAttrs.addAttribute(biojavaUri, "stopPosition", "stopPosition", CDATA, sHsp_hit_to);
                staxenv.listener.startElement(biojavaUri, "HitSequence", "HitSequence", hitAttrs);

                // pass the sequence symbol tokens over
                staxenv.listener.characters(sHsp_hseq.toCharArray(), 0, sHsp_hseq.length());

                // generate end of <biojava:HitSequence>
                staxenv.listener.endElement(biojavaUri, "HitSequence", biojavaUri + ":HitSequence");

            // generate end of <biojava::BlastLikeAlignment>    
            staxenv.listener.endElement(biojavaUri, "BlastLikeAlignment", biojavaUri + ":BlastLikeAlignment"); 

        // generate end of <biojava:Hsp>
        staxenv.listener.endElement(biojavaUri, "HSP", biojavaUri + ":HSP");
    }
}
