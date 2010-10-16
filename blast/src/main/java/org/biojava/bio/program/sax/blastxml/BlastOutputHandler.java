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
 * Class to parse NCBI Blast-XML output.
 * <p>
 * WARNING: when Blast is supplied with multiple sequences in a single
 * FASTA file, it generates horribly-formed XML code.  This occurs because
 * it generates a complete XML document per query sequence including <?XML> 
 * lines and concatenates all of these into one output file.  This document
 * is not well-formed and cannot be parsed by standard XML parsers.
 * <p>
 * You have various options to solve this.  You can parse the output file
 * externally, separating the component documents and feeding each individually to
 * this parser.  Alternatively, you could strip out the <?XML> lines and wrap the
 * entire output in a fake <blast_aggregator> element.  This can be
 * done in Linux (and other Unixen) with :-
 * <p>
 * <pre>
 * #!/bin/sh
 * # Converts a Blast XML output to something vaguely well-formed
 * # for parsing.
 * # Use: blast_aggregate <XML output> <editted file>
 *
 * # strips all &lt;?xml&gt; and &lt;!DOCTYPE&gt; tags
 * # encapsulates the multiple &lt;BlastOutput&gt; elements into &lt;blast_aggregator&gt;
 *
 * sed '/&gt;?xml/d' $1 | sed '/&lt;!DOCTYPE/d' | sed '1i\
 * &lt;blast_aggregator&gt;
 * $a\
 * &lt;/blast_aggregator&gt;' > $2
 *</pre>
 * <p>
 * The resultant file can then be parsed with the BlastAggregator object.
 *
 * @author David Huen
 */
class BlastOutputHandler
    extends StAXFeatureHandler
{
    // create static factory class that makes an instance
    // of this class.
    public final static StAXHandlerFactory BLASTOUTPUT_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new BlastOutputHandler(staxenv);
            }
        };

    // class variables
    private String program = null;
    private String version = null;
    private String databaseId = null;
    private String queryId = null;
    private String queryDef = null;

    /**
     * If set, the output is wrapped in a <biojava:BlastLikeDataSetCollection>.
     * <p>
     * If not, it is left as a <biojava:BlastLikeDataSet>.
     * <p>
     * Default is true.
     */
    boolean wrap = true;

    // constructor when this is a element class in a document
    public BlastOutputHandler(StAXFeatureHandler staxenv)
    {
        super(staxenv);
//        System.out.println("BlastOutputHandler staxenv " + staxenv);
        // initialise delegation
        initDelegation();
    }

    // constructor
    private void initDelegation()
    {
        // delegate handling of <BlastOutput_param>
//        super.addHandler(new ElementRecognizer.ByLocalName("BlastOutput_param"),
//            BlastOutputParamHandler.BLAST_OUTPUT_PARAM_HANDLER_FACTORY);

        // delegate handling of <BlastOutput_program>
        super.addHandler(new ElementRecognizer.ByLocalName("BlastOutput_program"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s) throws SAXException {
                            program = s.trim();

                            // at this point, I can set the sequence types
                            // note that the sequence type here is the form
                            // in which blast DISPLAYS its output, not the
                            // sequence type of either query or target.
                            if (program.equals("blastn")) {
                                querySequenceType = "dna";
                                hitSequenceType = "dna";
                            }
                            else if (program.equals("blastp")) {
                                querySequenceType = "protein";
                                hitSequenceType = "protein";
                            }
                            else if (program.equals("blastx")) {
                                // nucleotide query translated in all frames 
                                // against protein database.
                                querySequenceType = "protein";
                                hitSequenceType = "protein";
                            }
                            else if (program.equals("tblastn")) {
                                // protein query against dna database in all frames
                                // hit frame is displayed only, no query frame
                                // irrespective of frame, both sequences displayed
                                // in increasing seq DNA coordinates (by from-to).
                                querySequenceType = "protein";
                                hitSequenceType = "protein";
                            }
                            else if (program.equals("tblastx")) {
                                // dna query translated in all frames against
                                // dna database in all frames
                                // irrespective of frame, both sequences displayed
                                // in increasing seq DNA coordinates.
                                querySequenceType = "protein";
                                hitSequenceType = "protein";
                            }
                            else throw new SAXException("unknown BLAST program.");
                        }
                    };
                }
            }
        );

        // delegate handling of <BlastOutput_version>
        super.addHandler(new ElementRecognizer.ByLocalName("BlastOutput_version"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {

                        public void setStringValue(String s) throws SAXException {
                            version = s.trim();
                        }
                        
                        public void startElement(
                            String nsURI,
                            String localName,
                            String qName,
                            Attributes attrs,
                            DelegationManager dm)
                            throws SAXException
                        {
                            // now generate my own start element
                            super.startElement(nsURI, localName, qName, attrs, dm);
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

                            super.endElement(nsURI, localName, qName, handler);

                            // generate the start of <BlastLikeDataSet> here.

                            // generate attributes
                            AttributesImpl bldsAttrs = new AttributesImpl();

//                            System.out.println("program, version " + program + " " + version);
                            if ((program != null) && (version != null)) {
                                bldsAttrs.addAttribute(biojavaUri, "program", "program", CDATA, program);
                                bldsAttrs.addAttribute(biojavaUri, "version", "version", CDATA, version);
                            }

                            listener.startElement(biojavaUri, "BlastLikeDataSet", biojavaUri + ":BlastLikeDataSet", bldsAttrs);

                            // generate start of header
                            listener.startElement(biojavaUri, "Header", biojavaUri + ":Header", new AttributesImpl());

                            // we don't have raw output but it is compulsory
                            listener.startElement(biojavaUri, "RawOutput", biojavaUri + ":RawOutput", new AttributesImpl());
                            listener.endElement(biojavaUri, "RawOutput", biojavaUri + ":RawOutput");
                        }
                    };
                }
            }
        );

        // delegate handling of <BlastOutput_reference>
//        super.addHandler(new ElementRecognizer.ByLocalName("BlastOutput_reference"),
//            SearchPropertyHandler.SEARCH_PROPERTY_HANDLER_FACTORY);

        // delegate handling of <BlastOutput_db>
        super.addHandler(new ElementRecognizer.ByLocalName("BlastOutput_db"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s)  throws SAXException {
                            databaseId = s.trim();
                        }
                    };
                }
            }
        );

        // delegate handling of <BlastOutput_query-ID>
        super.addHandler(new ElementRecognizer.ByLocalName("BlastOutput_query-ID"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s)  throws SAXException {
                            queryId = s.trim();
                        }
                    };
                }
            }
        );

        // delegate handling of <BlastOutput_query-def>
        super.addHandler(new ElementRecognizer.ByLocalName("BlastOutput_query-def"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StringElementHandlerBase() {
                        public void setStringValue(String s)  throws SAXException {
                            queryDef = s.trim();
                        }
                    };
                }
            }
        );

        // delegate handling of <BlastOutput_query-len>
//        super.addHandler(new ElementRecognizer.ByLocalName(""),
//            SearchPropertyHandler.SEARCH_PROPERTY_HANDLER_FACTORY);

        // delegate handling of <BlastOutput_iterations>
        super.addHandler(new ElementRecognizer.ByLocalName("BlastOutput_iterations"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new BlastOutputIterationsHandler(staxenv) {
                        public void startElementHandler(
                            String nsURI,
                            String localName,
                            String qName,
                            Attributes attrs)
                            throws SAXException
                        {
                            // necessary as staxenv cannot be final and therefore
                            // staxenv.listener cannot be accessed from inner class
                            ContentHandler listener = getListener();
//                            System.out.println("about to start outputting QueryId. listener is " + listener);

                            // the DatabaseId and QueryId elements are generated
                            // in reversed order so I cannot generate them on-the-fly.

                            // generate <QueryId> if required
                            if (queryId != null) {
                                AttributesImpl queryAttrs = new AttributesImpl();

                                queryAttrs.addAttribute(biojavaUri, "id", "id", CDATA, queryId);
                                queryAttrs.addAttribute(biojavaUri, "metadata", "metadata", CDATA, "none");
                                listener.startElement(biojavaUri, "QueryId", biojavaUri + ":QueryId", queryAttrs);
                                listener.endElement(biojavaUri, "QueryId", biojavaUri + ":QueryId");
                            }

                            // generate <QueryDescription> if required
                            if (queryDef != null) {
                                listener.startElement(biojavaUri, "QueryDescription", biojavaUri + ":QueryDescription", new AttributesImpl());
                                listener.characters(queryDef.toCharArray(), 0, queryDef.length());
                                listener.endElement(biojavaUri, "QueryDescription", biojavaUri + ":QueryDescription");
                            }

                            if (databaseId != null) {
                                AttributesImpl dbAttrs = new AttributesImpl();

                                dbAttrs.addAttribute(biojavaUri, "id", "id", CDATA, databaseId);
                                dbAttrs.addAttribute(biojavaUri, "metadata", "metadata", CDATA, "none");
                                listener.startElement(biojavaUri, "DatabaseId", biojavaUri + ":DatabaseId", dbAttrs);
                                listener.endElement(biojavaUri, "DatabaseId", biojavaUri + ":DatabaseId");
                            }

                            // generate end of <Header>
                            listener.endElement(biojavaUri, "Header", biojavaUri + ":Header");

                            // now I generate my own start element: does nothing.
                            super.startElementHandler(nsURI, localName, qName, attrs);
                        }
                    };
                }
            }
        );
    }

    void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler handler)
             throws SAXException
    {
        // generate end of <biojava:BlastLikeDataSet>
        staxenv.listener.endElement(biojavaUri, "BlastLikeDataSet", biojavaUri + ":BlastLikeDataSet");
    }
}
