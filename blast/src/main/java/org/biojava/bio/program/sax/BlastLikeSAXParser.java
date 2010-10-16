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
package org.biojava.bio.program.sax;

import java.io.BufferedReader;
import java.io.IOException;

import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;

/**
 * A facade class allowing for direct SAX2-like parsing of the native
 * output from Blast-like bioinformatics software.  Because the parser is SAX2
 * compliant, application writers can simply pass XML ContentHandlers
 * to the parser in order to receive notifcation of SAX2 events.
 * <p>
 * The SAX2 events produced are as if the input to the parser was
 * an XML file validating against the biojava BlastLikeDataSetCollection DTD.
 * There is no requirement for an intermediate conversion of native output to
 * XML format. An application of the parsing framework, however, is to
 * create XML format files from native output files.
 * <p>
 * The biojava Blast-like parsing framework is designed to uses minimal 
 * memory,so that in principle, extremely large native outputs can be
 * parsed and XML ContentHandlers can listen only for small amounts of
 * information.
 * <p>
 * The framework currently supports parsing of native output from
 * the following bioinformatics programs. Please note that if
 * you are using different versions of NCBI or WU Blast to those
 * listed below, it is worth considering trying setting the parsing 
 * mode to Lazy, which means parsing will be attempted if the program
 * is recognised, regardless of version.
 * <ul>
 * <li>NCBI Blast version 2.0.11
 * <li>NCBI Blast version 2.2.2
 * <li>NCBI Blast version 2.2.3
 * <li>WU-Blast version 2.0a19mp-washu
 * <li>HMMER 2.1.1 hmmsearch 
 * </ul>
 * Planned addition support
 * <ul>
 * <li> Support for HMMER hmmpfam almost there but not fully tested
 * </ul>
 * <p>
 * <p>
 * <b>Notes to SAX driver writers</b>
 * <p>
 * The framework that this parser is built on is designed to be
 * extensible with support for both different pieces of software
 * (<i>i.e.</i> not just software that produces Blast-like output),
 * and multiple versions of programs.
 * <p>
 * This class inherits from the 
 * org.biojava.bio.program.sax.AbstractNativeAppSAXParser
 * abstract base class.  The abstract base class is a good place to
 * start looking if you want to write new native application SAX2 parsers.
 * This and releated classes have only package-level visibility.
 * Typically, application writers are expected to provide a facade class
 * in this package (similar to the current class) to allow
 * users access to functionality.
 * <p>
 * NB Support for InputSource is not complete due to the fact
 * that URLs are not resolved and cannot, therefore, be used
 * as an InputSource.  System pathnames, ByteStreams and CharacterStreams,
 * however, are all supported.
 * <p>
 *
 * Copyright &copy; 2000 Cambridge Antibody Technology.
 * 
 * <p>
 * Primary author -<ul>
 * <li>Simon Brocklehurst (CAT)
 * </ul>
 * Other authors  -<ul>
 * <li>Tim Dilks          (CAT)
 * <li>Colin Hardman      (CAT)
 * <li>Stuart Johnston    (CAT)
 * <li>Mathieu Wiepert    (Mayo Foundation)
 * <li>Travis Banks
 *</ul>
 *
 * @author Cambridge Antibody Technology (CAT)
 * @author Travis Banks
 * @version 1.0
 *
 * @see org.biojava.bio.program.BlastLikeToXMLConverter
 */
public class BlastLikeSAXParser extends AbstractNativeAppSAXParser {

    private BlastLikeVersionSupport oVersion  = new BlastLikeVersionSupport();
    private BlastSAXParser          oBlast; 

    private AttributesImpl          oAtts     = new AttributesImpl();
    private QName                   oAttQName = new QName(this);     
    private boolean                 tValidFormat  = false;

    private static final int        STARTUP                   = 0;
    private static final int        INSIDE_FILE               = 1;

    private String                  oStoredLine = null;

    /**
     * Initialises SAXParser, and sets default namespace prefix
     * to "biojava".
     */
    public BlastLikeSAXParser() {
        this.changeState(STARTUP);

        //centralised setting of namespace prefix
        //the setting is cascaded everywhere else
        this.setNamespacePrefix("biojava");
        this.addPrefixMapping("biojava","http://www.biojava.org");
        
        oVersion.setMode(BlastLikeVersionSupport.LAZY);
    }

    /**
     * <code>parse</code> initiates the parsing operation.
     *
     * @param poSource an <code>InputSource</code>.
     * @exception IOException if an error occurs.
     * @exception SAXException if an error occurs.
     */
    public void parse(InputSource poSource ) 
    throws IOException, SAXException {

        BufferedReader            oContents;
        String                    oLine;

        this.changeState(STARTUP);

        //Use method form superclass
        oContents = this.getContentStream(poSource);
        //This sets contentHandler document for XSLT
        this.getContentHandler().startDocument();
        
        try {
            // loop over file
            oLine = oContents.readLine();
            while (oLine != null) {
                //interpret line and send messages accordingly        
                this.interpret(oContents,oLine);
                //do extra interpretation of lines reached by subparser
                //objects
                if (iState == INSIDE_FILE) {
                    oLine = oStoredLine;
                    if (oStoredLine != null) {
                        this.interpret(oContents,oLine);
                    }
                } else {
                    oLine = oContents.readLine();
                }

            } // end while
        } catch (IOException x) {
            System.out.println(x.getMessage());
            System.out.println("File read interrupted");
        } // end try/catch

        //at end of file...
        oContents.close();

        if (!tValidFormat) {
            throw (new SAXException("Could not recognise the format " +
            "of this file as one supported by the framework."));
        }

        this.endElement(new QName(this,
                this.prefix("BlastLikeDataSetCollection")));
    }

    /**
     * This is the default, parsing will be attempted only if both
     * the program e.g. NCBI BlastP, and a particular version 
     * are recognised as bsing supported.
     *
     */
    public void setModeStrict() {
        oVersion.setMode(BlastLikeVersionSupport.STRICT);
    }

    /**
     * Setting the mode to lazy means that, if the program is recognised,
     * e.g. WU-TBlastX, then parsing will be attempted even if 
     * the particular version is not recognised. Using this option
     * is more likely to result in erroneous parsing than if the
     * strict mode is used.
     *
     */
    public void setModeLazy() {
        oVersion.setMode(BlastLikeVersionSupport.LAZY);
    }

    /**
     * Deal with line according to state parser is in.
     *
     * @param poLine     A line of Blast output
     */
    private void interpret(BufferedReader poContents, String poLine)
    throws SAXException {
        //For a brand new collection,
        //check for the start of a new BlastDataSet
        if (iState == STARTUP) {
            //look for characteristic of start of dataset
            if (oVersion.isStartOfDataSet(poLine)) {
                //For GCG, oVersion is set as an indicator to get 
                //program info from the second line
                if (oVersion.getProgram() == BlastLikeVersionSupport.GCG) { 
                    //if GCG, skip to next line to get program info
                    try {
                        poLine = poContents.readLine ();
                    } catch (java.io.IOException x) {
                        System.out.println(x.getMessage());
                        System.out.println("File read interrupted");
                        throw (new SAXException("Error parsing GCG File"));
                    } // end try/catch
                }

                tValidFormat = oVersion.assignProgramAndVersion(poLine);

                if (!oVersion.isSupported()) {
                    throw (new SAXException("Program "
                            + oVersion.getProgramString()
                            + " Version "
                            + oVersion.getVersionString()
                            + " is not supported by the biojava blast-like "
                            + "parsing framework"));
                }

                oAtts.clear();
                oAttQName.setQName("xmlns");
                //check if namespace configuration means attribute
                //should not be reported.
                if (!oAttQName.getLocalName().equals("")) {
                    oAtts.addAttribute(oAttQName.getURI(),
                            oAttQName.getLocalName(),
                            oAttQName.getQName(),
                            "CDATA","");
                }

                oAttQName.setQName("xmlns:biojava");
                //check if namespace configuration means attribute
                //should not be reported.
                if (!oAttQName.getLocalName().equals("")) {
                    oAtts.addAttribute(oAttQName.getURI(),
                            oAttQName.getLocalName(),
                            oAttQName.getQName(),
                            "CDATA","http://www.biojava.org");
                }

                this.startElement(new QName(this,
                        this.prefix("BlastLikeDataSetCollection")),
                        (Attributes)oAtts);

                this.onNewDataSet(poContents,poLine);
                return;
            }
        }   //End check for the start of a new BlastDataSet

        if (iState == INSIDE_FILE) {
            //look for characteristic of start of dataset
            if (oVersion.isStartOfDataSet(poLine)) {

                tValidFormat = oVersion.assignProgramAndVersion(poLine);

                this.onNewDataSet(poContents,poLine);

                return;
            }
        }   //End check for the start of a new BlastDataSet
    }

    /**
     *
     * When this method is called, the line will look something line:
     *
     * BLASTN 2.0.11 [Jan-20-2000]
     *
     * The above would be parsed to program blastn, and version number.
     *
     * @param poLine     -
     */
    private void onNewDataSet(BufferedReader poContents, String poLine)
    throws SAXException {

        //choose according to version...
        oBlast = new BlastSAXParser(oVersion,this.getNamespacePrefix());
        String oLine="";

        //Parse Contents stream up to end of a single BlastDataSet.
        oBlast.setContentHandler(oHandler);
        while(oLine!=null) {
        	oLine = oBlast.parse(poContents,poLine);
        }

        this.changeState(INSIDE_FILE);

        //now make sure to interpret the line the BlastSAXParser returned from
        //in top-level parse method.
        if (oLine != null) {
            oStoredLine = oLine;
            return;
            //           this.interpret(poContents,oLine);
        } else {
            //here if at the EOF
            oStoredLine = null;
            return;
        }
    }
}


