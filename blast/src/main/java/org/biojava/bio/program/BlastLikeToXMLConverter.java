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
package org.biojava.bio.program;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;

import org.biojava.bio.program.sax.BlastLikeSAXParser;
import org.biojava.bio.program.xml.SimpleXMLEmitter;
import org.xml.sax.ContentHandler;
import org.xml.sax.InputSource;
import org.xml.sax.XMLReader;

/**
 * <p>
 * A class that converts the raw output from a variety of bioinformatics
 * software and converts it to XML that will validate against the
 * biojava:BlastLikeDataSetCollection DTD.
 * <p>
 * For applications supported, please the documentation for the
 * BlastLikeSAXParser.
 * <p>
 * Examination of the source code of this application also serves as
 * demonstration of the simplicity of using the biojava blast-like SAX2 
 * parsing framework.  The main functionality of the application is
 * simply built from the following code, <i>viz.</i>:
 * <pre>
 *
 *      !**
 *       * The following code creates a parser for native output
 *       * from BlastLike programs. That is,
 *       * Create a SAX2 Parser that takes the native output
 *       * from blast-like bioinformatics software.
 *       *!
 *       <font color="#0000FF">
 *        XMLReader oParser =
 *       (XMLReader) new BlastLikeSAXParser();
 *       </font>
 *     !**
 *       * Namespace support controls the way in which
 *       * XML elements are reported. In XML, when an element
 *       * looks something like <biojava:Hit> then,
 *       * the part before the colon, i.e. biojava is the namespace,
 *       * and the part after the colon i.e. Hit is the Local name.
 *       * The full "biojava:Hit" name is termed the Qualified Name (QNames).
 *       * By default SAX2 parsers report Local Names, in this
 *       * example, we decided we wanted to make the parser report QNames.
 *       *
 *       * If you don't want to change default namespace support, you
 *       * can ignore the next piece of code.
 *       *
 *       * Dynamically change configuration of the parser
 *       * in regard of Namespace support. Here,
 *       * the xml.org/features/namespaces feature is simply "reset"
 *       * to its default value for SAX2.
 *       * The xml.org/features/namespaces-prefixes feature is
 *       * also set to true.  This is to ensure that xmlns attributes
 *       * are reported by the parser. These are required because we want
 *       * to configure the XMLEmitter to output qualified names (see below).
 *       *!
 *      <font color="#0000FF">
 *      try {
 *      oParser.setFeature("http://xml.org/sax/features/namespaces",true);
 *      oParser.setFeature(
 *              "http://xml.org/sax/features/namespace-prefixes",true);
 *
 *      } catch (Exception e) {
 *      }
 *      </font>
 *
 *      !**
 *       * Having selected the parser, we now want to
 *       * choose an object to deal with the SAX2 events
 *       * that the parser produces. This is the class
 *       * that you would normally write yourself to deal
 *       * with particular events you are interested in.
 *       * This class implements the ContentHandler - usually,
 *       * you would inherit from a SAX2 helper class that
 *       * implements this interface for you.
 *       *
 *       * Create an XML ContentHandler. This
 *       * implementation of the DocumentHandler
 *       * interface simply outputs nicely formatted
 *       * XML. Passing a true value to the SimpleXMLEmitter
 *       * constructor instructs the ContentHandler
 *       * to take QNames from the SAXParser, rather
 *       * than LocalNames.
 *       *
 *      <font color="#0000FF">
 *      ContentHandler oHandler  = 
 *      (ContentHandler) new SimpleXMLEmitter(true);
 *      </font>
 *
 *      !**
 *       * Now, link the Parser and the ContentHandler.
 *       *
 *       * Give the parser a reference to the ContentHandler
 *       * so that it can send SAX2 mesagges.
 *       *!
 *      <font color="#0000FF">
 *      oParser.setContentHandler(oHandler);
 *      </font>
 *      !**
 *       * Finally, parse your Blast-like output.
 *       *
 *       * Now make the Parser parse the output from the
 *       * blast-like software and emit XML as specificed
 *       * by the ContentHandler.
 *       *!
 *      <font color="#0000FF">
 *      oParser.parse(oInput);  
 *      </font>
 * </pre>
 *
 * <p>
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
 *</ul>
 *
 *
 * @author Cambridge Antibody Technology (CAT)
 * @version 1.0
 * 
 * @see BlastLikeSAXParser
 * @see SimpleXMLEmitter
 */
public class BlastLikeToXMLConverter {

    private String            oInput;
    private XMLReader         oParser;
    private boolean           tStrict         = true;

    /**
     * Creates a new <code>BlastToXMLConverter</code> instance.
     *
     */
    public BlastLikeToXMLConverter(String poInput) {
    oInput = poInput;
    }

    public void convert() throws java.io.IOException,
                                 org.xml.sax.SAXException {

    //Access functionality of biojava classes through
    //standard org.xml.sax interfaces...

    /**
     * Create a SAX Parser that takes the native output
     * from blast-like bioinformatics software.
     */
    oParser = (XMLReader) new BlastLikeSAXParser();

    if (tStrict) {
        ((BlastLikeSAXParser) oParser).setModeStrict();
    } else {
        ((BlastLikeSAXParser) oParser).setModeLazy();
    }
    /**
     * Dynamically change configuration of the parser
     * in regard of Namespace support. Here,
     * the xml.org/sax/features/namespaces feature is simply "reset"
     * to its default value for SAX2.
     * The xml.org/sax/features/namespaces-prefixes feature is
     * also set to true.  This is to ensure that xmlns attributes
     * are reported by the parser. These are required because we want
     * to configure the XMLEmitter to output qualified names (see below).
     */
    try {
        oParser.setFeature("http://xml.org/sax/features/namespaces",true);
        oParser.setFeature("http://xml.org/sax/features/namespace-prefixes",
                   true);

    } catch (Exception e) {
        //If an illegal conmbination of features is chosen,
        //roll back to default settings. Output a warning,
        //even though this might mess up the output...
        System.out.println("WARNING: ignoring attempt to set illegal " +
                   "combination of parser features");
        System.out.println(e);
    }
    /**
     * Create an XML ContentHandler. This
     * implementation of the DocumentHandler
     * interface simple outputs nicely formatted
     * XML. Passing a true value to the SimpleXMLEmitter
     * constructor instructs the ContentHandler
     * to take QNames from the SAXParser, rather
     * than LocalNames.
     */
    ContentHandler oHandler  = 
        (ContentHandler) new SimpleXMLEmitter(true);

    /**
     * Give the parser a reference to the ContentHandler
     * so that it can send SAX2 mesagges.
     */
    oParser.setContentHandler(oHandler);
    /**
     * Now make the Parser parse the output from the
     * blast-like software and emit XML as specificed
     * by the DocumentHandler.
     */
    //Test direct specification of URI
    //oParser.parse(oInput);  

    
    
    //Test direct specification of URI via InputSource
    //oParser.parse(new InputSource(oInput));  



    FileInputStream           oInputFileStream;
    BufferedReader            oContents;

    //Test parsing using ByteSteam as InputSource
        // Open file and read all lines from file sequentially
//         try{
//             oInputFileStream = new FileInputStream(oInput);
//             // create input stream

//      oParser.parse(new InputSource(oInputFileStream));

//         } catch (java.io.FileNotFoundException x) {
//             System.out.println(x.getMessage());
//             System.out.println("Couldn't open file");
//             System.exit(0);
//         }


    //Test parsing using CharacterStream as InputSource
        // Open file and read all lines from file sequentially
        try{
             oInputFileStream = new FileInputStream(oInput);
             // create input stream
             oContents = new
                 BufferedReader(new InputStreamReader(oInputFileStream));

        oParser.parse(new InputSource(oContents));

         } catch (java.io.FileNotFoundException x) {
             System.out.println(x.getMessage());
             System.out.println("Couldn't open file");
             System.exit(0);
     }


    System.out.println();
    
    }

    public void setModeStrict() {
    tStrict = true;
    }
    public void setModeLazy() {
    tStrict = false;

    }


}
