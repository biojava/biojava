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


import org.xml.sax.ContentHandler;
import org.xml.sax.XMLReader;

/**
 * Test application for use by SAX Parser writers.  Allows simple functions
 * of a SAX Parser to be developed/debugged without being dependent on other
 * packages.This is useful in complex development/build environments
 * <p>
 * Copyright &copy; 2001 Cambridge Antibody Technology.
 * 
 * <p>
 * Primary author -<ul>
 * <li>Simon Brocklehurst (CAT)
 * </ul>
 * Other authors  -<ul>
 * <li>Neil Benn          (CAT)
 * <li>Lawrence Bower     (CAT)
 * <li>Derek Crockford    (CAT)
 * <li>Tim Dilks          (CAT)
 * <li>Colin Hardman      (CAT)
 * <li>Stuart Johnston    (CAT)
 *</ul>
 *
 * @author Cambridge Antibody Technology (CAT)
 * @version 1.0
 *
 */
class GenericSAXParserTest {

    /**
     * Given the name of a SAXParser and a pathname, produce XML
     *
     * @param args a <code>String[]</code> representation of a pathname
     * @exception Exception if an error occurs
     */
    public static void main(String[] args) throws Exception {

	String oInput = null;
	String oClassName = null;

        // Catch wrong number of arguments or help requests

        if ( (args.length != 2)        ||
	     (args[0].equals("-help")) || 
	     (args[0].equals("-h")) ) {

	    System.out.println();
	    System.out.println();
	    System.out.println(
             "Given the name of a SAXParser (i.e. the name of the class, and");
	    System.out.println(
             "and the pathname of an example file you wish to parse, the");
	    System.out.println(
             "application outputs an XML representation of the data");
	    System.out.println();
	    System.out.println(
             "Usage: java GenericSAXParserTest <classname> <pathname>");

	    System.out.println();
	    System.out.println(
             "For example java " +
             "org.biojava.bio.program.sax.GenericSAxParserTest" +
	     "org.biojava.bio.program.sax.BlastLikeSAXParser");
	    System.out.println();
	    System.out.println();

            System.exit(1);
        }

	if (args.length == 2) {
	    oClassName = args[0];
	    oInput = args[1];
	}

        //Now the actual application...

	// Get hold of the chosen SAXParser

	XMLReader oChosenSAXParser = null;

	try{

	    oChosenSAXParser = 
		(XMLReader)Class.forName(oClassName).newInstance();

	}catch(InstantiationException ie){
	    ie.printStackTrace();
	}catch(IllegalAccessException iae){
	    iae.printStackTrace();
	}catch(ClassNotFoundException cnfe){
	    cnfe.printStackTrace();
	}

	/**
	 * Create a SAX Parser
	 */
	XMLReader oParser = oChosenSAXParser;

	/**
	 * Create an XML ContentHandler. This
	 * implementation of the DocumentHandler
	 * interface simply outputs nicely formatted
	 * XML. Passing a true value to the SimpleXMLEmitter
	 * constructor instructs the ContentHandler
	 * to take QNames from the SAXParser, rather
	 * than LocalNames.
	 */

	ContentHandler oHandler  = 
	    (ContentHandler) new SimpleXMLEmitterTestHelper();

	/*
	 * Give the parser a reference to the ContentHandler
	 * so that it can send SAX2 mesagges.
	 */
	oParser.setContentHandler(oHandler);

	/*
	* Parsing using file as input
        * Open file and read all lines from file sequentially
	*/
        oParser.parse(oInput);

	System.out.println();
	System.out.println("---------------------END------------------");
    }
}
