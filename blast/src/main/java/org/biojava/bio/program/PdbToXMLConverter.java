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

import org.biojava.bio.program.sax.PdbSAXParser;
import org.biojava.bio.program.xml.SimpleXMLEmitter;
import org.xml.sax.ContentHandler;
import org.xml.sax.XMLReader;

/**
 * <p>
 * A class that converts Protein Data Bank (PDB) to
 * XML that will validate against the biojava:MacromolecularStructure DTD.
 * <p>
 * <b>Note this code is experimental and subject to change without notice.
 * </b>
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
 *</ul>
 *
 * This code was first released to the biojava.org project, July 2000.
 *
 * @author Cambridge Antibody Technology (CAT)
 * @version 0.1
 * 
 * @see org.biojava.bio.program.sax.BlastLikeSAXParser
 * @see SimpleXMLEmitter
 */
public class PdbToXMLConverter {

    private String            oInput;
    private XMLReader         oParser;

    /**
     * Creates a new <code>BlastToXMLConverter</code> instance.
     *
     */
    public PdbToXMLConverter(String poInput) {
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
	oParser = (XMLReader) new PdbSAXParser();


	/**
	 * Dynamically change configuration of the parser
	 * in regard of Namespace support. Here,
	 * the xml.org/features/namespaces feature is simply "reset"
	 * to its default value for SAX2.
	 * The xml.org/features/namespaces-prefixes feature is
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
	oParser.parse(oInput);  

	System.out.println();
	
    }

}
