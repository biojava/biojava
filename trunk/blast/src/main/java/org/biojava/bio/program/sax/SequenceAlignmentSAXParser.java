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

import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;

/**
 * A SAX2 parser for dealing with a sequence alignments.  The format
 * of any given alignment is automatically detected (e.g. ClustalW,
 * Needle).
 *
 * Supported alignment formats are:
 * <ul>
 * <li>ClustalW
 * <li>Needle
 * </ul>
 *
 * Copyright &copy; 2000-2002 Cambridge Antibody Technology.
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
public class SequenceAlignmentSAXParser extends AbstractNativeAppSAXParser {


    private static final int        CLUSTALW           = 2;
    private static final int        NEEDLE             = 3;

    private int                     iAlignmentType     = -1;
    /**
     * Initialises internal state
     * Sets namespace prefix to "biojava"
     */
    public SequenceAlignmentSAXParser() {
	this.setNamespacePrefix("biojava");
    }

    /**
     * Describe 'parse' method here.
     *
     * @param poSource	 -
     */
    public void parse(InputSource poSource ) 
	throws IOException,SAXException {

	BufferedReader            oContents;
	String                    oLine = null;

	//Use method form superclass
	oContents = this.getContentStream(poSource);

	oContents.mark(500);

	oLine = null;
	try {
	    oLine = oContents.readLine();
	} catch (java.io.IOException x) {
	    System.out.println(x.getMessage());
	    System.out.println("Stream read interupted");
	} // end try/catch

	//at end of stream...

        //System.out.println(oLine);
	//Choose SAX Parser

	XMLReader oChosenSAXParser = null;

	if (oLine.startsWith("CLUSTAL W")) {
	    iAlignmentType = CLUSTALW;
	}

	if (oLine.startsWith("Global: ")) {
	    iAlignmentType = NEEDLE;

	}

	switch (iAlignmentType) {
	case CLUSTALW:
	    //	    System.out.println("FOUND CLUSTALW");
	    oChosenSAXParser = new ClustalWAlignmentSAXParser();
	    break;
	case NEEDLE:
	    //System.out.println("FOUND NEEDLE");
	    oChosenSAXParser = new NeedleAlignmentSAXParser();
	    break;
	default:
	    //	    System.out.println("ALIGNMENT TYPE NOT FOUND");
	    //TODO SHOULD THROW AN EXCEPTION!
	    break;
	}

	oContents.reset();

	XMLReader oParser = oChosenSAXParser;

	oParser.setContentHandler(oHandler);

	/*
	 * Parse stream with appropriate parser
         */
        oParser.parse(new InputSource(oContents));
	
	oContents.close();

    }

}
