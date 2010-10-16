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
import java.util.StringTokenizer;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;

/**
 * A reusable class for parsing Hmmr Domain Summary Section
 *
 * Primary author -
 *                 Colin Hardman      (CAT)
 * Other authors  -
 *                 Tim Dilks          (CAT)
 *                 Simon Brocklehurst (CAT)
 *                 Stuart Johnston    (CAT)
 *                 Lawerence Bower    (CAT)
 *                 Derek Crockford    (CAT)
 *                 Neil Benn          (CAT)
 *
 * Copyright 2001 Cambridge Antibody Technology Group plc.
 * 
 *
 * This code released to the biojava project, May 2001
 * under the LGPL license.
 *
 * @author Cambridge Antibody Technology Group plc
 * @version 0.1
 *
 */
class DomainSectionSAXParser extends AbstractNativeAppSAXParser {


    private BufferedReader       oContents;
    private AttributesImpl       oAtts              = new AttributesImpl();
    private QName                oAttQName          = new QName(this);     
    private String               oLine;

    /**
     * Creates a new domain section parser.
     *
     * @param poVersion <code>BlastLikeVersionSupport</code>
     * @param poNamespacePrefix - the namespace prefix 
     */
    DomainSectionSAXParser( BlastLikeVersionSupport poVersion,
			    String poNamespacePrefix ) {
	this.setNamespacePrefix(poNamespacePrefix);
	//For XSLT Parser Compliance
	this.addPrefixMapping("biojava","http://www.biojava.org");

    }

    /**
     * Parse the buffer from the current position, until domain hits
     * are parsed.
     *
     * @param poContents <code>BufferedReader</code> to parse.
     * @param poLine - the value of the current line.
     * @exception SAXException if an error occurs
     */
    public void parse( BufferedReader poContents, String poLine ) 
	throws SAXException {

	oContents = poContents;

	try {
	    oLine = poLine;

	    while (! oLine.startsWith( "----" )) { // skip headers
		oLine = oContents.readLine();
	    }
	    oLine = oContents.readLine();

	    if ( oLine.trim().equals( "[no hits above thresholds]" ) ) {
		oLine = oContents.readLine();
		return;
	    }

	    while (!(oLine.trim().equals(""))) {

		//		System.out.println( "->" + oLine + "<-" );

		StringTokenizer st = new StringTokenizer( oLine );
		oAtts.clear();
		oAttQName.setQName("modelId");
		oAtts.addAttribute(oAttQName.getURI(),
				   oAttQName.getLocalName(),
				   oAttQName.getQName(),
				   "CDATA", st.nextToken() );

		StringTokenizer st2 = new StringTokenizer
		    ( st.nextToken(), "/");

		oAttQName.setQName("domainPosition");
		oAtts.addAttribute(oAttQName.getURI(),
				   oAttQName.getLocalName(),
				   oAttQName.getQName(),
				   "CDATA", st2.nextToken() );

		oAttQName.setQName("sequenceFrom");
		oAtts.addAttribute(oAttQName.getURI(),
				   oAttQName.getLocalName(),
				   oAttQName.getQName(),
				   "CDATA", st.nextToken() );

		oAttQName.setQName("sequenceTo");
		oAtts.addAttribute(oAttQName.getURI(),
				   oAttQName.getLocalName(),
				   oAttQName.getQName(),
				   "CDATA", st.nextToken() );
		

		String oTemp       = st.nextToken();

		oAttQName.setQName("startPositionOfSequence");

		if(oTemp.charAt(0) == '.') {
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "NMTOKEN", "internal" );

		}else{
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "NMTOKEN", "end" );
		}
		oAttQName.setQName("endPositionOfSequence");

		if(oTemp.charAt(1) == '.') {
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "NMTOKEN", "internal" );
		}else{
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "NMTOKEN", "end" );
		}

		oAttQName.setQName("hmmFrom");
		oAtts.addAttribute(oAttQName.getURI(),
				   oAttQName.getLocalName(),
				   oAttQName.getQName(),
				   "CDATA", st.nextToken() );

		oAttQName.setQName("hmmTo");
		oAtts.addAttribute(oAttQName.getURI(),
				   oAttQName.getLocalName(),
				   oAttQName.getQName(),
				   "CDATA", st.nextToken() );

		oTemp       = st.nextToken();

		oAttQName.setQName("startPositionOfModel");

		if(oTemp.charAt(0) == '.') {
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "NMTOKEN", "internal" );

		}else{
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "NMTOKEN", "end" );
		}
		oAttQName.setQName("endPositionOfModel");

		if(oTemp.charAt(1) == '.') {
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "NMTOKEN", "internal" );
		}else{
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "NMTOKEN", "end" );
		}


		oAttQName.setQName("score");
		oAtts.addAttribute(oAttQName.getURI(),
				   oAttQName.getLocalName(),
				   oAttQName.getQName(),
				   "CDATA", st.nextToken() );

		oAttQName.setQName("expectValue");
		oAtts.addAttribute(oAttQName.getURI(),
				   oAttQName.getLocalName(),
				   oAttQName.getQName(),
				   "CDATA", st.nextToken() );

		this.startElement(new QName(this,this.prefix("DomainHit")),
				  (Attributes)oAtts);
		this.endElement(new QName(this,this.prefix("DomainHit")));
		
		oLine = oContents.readLine();
	    } // end while in domain summary

	} catch (java.io.IOException x) {
	    System.out.println(x.getMessage());
	    System.out.println("File read interupted");
	} // end try/catch
    } 
}
