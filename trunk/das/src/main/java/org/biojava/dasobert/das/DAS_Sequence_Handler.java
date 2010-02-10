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
 * Created on 19.03.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.dasobert.das;


import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.Attributes            ;
import java.util.logging.*                ;

/**
 * a class that parses the XML response of a DAS - sequence command.
 * @author Andreas Prlic
 *
 */
public class DAS_Sequence_Handler extends DefaultHandler {

	StringBuffer sequence ;
	int length ;

	int maxLength;
	String version;

	boolean dna_flag; 
	/**
	 * 
	 */
	public DAS_Sequence_Handler() {
		super();

		sequence = new StringBuffer() ;
		length = 0;
		dna_flag = false ;
		maxLength = -1;
		version = "";
	}



	/** set a maximum length of sequence that should be loaded
	 * default: -1. if -1 no length restriction is being supplied
	 * @return the maximum length or -1 if no restriction
	 */
	public int getMaxLength() {
		return maxLength;
	}



	/** set a maximum length of sequence that should be loaded
	 * default: -1. if -1 no length restriction is being supplied
	 * @param maxLength the maximum length or -1 if unrestricted
	 */
	public void setMaxLength(int maxLength) {
		this.maxLength = maxLength;
	}




	public void startElement (String uri, String name, String qName, Attributes atts){

		if ( qName.equals("SEQUENCE")){
			version = atts.getValue("version");
			String lenstr 	= atts.getValue("stop");
			length = Integer.parseInt(lenstr);
			dna_flag = true ;
		}

	}

	public void characters (char ch[], int start, int length){		

		if (maxLength > 0)
			if ( sequence.length() > maxLength)
				return;

		if (dna_flag) 
			for (int i = start; i < start + length; i++) {

				// all sorts of characters can be found in "seqeunces" ... ignore them...
				switch (ch[i]) {
				case '\\':
					//System.out.print("\\\\");
					break;
				case '"':
					//System.out.print("\\\"");
					break;
				case '\n':
					//System.out.print("\\n");
					break;
				case '\r':
					//System.out.print("\\r");
					break;
				case '\t':
					//System.out.print("\\t");
					break;
				case ' ':
					break;
				default:
					sequence = sequence.append(ch[i]);

				break;
				}
			}

	}

	public String get_sequence() {

		if ( maxLength < 0) {		
			if ( length != sequence.length()) {	
				Logger logger  = Logger.getLogger("org.biojava.spice");
				logger.warning("Sequence does not match specified length!");
			}
		}

		return sequence.toString();
	}

	public String getVersion() {
		return version;
	}

}
