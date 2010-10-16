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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;

/**
 * A reusable class for parsing Detail
 * sections of Blast-like programs:
 *        oNCBI Blast
 *
 * Primary author -
 *                 Simon Brocklehurst (CAT)
 * Other authors  -
 *                 Tim Dilks          (CAT)
 *                 Colin Hardman      (CAT)
 *                 Stuart Johnston    (CAT)
 *                 Mathieu Wiepert    (Mayo Foundation)
 *                 Travis Banks		  (AAFC)
 *
 * Copyright 2000 Cambridge Antibody Technology Group plc.
 * 
 *
 * This code released to the biojava project, May 2000
 * under the LGPL license.
 *
 * @author Cambridge Antibody Technology Group plc
 * @author Greg Cox
 * @author Travis Banks
 * @version 0.1
 *
 */
final class HitSectionSAXParser extends AbstractNativeAppSAXParser {

	private BlastLikeAlignmentSAXParser oAlignmentParser;
	private BlastLikeVersionSupport     oVersion;

	private BufferedReader       oContents;
	private AttributesImpl       oAtts              = new AttributesImpl();
	private QName                oAttQName          = new QName(this);
	private char[]               aoChars;
	private char[]               aoLineSeparator;
	private ArrayList<String>    oGlobalEndSignals;
	private ArrayList            oBuffer            = new ArrayList();
	private ArrayList<String>    oAlignmentBuffer   = new ArrayList<String>();
	private StringBuffer         oStringBuffer      = new StringBuffer();
	private StringBuffer         oDescription       = new StringBuffer();
	private String               oLine;
	private HashMap              oMap               = new HashMap();
	private String[]             aoKeys;
	private String[]             aoArrayType        = new String[1];
	private boolean              tClearOfWarning    = true;

	private static final int STARTUP                = 0;
	private static final int DONE                   = 1;
	private static final int CAPTURING_HIT_SUMMARY  = 2;
	private static final int IN_HSP_COLLECTION      = 3;
	private static final int ON_FIRST_HSP           = 4;
	private static final int IN_HSP_SUMMARY         = 5;
	private static final int IN_ALIGNMENT           = 6;


	HitSectionSAXParser(BlastLikeVersionSupport poVersion,
			String poNamespacePrefix) {
		oVersion = poVersion;
		this.setNamespacePrefix(poNamespacePrefix);
		//For XSLT Parser Compliance
		this.addPrefixMapping("biojava","http://www.biojava.org");

		this.changeState(STARTUP);
		aoLineSeparator = System.getProperty("line.separator").toCharArray();
	}

	public void parse(BufferedReader poContents, String poLine, ArrayList<String> poEndSignals) throws SAXException {
		oLine = null;
		oContents = poContents;
		setGlobalEndSignal(poEndSignals);
		//return immediately if this is not the start
		//of a hit...
		if (!poLine.startsWith(">")) return;

		try {

			oLine = poLine;
			while ((oLine != null) &&
					(!this.matchesGlobalEndSignal(oLine)) &&
					(!(iState == DONE)) )
			{
				//interpret line and send messages accordingly
				this.interpret(oLine);
				//check for End again cos stream read elsewhere

				if (this.matchesGlobalEndSignal(oLine)) {
					this.changeState(DONE);
					oContents.reset();
					break;
				}

				oContents.mark(10000000);
				oLine = oContents.readLine();

			} // end while

		} catch (java.io.IOException x) {
			System.out.println(x.getMessage());
			System.out.println("File read interupted");
		} // end try/catch
	}
	/**
	 * Typically parse a single line, and return
	 *
	 * @param poLine     -
	 * @exception SAXException thrown if
	 * @exception  thrown if
	 */
	private void interpret(String poLine) throws SAXException {


		if (poLine.startsWith(">")) {
			//start of hit, accumulate title into buffer.
			//omit intial ">" character
			oStringBuffer.setLength(0);

			oStringBuffer.append(poLine.substring(1));
			this.changeState(CAPTURING_HIT_SUMMARY);
			return;
		}

		if (iState == CAPTURING_HIT_SUMMARY) {

			if (poLine.trim().startsWith("Length =")) {
				//here when HitSummary is complete

				//get sequenceLength, and then startElement
				StringTokenizer oSt = new StringTokenizer(poLine);

				//zip through tokens up to last one
				int iTmpTokenCount = oSt.countTokens() - 1;
				for (int i = 0; i < iTmpTokenCount; i++) {
					oSt.nextToken();
				}
				//last token is the length
				String oLength = oSt.nextToken();

				oAtts.clear();
				oAttQName.setQName("sequenceLength");
				oAtts.addAttribute(oAttQName.getURI(),
						oAttQName.getLocalName(),
						oAttQName.getQName(),
						"CDATA",oLength);

				this.startElement(new QName(this,this.prefix("Hit")),
						(Attributes)oAtts);

				//Here, oStringBuffer contains ID + Description
				oSt = new StringTokenizer(oStringBuffer.substring(0));

				int iCount = oSt.countTokens();

				String oId = oSt.nextToken(); //get Id

				oAtts.clear();
				oAttQName.setQName("id");
				oAtts.addAttribute(oAttQName.getURI(),
						oAttQName.getLocalName(),
						oAttQName.getQName(),
						"CDATA",oId);

				oAttQName.setQName("metaData");
				oAtts.addAttribute(oAttQName.getURI(),
						oAttQName.getLocalName(),
						oAttQName.getQName(),
						"CDATA","none");

				this.startElement(new QName(this,this.prefix("HitId")),
						(Attributes)oAtts);

				this.endElement(new QName(this,this.prefix("HitId")));

				oDescription.setLength(0);

				if (iCount > 0) {
					//deal with hit description if one available

					while (oSt.hasMoreTokens()) {
						//construct description
						oDescription.append(oSt.nextToken() + " ");
						//System.out.println(oDescription);
					}
					oAtts.clear();
					this.startElement(new QName(this,this.prefix("HitDescription")),
							(Attributes)oAtts);

					aoChars = oDescription.substring(0).trim().toCharArray();
					this.characters(aoChars,0,aoChars.length);

					this.endElement(new QName(this,this.prefix("HitDescription")));

				} //end if there is a hit description

				//Here when we have HitId and HitDescription

				this.changeState(IN_HSP_COLLECTION);

				return;
			} else {
				//here if collating multi-line hit descriptions
				oStringBuffer.append(" " + poLine.trim());
				return;
			}
		} //end capturingHitSummary

		//parse HSPs
		if (iState == IN_HSP_COLLECTION) {
			//Look for start of a new HSP
			if (poLine.trim().startsWith("Score")) {
				//here if on a new HSP
				oAtts.clear();
				this.startElement(new QName(this,this.prefix("HSPCollection")),
						(Attributes)oAtts);

				//Note, this method will have changed
				//the State when it returns
				this.firstHSPEvent(poLine);
				this.endElement(new QName(this,this.prefix("HSPCollection")));
				this.endElement(new QName(this,this.prefix("Hit")));
				this.changeState(CAPTURING_HIT_SUMMARY);
			}
		}
	}

	/**
	 * Deal with parsing of all HSPs in a Hit.
	 * Continue until a new Hit is reached...
	 *
	 * @param poLine     The first line of the HSP
	 *
	 */
	private void firstHSPEvent(String poLine) throws SAXException {

		this.changeState(ON_FIRST_HSP);


		try {

			oLine = poLine;
			while ((oLine != null) &&
					(!oLine.trim().startsWith(">")) &&
					(!this.matchesGlobalEndSignal(oLine)) )

			{
				//interpret line and send messages accordingly
				this.interpretHSP(oLine);
				oLine = oContents.readLine();
			} // end while

			//output final HSP of collection
			if (!(iState == ON_FIRST_HSP)) {
				//output previous HSP-related data
				this.outputHSPInfo();
				this.endElement(new QName(this,this.prefix("HSP")));
			}

		} catch (java.io.IOException x) {
			System.out.println(x.getMessage());
			System.out.println("File read interupted");
		} // end try/catch

		//Here at the end of dealing with HSPCollection
		//Could go on to next hit, or on to trailer...

		if(oLine==null) {
			return;
		}
		
		if (oLine.startsWith(">")) {
			//here when a new Hit is starting...


			//start of new hit, accumulate title into buffer.
			//omit intial ">" character

			oStringBuffer.setLength(0);

			oStringBuffer.append(oLine.substring(1));

			return;
		}

		if (this.matchesGlobalEndSignal(oLine)) {
			//here when we've hit the trailer...

			//this.endElement(this.prefix("HSP"));
			this.changeState(DONE);
			return;
		}

	}
	/**
	 * Deal with a line of an HSP
	 *
	 * @param poLine     A String representation of the line
	 */
	private void interpretHSP(String poLine) throws SAXException {


		//System.out.println("HSPLine:>".concat(poLine));
		//System.out.println("GlobalState:>".concat(iState));

		if (!tClearOfWarning) {
			//look for white space to indicate we're passed a multi-line
			//warning (in WU-BLAST);
			if (poLine.trim().equals("")) {
				//here when clear
				tClearOfWarning = true;
			}
			return;
		}


		//ignore Minus Strand HSP and Plus Strand HSP (WuBlast)

		if (poLine.trim().toLowerCase().startsWith("minus strand")) return;
		if (poLine.trim().toLowerCase().startsWith("plus strand")) return;

		if (poLine.trim().toLowerCase().startsWith("warning")) {
			tClearOfWarning = false;
			return;
		}



		if (poLine.trim().startsWith("Score")) {
			if (!(iState == ON_FIRST_HSP)) {
				//output previous HSP-related data
				this.outputHSPInfo();
				this.endElement(new QName(this,this.prefix("HSP")));

			}
			oAtts.clear();
			this.startElement(new QName(this,this.prefix("HSP")),
					(Attributes)oAtts);

			//Start accumulating all HSP summary information
			//into buffer...

			oStringBuffer.setLength(0);
			oStringBuffer.append(poLine);

			//and raw info

			oBuffer.clear();
			oBuffer.add(poLine);
			this.changeState(IN_HSP_SUMMARY);

			return;
		}
		//continue to accumulate summary info
		//until an alignment is reached...
		if (iState == IN_HSP_SUMMARY) {
			//check for end of summary (Query: is end signal)
			if (poLine.startsWith("Query:")) {

				//System.out.println(oStringBuffer);

				//at this point, all available summary info
				//complete for current HSP (may need
				//extra info derived from alignment, so
				//so don't output HSPSummary element info yet

				//Put available HSPSummary info into a Map
				HSPSummaryHelper.parse(oStringBuffer.substring(0),oMap,
						oVersion);

				//really need to get alignment parsed before outputing
				//suummary info - not all programs output
				//alignment size (e.g. DBA).

				//change state 'cos hit a Blast-like alignment
				this.changeState(IN_ALIGNMENT);

				//get information for first alignment line

				oAlignmentBuffer.clear();
				oAlignmentBuffer.add(poLine);

				return;
			} //end if found first line of alignment
			//append summary

			//ignore blank lines
			if (poLine.trim().equals("")) return;



			oBuffer.add(poLine); //keep raw info

			oStringBuffer.append(", ");
			oStringBuffer.append(poLine);

			return;
		} //end if state is inHSPSummary


		//keep appending alignment info
		if (iState == IN_ALIGNMENT) {
			//ignore blank lines
			if (poLine.trim().equals("")) return;
			oAlignmentBuffer.add(poLine);
			return;
		}

	}


	/**
	 * Describe 'outputHSPInfo' method here.
	 *
	 * @param nil    -
	 * @exception SAXException thrown if
	 * @exception  thrown if
	 */
	private void outputHSPInfo() throws SAXException {
		//Output HSP Summary info

		//detailed info
		aoKeys = (String[])(oMap.keySet().toArray(aoArrayType));

		oAtts.clear();

		for (int i = 0; i < aoKeys.length; i++) {
			if ( (aoKeys[i].equals("queryFrame"))  ||
					(aoKeys[i].equals("hitFrame"))    ||
					(aoKeys[i].equals("queryStrand")) ||
					(aoKeys[i].equals("hitStrand")) ) {
				//nametoken if an enumeration


				oAttQName.setQName(aoKeys[i]);
				oAtts.addAttribute(oAttQName.getURI(),
						oAttQName.getLocalName(),
						oAttQName.getQName(),
						"NMTOKEN",(String)oMap.get(aoKeys[i]));
			} else {
				//CDATA if regular attribute

				oAttQName.setQName(aoKeys[i]);
				oAtts.addAttribute(oAttQName.getURI(),
						oAttQName.getLocalName(),
						oAttQName.getQName(),
						"CDATA",(String)oMap.get(aoKeys[i]));
			}
			//System.out.print(aoKeys[i] + ": ");
			//System.out.println(oMap.get(aoKeys[i]));
		}


		this.startElement(new QName(this,this.prefix("HSPSummary")),
				(Attributes)oAtts);
		//Raw HSPSummary Data
		oAtts.clear();
		oAttQName.setQName("xml:space");
		oAtts.addAttribute(oAttQName.getURI(),
				oAttQName.getLocalName(),
				oAttQName.getQName(),
				"NMTOKEN","preserve");
		this.startElement(new QName(this,this.prefix("RawOutput")),
				(Attributes)oAtts);

		int iTmpBufferSize = oBuffer.size();
		for (int i = 0; i < iTmpBufferSize;i++) {
			//aoChars = ((String)oBuffer.get(i)).trim().toCharArray();
			aoChars = ((String)oBuffer.get(i)).toCharArray();
			this.characters(aoChars,0,aoChars.length);
			this.characters(aoLineSeparator,0,1);

		}
		this.endElement(new QName(this,this.prefix("RawOutput")));

		this.endElement(new QName(this,this.prefix("HSPSummary")));

		//Output Alignment info via delegation to
		//a BlastLikeAlignmentSAXParser

		oAlignmentParser =
			new BlastLikeAlignmentSAXParser(this.getNamespacePrefix());

		oAlignmentParser.setContentHandler(oHandler);

		oAlignmentParser.parse(oAlignmentBuffer);
	}

	private void setGlobalEndSignal(ArrayList<String> oGlobalEndSignal) {
		this.oGlobalEndSignals = oGlobalEndSignal;
	}


	private boolean matchesGlobalEndSignal(String s) {
		if(s==null) {
			return true;
		}
		for(String signal: this.oGlobalEndSignals) {
			if(s.trim().startsWith(signal)) {
				return true;
			}
		}
		return false;
	}
}
