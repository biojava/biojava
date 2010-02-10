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

package org.biojava.bio.program.gff;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.utils.ParserException;
import org.biojava.utils.SmallMap;

/**
 * Parse a stream of GFF text into a stream of records and comments.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Keith James (docs)
 */
public class GFFParser {
    private GFFErrorHandler errors = GFFErrorHandler.ABORT_PARSING;

    /**
     * Set the error handler used by this parser.
     */

    public void setErrorHandler(GFFErrorHandler errors) {
	this.errors = errors;
    }

    /**
     * Find the error handler used by this parser.
     */

    public GFFErrorHandler getErrorHandler() {
	return errors;
    }

    /**
     * Informs <span class="arg">handler</span> of each line of
     * gff read from <span class="arg">bReader</span>.  This form
     * of the method should only be used if no locator string is
     * available for the resource being parsed.
     *
     * @param bReader the <span class="type">BufferedReader</span> to parse
     * @param handler the <span class="type">GFFDocumentHandler</span> that will
     *                listen for 'stuff'
     *
     * @throws <span class="type">IOException</span> if for any reason
     *         <span class="arg">bReader</span> throws one
     * @throws <span class="type">BioException</span> if
     *         <span class="arg">handler</span> can not correct a parse error
     */

    public void parse(BufferedReader bReader, GFFDocumentHandler handler)
	throws IOException, BioException, ParserException
    {
	parse(bReader, handler, "unknown:");
    }

    /**
     * Informs <span class="arg">handler</span> of each line of
     * GFF read from <span class="arg">bReader</span>
     *
     * @param bReader the <span class="type">BufferedReader</span> to parse
     * @param handler the <span class="type">GFFDocumentHandler</span> that will
     *                listen for 'stuff'
     *
     * @throws <span class="type">IOException</span> if for any reason
     *         <span class="arg">bReader</span> throws one
     * @throws <span class="type">BioException</span> if
     *         <span class="arg">handler</span> can not correct a parse error
     */

    public void parse(BufferedReader bReader, GFFDocumentHandler handler, String locator)
	throws IOException, BioException, ParserException
    {
	handler.startDocument(locator);
	ArrayList aList = new ArrayList();
	int lineNum = 0;
	for(String line = bReader.readLine(); line != null; line = bReader.readLine()) {
	    ++lineNum;

	    try {
		aList.clear();
		if(line.startsWith("#")) {
		    handler.commentLine(line.substring(1));
		} else if (line.length() == 0) {
		} else {
		    StringTokenizer st = new StringTokenizer(line, "\t", false);
		    while(st.hasMoreTokens() && aList.size() < 8) {
			String token = st.nextToken();
			aList.add(token);
		    }

                    if(aList.size() < 7) {
                      throw new ParserException(
                        "Line doesn't look like GFF",
                        locator,
                        lineNum,
                        line );
                    }

		    String rest = null;
		    String comment = null;
		    if(st.hasMoreTokens()) {
			try {
			    rest = st.nextToken(((char) 0) + "");
			} catch (NoSuchElementException nsee) {
			}
		    }
		    if(rest != null) {
			int ci = rest.indexOf("#");
			if (ci != -1) {
			    comment = rest.substring(ci);
			    rest = rest.substring(0, ci);
			}
		    }

		    GFFRecord record = createRecord(handler, aList, rest, comment);
		    handler.recordLine(record);
		}
	    } catch (ParserException ex) {
		throw new ParserException(ex.getMessage(),
					  locator,
					  lineNum,
					  line);
	    } catch (IgnoreRecordException ex) {
		// Silently skip any more work on this record
	    }
	}
	handler.endDocument();
    }

  /**
   * Actually turns a list of tokens, some value string and a comment into a
   * <span class="type">GFFRecord</span> and informs
   * <span class="arg">handler</span>.
   *
   * @param handler a <span class="type">GFFDocumentHandler</span> to inform of
   *                any parse errors, and the completed <span class="type">GFFRecord</span>
   * @param aList   a <span class="type">List</span> containing the 8 mandatory GFF columns
   * @param rest    a <span class="type">String</span> representing the unparsed
   *                attribute-value text, or <span class="kw">null</span> if there is none
   * @param comment a <span class="type">String</span> containing the comment (without the
   *                leading '<code>#</code>' character.
   * @throws <span class="type">BioException</span> if <span class="arg">handler</span>
   *         could not correct a parse error
   */
    protected GFFRecord createRecord(GFFDocumentHandler handler,
				     List aList,
				     String rest,
				     String comment)
	throws BioException, ParserException, IgnoreRecordException
    {
	SimpleGFFRecord record = new SimpleGFFRecord();

	record.setSeqName((String) aList.get(0));
	record.setSource((String) aList.get(1));
	record.setFeature((String) aList.get(2));

	int start = -1;
	try {
	    start = Integer.parseInt( (String) aList.get(3));
	} catch (NumberFormatException nfe) {
	    start = errors.invalidStart((String) aList.get(3));
	}
	record.setStart(start);

	int end = -1;
	try {
	    end = Integer.parseInt( (String) aList.get(4));
	} catch (NumberFormatException nfe) {
	    end = errors.invalidEnd((String) aList.get(3));
	}
	record.setEnd(end);

	String score = (String) aList.get(5);
	if(score == null     ||
	   score.equals("")  ||
	   score.equals(".") ||
	   score.equals("0")
	   ) 
	{
	    record.setScore(GFFTools.NO_SCORE);
	} else {
	    double sc = 0.0;
	    try {
		sc = Double.parseDouble(score);
	    } catch (NumberFormatException nfe) {
		sc = errors.invalidScore(score);
	    }
	    record.setScore(sc);
	}

	String strand = (String) aList.get(6);
	if(strand == null || strand.equals("") || strand.equals(".")) {
	    record.setStrand(StrandedFeature.UNKNOWN);
	} else {
	    if(strand.equals("+")) {
		record.setStrand(StrandedFeature.POSITIVE);
	    } else if(strand.equals("-")) {
		record.setStrand(StrandedFeature.NEGATIVE);
	    } else {
		record.setStrand(errors.invalidStrand(strand));
	    }
	}

	String frame = (String) aList.get(7);
	if(frame.equals(".")) {
	    record.setFrame(GFFTools.NO_FRAME);
	} else {
	    int fr = 0;
	    try {
		fr = Integer.parseInt(frame);
	    } catch (NumberFormatException nfe) {
		fr = errors.invalidFrame(frame);
	    }
	    record.setFrame(fr);
	}

	if (rest != null)
	    record.setGroupAttributes(parseAttribute(rest));
	else
	    record.setGroupAttributes(new SmallMap());
	record.setComment(comment);

	return record;
    }

    /**
     * Parse <span class="arg">attValList</span> into a
     * <span class="type">Map</span> of attributes and value lists.
     * <p>
     * The resulting <span class="type">Map</span> will have
     * <span class="type">String</span> keys, with
     * <span class="type">List</span> values. If there are no values
     * associated with a key, then it will have an empty
     * <span class="type">List</span>, not <span class="kw">null</span> as
     * its value.
     *
     * @param attValList  the <span class="type">String</span> to parse
     * @return a <span class="type">Map</span> of parsed attributes and value lists
     */

    protected Map parseAttribute(String attValList) {
	Map attMap = new SmallMap();

	StringTokenizer sTok = new StringTokenizer(attValList, ";", false);
	while(sTok.hasMoreTokens()) {
	    String attVal = sTok.nextToken().trim();
	    String attName;
	    List valList = new ArrayList();
	    int spaceIndx = attVal.indexOf(" ");
	    if(spaceIndx == -1) {
		attName = attVal;
	    } else {
		attName = attVal.substring(0, spaceIndx);
		attValList = attVal.substring(spaceIndx).trim();
		while(attValList.length() > 0) {
		    if(attValList.startsWith("\"")) {
			// System.out.println("Quoted");
			int quoteIndx = 0;
			do {
			    quoteIndx++;
			    quoteIndx = attValList.indexOf("\"", quoteIndx);
			} while(quoteIndx != -1 && attValList.charAt(quoteIndx-1) == '\\');
			if(quoteIndx > 0){
                          valList.add(attValList.substring(1, quoteIndx));
			  attValList = attValList.substring(quoteIndx+1).trim();
			}else{
                          valList.add(attValList);
                          attValList = "";
			}
		    } else {
			spaceIndx = attValList.indexOf(" ");
			if(spaceIndx == -1) {
			    valList.add(attValList);
			    attValList = "";
			} else {
			    valList.add(attValList.substring(0, spaceIndx));
			    attValList = attValList.substring(spaceIndx).trim();
			}
		    }
		}
	    }
	    attMap.put(attName, valList);
	}

	return attMap;
    }
}
