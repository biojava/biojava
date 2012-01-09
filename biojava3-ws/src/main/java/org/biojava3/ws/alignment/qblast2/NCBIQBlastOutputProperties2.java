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
 * Created on 2011-11-20
 *
 */
package org.biojava3.ws.alignment.qblast2;

import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.ALIGNMENTS;
import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.ALIGNMENT_VIEW;
import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.DESCRIPTIONS;
import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.FORMAT_TYPE;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.biojava3.ws.alignment.RemotePairwiseAlignmentOutputProperties;
import org.biojava3.ws.alignment.qblast2.enums.BlastOutputAlignmentFormat;
import org.biojava3.ws.alignment.qblast2.enums.BlastOutputFormat2;

/**
 * This class implements {@code RemotePairwiseAlignmentOutputProperties} by
 * adding several convenient methods used to wrap the addition of Blast output
 * parameters.
 * <p>
 * The constructor for this class assigns default output format values. To set
 * the options either use {@link #setOutputOption(String, String)} method or the
 * wrapper methods.
 * 
 * @author Sylvain Foisy, Diploide BioIT
 * @author Gediminas Rimsa
 */
public class NCBIQBlastOutputProperties2 implements RemotePairwiseAlignmentOutputProperties {
	private static final long serialVersionUID = -3940778840861986573L;

	private Map<String, String> param = new HashMap<String, String>();

	/**
	 * This constructor builds the parameters for the output of the GET command
	 * sent to the QBlast service with default values:
	 * 
	 * <pre>
	 * FORMAT_TYPE = XML;
	 * ALIGNMENT_VIEW = Pairwise;
	 * DESCRIPTIONS = 100;
	 * ALIGNMENTS = 100;
	 * </pre>
	 */
	public void NCBIQBlastOutputProperties() {
		param.put(FORMAT_TYPE, BlastOutputFormat2.XML.toString());
		param.put(ALIGNMENT_VIEW, BlastOutputAlignmentFormat.PAIRWISE.toString());
		param.put(DESCRIPTIONS, "100");
		param.put(ALIGNMENTS, "100");
	}

	/**
	 * @return stream output format - a String with the value of key FORMAT_TYPE
	 */
	public String getOutputFormat() {
		return param.get(FORMAT_TYPE);
	}

	/**
	 * Sets the stream output format to get from the QBlast service
	 * <p>
	 * If {@code HTML} format is selected, also adds parameters (which are
	 * removed if another output format is chosen):
	 * 
	 * <pre>
	 * NOHEADER = true;
	 * SHOW_OVERVIEW = false;
	 * SHOW_LINKOUT = false;
	 * </pre>
	 * 
	 * @param formatType : one of the output format types defined in enum
	 */
	public void setOutputFormat(BlastOutputFormat2 formatType) {
		param.put(FORMAT_TYPE, formatType.getValue());
		if (BlastOutputFormat2.HTML.equals(formatType)) {
			// add default parameters associated with HTML
			param.put("NOHEADER", "true");
			param.put("SHOW_OVERVIEW", "false");
			param.put("SHOW_LINKOUT", "false");
		} else {
			// remove default parameters associated with HTML
			param.remove("NOHEADER");
			param.remove("SHOW_OVERVIEW");
			param.remove("SHOW_LINKOUT");
		}
	}

	/**
	 * @return alignment output format - a String with the value of key
	 *         ALIGNMENT_VIEW
	 */
	public String getAlignmentOutputFormat() {
		return param.get(ALIGNMENT_VIEW);
	}

	/**
	 * Sets the alignment output format to get from the QBlast service
	 * 
	 * @param alignmentFormat : one of available alignment types
	 */
	public void setAlignmentOutputFormat(BlastOutputAlignmentFormat alignmentFormat) throws Exception {
		param.put(ALIGNMENT_VIEW, alignmentFormat.getValue());
	}

	/**
	 * @return number of descriptions fetched - an int with the value of the key
	 *         DESCRIPTIONS
	 */
	public int getDescriptionNumber() {
		return Integer.parseInt(param.get(DESCRIPTIONS));
	}

	/**
	 * Sets the number of descriptions to fetch
	 * 
	 * @param number : an int with the required number of descriptions to fetch
	 */
	public void setDescriptionNumber(int number) {
		param.put(DESCRIPTIONS, Integer.toString(number));
	}

	/**
	 * @return number of alignments fetched - an int with the value of the key
	 *         ALIGNMENTS
	 */
	public int getAlignmentNumber() {
		return Integer.parseInt(param.get(ALIGNMENTS));
	}

	/**
	 * Set the number of alignments to fetch
	 * 
	 * @param number : an int with the required number of alignments to fetch
	 */
	public void setAlignmentNumber(int number) {
		param.put(ALIGNMENTS, Integer.toString(number));
	}

	/**
	 * @throws IllegalArgumentException if output options map does not contain
	 *             given key
	 */
	@Override
	public String getOutputOption(String key) {
		if (param.containsKey(key)) {
			return param.get(key);
		} else {
			throw new IllegalArgumentException("The key named " + key + " is not set in this RemoteQBlastOutputProperties object");
		}
	}

	@Override
	public void setOutputOption(String key, String value) {
		param.put(key, value);
	}

	@Override
	public Set<String> getOutputOptions() {
		return param.keySet();
	}
}