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
package org.biojava3.ws.alignment.qblast;

import static org.biojava3.ws.alignment.qblast.BlastOutputParameterEnum.ALIGNMENTS;
import static org.biojava3.ws.alignment.qblast.BlastOutputParameterEnum.ALIGNMENT_VIEW;
import static org.biojava3.ws.alignment.qblast.BlastOutputParameterEnum.DESCRIPTIONS;
import static org.biojava3.ws.alignment.qblast.BlastOutputParameterEnum.FORMAT_TYPE;
import static org.biojava3.ws.alignment.qblast.BlastOutputParameterEnum.NOHEADER;
import static org.biojava3.ws.alignment.qblast.BlastOutputParameterEnum.SHOW_LINKOUT;
import static org.biojava3.ws.alignment.qblast.BlastOutputParameterEnum.SHOW_OVERVIEW;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.biojava3.ws.alignment.RemotePairwiseAlignmentOutputProperties;

/**
 * This class wraps a QBlast output parameter {@code Map} by adding several convenient parameter addition methods. Other
 * QBlast URL API parameters should be added using {@link #setOutputOption(BlastOutputParameterEnum, String)}
 * 
 * @author Sylvain Foisy, Diploide BioIT
 * @author Gediminas Rimsa
 */
public class NCBIQBlastOutputProperties implements RemotePairwiseAlignmentOutputProperties {
	private static final long serialVersionUID = -9202060390925345163L;

	private Map<BlastOutputParameterEnum, String> param = new HashMap<BlastOutputParameterEnum, String>();

	/**
	 * This constructor builds the parameters for the output of the GET command sent to the QBlast service with default
	 * values:
	 * 
	 * <pre>
	 * FORMAT_TYPE = XML;
	 * ALIGNMENT_VIEW = Pairwise;
	 * DESCRIPTIONS = 100;
	 * ALIGNMENTS = 100;
	 * </pre>
	 */
	public NCBIQBlastOutputProperties() {
		setOutputOption(FORMAT_TYPE, BlastOutputFormatEnum.XML.name());
		setOutputOption(ALIGNMENT_VIEW, BlastOutputAlignmentFormatEnum.Pairwise.name());
		setOutputOption(DESCRIPTIONS, "100");
		setOutputOption(ALIGNMENTS, "100");
	}

	/**
	 * This method forwards to {@link #getOutputOption(BlastOutputParameterEnum)}. Consider using it instead.
	 */
	@Override
	public String getOutputOption(String key) {
		return getOutputOption(BlastOutputParameterEnum.valueOf(key));
	}

	/**
	 * This method forwards to {@link #setOutputOption(BlastOutputParameterEnum, String)}. Consider using it instead.
	 */
	@Override
	public void setOutputOption(String key, String val) {
		setOutputOption(BlastOutputParameterEnum.valueOf(key), val);
	}

	/**
	 * Gets the value of specified parameter or {@code null} if it is not set
	 */
	public String getOutputOption(BlastOutputParameterEnum key) {
		return param.get(key);
	}

	/**
	 * Sets the value of specified output parameter
	 */
	public void setOutputOption(BlastOutputParameterEnum key, String value) {
		param.put(key, value);
	}

	/**
	 * Gets output parameters, which are currently set
	 */
	@Override
	public Set<String> getOutputOptions() {
		Set<String> result = new HashSet<String>();
		for (BlastOutputParameterEnum parameter : param.keySet()) {
			result.add(parameter.name());
		}
		return result;
	}

	/**
	 * Removes given parameter
	 */
	public void removeOutputOption(BlastOutputParameterEnum key) {
		param.remove(key);
	}

	/**
	 * @return stream output format - a String with the value of key FORMAT_TYPE
	 */
	public String getOutputFormat() {
		return getOutputOption(FORMAT_TYPE);
	}

	/**
	 * Sets the stream output format to get from the QBlast service
	 * <p/>
	 * If {@code HTML} format is selected, also adds the following parameters (which are removed if another output
	 * format is chosen):
	 * 
	 * <pre>
	 * NOHEADER = true;
	 * SHOW_OVERVIEW = false;
	 * SHOW_LINKOUT = false;
	 * </pre>
	 * 
	 * @param formatType : one of the output format types defined in enum
	 */
	public void setOutputFormat(BlastOutputFormatEnum formatType) {
		setOutputOption(FORMAT_TYPE, formatType.name());
		if (BlastOutputFormatEnum.HTML.equals(formatType)) {
			// add default parameters associated with HTML
			setOutputOption(NOHEADER, "true");
			setOutputOption(SHOW_OVERVIEW, "false");
			setOutputOption(SHOW_LINKOUT, "false");
		} else {
			// remove default parameters associated with HTML
			removeOutputOption(NOHEADER);
			removeOutputOption(SHOW_OVERVIEW);
			removeOutputOption(SHOW_LINKOUT);
		}
	}

	/**
	 * @return alignment output format - a String with the value of key ALIGNMENT_VIEW
	 */
	public String getAlignmentOutputFormat() {
		return getOutputOption(ALIGNMENT_VIEW);
	}

	/**
	 * Sets the alignment output format to get from the QBlast service
	 * 
	 * @param alignmentFormat : one of available alignment types
	 */
	public void setAlignmentOutputFormat(BlastOutputAlignmentFormatEnum alignmentFormat) {
		setOutputOption(ALIGNMENT_VIEW, alignmentFormat.name());
	}

	/**
	 * @return number of descriptions fetched - an int with the value of the key DESCRIPTIONS
	 */
	public int getDescriptionNumber() {
		return Integer.parseInt(getOutputOption(DESCRIPTIONS));
	}

	/**
	 * Sets the number of descriptions to fetch
	 * 
	 * @param number : an int with the required number of descriptions to fetch
	 */
	public void setDescriptionNumber(int number) {
		setOutputOption(DESCRIPTIONS, Integer.toString(number));
	}

	/**
	 * @return number of alignments fetched - an int with the value of the key ALIGNMENTS
	 */
	public int getAlignmentNumber() {
		return Integer.parseInt(getOutputOption(ALIGNMENTS));
	}

	/**
	 * Set the number of alignments to fetch
	 * 
	 * @param number : an int with the required number of alignments to fetch
	 */
	public void setAlignmentNumber(int number) {
		setOutputOption(ALIGNMENTS, Integer.toString(number));
	}

}