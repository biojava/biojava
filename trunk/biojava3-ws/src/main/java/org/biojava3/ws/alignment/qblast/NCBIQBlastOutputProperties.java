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

package org.biojava3.ws.alignment.qblast;

import java.util.HashMap;
import java.util.Set;

import org.biojava3.ws.alignment.RemotePairwiseAlignmentOutputProperties;

/**
 * The actual implementation of the RemotePairwiseAlignmentOutputProperties 
 * interface for the QBlast service.
 * 
 * The constructor for this class builds an object with default format values. Any modification will 
 * either use the generic setOutputOption method or use the wrapper methods that are actually 
 * build around the generic method.

 * @author Sylvain Foisy, Diploide BioIT
 * @since Biojava 3
 *
 */
public class NCBIQBlastOutputProperties implements
		RemotePairwiseAlignmentOutputProperties {

	private static final long serialVersionUID = 1L;
	private HashMap<String, String> out = new HashMap<String, String>();
	private String outFormat = "FORMAT_TYPE=Text";
	private String alignFormat = "ALIGNMENT_VIEW=Pairwise";

	private int descNumbers = 100;
	private int alignNumbers = 100;

	/**
	 * This constructor build the parameters for the default output of the GET command sent to the QBlast service.
	 * Here are the default values:
	 * 
	 * FORMAT_TYPE = Text;
	 * ALIGNMENT_VIEW = Pairwise
	 * DESCRIPTIONS = 100;
	 * ALIGNMENTS = 100;
	 * 
	 */
	public void NCBIQBlastOutputProperties() {
		this.out.put("FORMAT_TYPE", outFormat);
		this.out.put("ALIGNMENT_VIEW", alignFormat);
		this.out.put("DESCRIPTIONS", "DESCRIPTIONS=" + descNumbers);
		this.out.put("ALIGNMENTS", "ALIGNMENTS=" + alignNumbers);
	}

	/**
	 * Simply returns stream output format for the actual RemoteQBlastOutputProperties object
	 * 
	 * @return a String with the value of key FORMAT_TYPE.
	 */
	public String getOutputFormat() {
		return this.out.get("FORMAT_TYPE");
	}

	/**
	 * This method is use to set the stream output format to get from the QBlast service
	 * 
	 * @param rf :an enum from RemoteQBlastOutputFormat
	 * @throws Exception if the enum is neither of RemoteQBlastOutputFormat.TEXT/XML/HTML
	 */
	public void setOutputFormat(NCBIQBlastOutputFormat rf)
			throws Exception {
		switch (rf) {
		case TEXT:
			this.outFormat = "FORMAT_TYPE=Text";
			this.out.put("FORMAT_TYPE", outFormat);
			break;
		case XML:
			this.outFormat = "FORMAT_TYPE=XML";
			this.out.put("FORMAT_TYPE", outFormat);
			break;
		case HTML:
			this.outFormat = "FORMAT_TYPE=HTML";
			this.out.put("FORMAT_TYPE", outFormat);
			break;
		default:
			throw new Exception(
					"Unacceptable selection of format type. Only values text / XML / HTML accepted");
		}
	}

	/**
	 * Method that returns the alignment output format for this actual RemoteQBlastOutputProperties object
	 * 
	 * @return a String with the value of key ALIGNMENT_VIEW
	 */
	public String getAlignmentOutputFormat() {
		return this.out.get("ALIGNMENT_VIEW");
	}

	/**
	 * This method is use to set the alignment output format to get from the QBlast service
	 * 
	 * @param rf :an enum from RemoteQBlastOutputFormat
	 * @throws Exception if the enum is neither of RemoteQBlastOutputFormat.PAIRWISE/QUERY_ANCHORED/QUERY_ANCHORED_NO_IDENTITIES/FLAT_QUERY_ANCHORED
	 *         FLAT_QUERY_ANCHORED_NO_IDENTITIES/TABULAR
	 */
	public void setAlignmentOutputFormat(NCBIQBlastOutputFormat rf)
			throws Exception {
		switch (rf) {
		case PAIRWISE:
			this.alignFormat = "ALIGNMENT_VIEW=Pairwise";
			this.out.put("ALIGNMENT_VIEW", alignFormat);
			break;
		case QUERY_ANCHORED:
			this.alignFormat = "ALIGNMENT_VIEW=QueryAnchored";
			this.out.put("ALIGNMENT_VIEW", alignFormat);
			break;
		case QUERY_ANCHORED_NO_IDENTITIES:
			this.alignFormat = "ALIGNMENT_VIEW=QueryAnchoredNoIdentities";
			this.out.put("ALIGNMENT_VIEW", alignFormat);
			break;
		case FLAT_QUERY_ANCHORED:
			this.alignFormat = "ALIGNMENT_VIEW=FlatQueryAnchored";
			this.out.put("ALIGNMENT_VIEW", alignFormat);
			break;
		case FLAT_QUERY_ANCHORED_NO_IDENTITIES:
			this.alignFormat = "ALIGNMENT_VIEW=FlatQueryAnchoredNoIdentities";
			this.out.put("ALIGNMENT_VIEW", alignFormat);
			break;
		case TABULAR:
			this.alignFormat = "ALIGNMENT_VIEW=Tabular";
			this.out.put("ALIGNMENT_VIEW", alignFormat);
			break;
		default:
			throw new Exception(
					"Unacceptable selection of alignment type. Only values Pairwise / QueryAnchored / QueryAnchoredNoIdentities / FlatQueryAnchored / FlatQueryAnchoredNoIdentities / Tabular accepted");
		}
	}

	/**
	 * A method that simply returns the number of descriptions fetched with this RemoteQBlastOutputProperties object.
	 * 
	 * @return an int with the value of the key DESCRIPTIONS
	 */
	public int getDescriptionNumber() {
		String val = this.out.get("DESCRIPTIONS");
		String vals[] = val.split("=");
		int i = Integer.parseInt(vals[1]);

		return i;
	}

	/**
	 * A method to set the number of descriptions to fetch to the GET command.
	 * 
	 * @param i :an int with the required number of descriptions to fetch.
	 */
	public void setDescriptionNumber(int i) {
		this.descNumbers = i;
		this.out.put("DESCRIPTIONS", "DESCRIPTIONS=" + descNumbers);
	}

	/**
	 * A method that simply returns the number of alignments fetched with this RemoteQBlastOutputProperties object.
	 * 
	 * @return an int with the value of the key ALIGNMENTS.
	 */	
	public int getAlignmentNumber() {
		String val = this.out.get("ALIGNMENTS");
		String vals[] = val.split("=");
		int i = Integer.parseInt(vals[1]);

		return i;
	}

	/**
	 * A method to set the number of alignments to fetch to the GET command.
	 * 
	 * @param i :an int with the required number of alignments to fetch.
	 */

	public void setAlignmentNumber(int i) {
		this.alignNumbers = i;
		this.out.put("ALIGNMENTS", "ALIGNMENTS=" + alignNumbers);
	}

	/**
	 * Method that returns any value associated to any key for this RemoteQBlastOutputProperties object.
	 * 
	 */
	public String getOutputOption(String o) throws Exception{
		if(out.containsKey(o)){
			return this.out.get(o);}
		else{
			throw new Exception("The key named "+o+" is not set in this RemoteQBlastOutputProperties object");
		}
	}

	public void setOutputOption(String o, String v) {
		this.out.put(o, v);
	}

	public Set<String> getOutputOptions() {
		return out.keySet();
	}

}