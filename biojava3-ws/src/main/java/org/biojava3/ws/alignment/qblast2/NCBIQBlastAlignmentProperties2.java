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

import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.DATABASE;
import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.EXPECT;
import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.GAPCOSTS;
import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.MATRIX_NAME;
import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.MEGABLAST;
import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.OTHER_ADVANCED;
import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.PROGRAM;
import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.QUERY_FROM;
import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.QUERY_TO;
import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.WORD_SIZE;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.biojava3.ws.alignment.RemotePairwiseAlignmentProperties;
import org.biojava3.ws.alignment.qblast2.enums.BlastMatrix;
import org.biojava3.ws.alignment.qblast2.enums.BlastProgram;

/**
 * This class implements {@code RemotePairwiseAlignmentProperties} by adding several convenient methods used to wrap the
 * addition of Blast alignment parameters.
 * <p>
 * Required parameters are PROGRAM and DATABASE, other parameters are optional
 * <p>
 * Any other QBlast URL API parameters can be added using {@link #setAlignementOption(String, String)}
 * 
 * @author Sylvain Foisy, Diploide BioIT
 * @author Gediminas Rimsa
 */
public class NCBIQBlastAlignmentProperties2 implements RemotePairwiseAlignmentProperties {

	private static final long serialVersionUID = -7710581294480758153L;

	private Map<String, String> param = new HashMap<String, String>();

	/**
	 * @return the name of {@linkplain BlastProgram} used for blast run or value of PROGRAM parameter if it does not
	 *         match any of enumerated BlastProgram types
	 * @see #getBlastProgram()
	 */
	public String getBlastProgramName() {
		try {
			return getBlastProgram().getValue();
		} catch (IllegalArgumentException e) {
			return param.get(PROGRAM);
		}
	}

	/**
	 * @return {@linkplain BlastProgram} used for blast run
	 * @throws IllegalArgumentException if BlastProgram with given name does not exist
	 * @see #getBlastProgramName()
	 */
	public BlastProgram getBlastProgram() {
		BlastProgram program = BlastProgram.get(param.get(PROGRAM));
		if (BlastProgram.BLASTN == program && param.get(MEGABLAST).equals("on")) {
			return BlastProgram.MEGABLAST;
		}
		return program;
	}

	/**
	 * Sets the program to be used with blastall
	 * 
	 * @param program : one of blastall programs
	 */
	public void setBlastProgram(BlastProgram program) {
		if (BlastProgram.MEGABLAST != program) {
			param.put(PROGRAM, program.getValue());
			param.remove(MEGABLAST);
		} else {
			param.put(PROGRAM, BlastProgram.BLASTN.getValue());
			param.put(MEGABLAST, "on");
		}
	}

	/**
	 * @return name of database used with blastall
	 */
	public String getBlastDatabase() {
		return param.get(DATABASE);
	}

	/**
	 * Sets the database to be used with blastall
	 * <p>
	 * A quite exhaustive list of the databases available for QBlast requests can be found here:
	 * <p>
	 * http://&lt to_be_completed &gt
	 * <p>
	 * Blastall equivalent: -d
	 * 
	 * @param db : a valid name to a NCBI blastable database
	 */
	public void setBlastDatabase(String database) {
		param.put(DATABASE, database);
	}

	/**
	 * @return double value of EXPECT parameter used for blast run
	 */
	public double getBlastExpect() {
		if (param.containsKey(EXPECT)) {
			return Double.parseDouble(param.get(EXPECT));
		}
		return 10;
	}

	/**
	 * Sets the EXPECT parameter to be use with blastall
	 * <p>
	 * Example: if you want a EXPECT of 1e-10, pass {@code Double.parseDouble("1e-10")} as a parameter
	 * <p>
	 * Blastall equivalent: -e
	 * 
	 * @param expect: a double value of EXPECT parameter
	 */
	public void setBlastExpect(double expect) {
		param.put(EXPECT, Double.toString(expect));
	}

	/**
	 * Returns the value of the WORD_SIZE parameter used for this blast run
	 * 
	 * @return int value of WORD_SIZE used by this search
	 * @throws IllegalArgumentException when program type is not set and program type is not supported
	 */
	public int getBlastWordSize() {
		if (param.containsKey(WORD_SIZE)) {
			return Integer.parseInt(param.get(WORD_SIZE));
		}

		// return default word size value
		try {
			BlastProgram programType = getBlastProgram();
			switch (programType) {
			case BLASTN:
				return 11;
			case MEGABLAST:
				return 28;
			default: // blastp / blastx / tbastn / tblastx
				return 3;
			}
		} catch (IllegalArgumentException e) {
			throw new IllegalArgumentException("Blast program " + getBlastProgramName() + " is not supported.", e);
		}
	}

	/**
	 * Sets the WORD_SIZE parameter to be use with blastall
	 * <p>
	 * <b>WARNING!!</b> At this point, the method does not verify the validity of your choice; for example, word size of
	 * greater than 5 with blastp returns error messages from QBlast. Word size range depends on the algorithm chosen.
	 * <p>
	 * More at http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node74.html
	 * <p>
	 * Blastall equivalent: -W
	 * 
	 * @param word: an int used to set WORD_SIZE
	 */
	public void setBlastWordSize(int word) {
		param.put(WORD_SIZE, Integer.toString(word));
	}

	/**
	 * Returns the value for the GAP_CREATION parameter (first half of GAPCOSTS parameter)
	 * 
	 * @return an integer value for gap creation used by this search, -1 if not set or not a number
	 */
	public int getBlastGapCreation() {
		String gapCosts = param.get(GAPCOSTS);
		try {
			String gapCreation = gapCosts.split("+")[0];
			return Integer.parseInt(gapCreation);
		} catch (Exception e) {
			return -1;
		}
	}

	/**
	 * Returns the value for the gap extension parameter (second half of GAPCOSTS parameter)
	 * 
	 * @return an integer for the value for gap extension used by this search, -1 if not set or not a number
	 */
	public int getBlastGapExtension() {
		String gapCosts = param.get(GAPCOSTS);
		try {
			String gapExtension = gapCosts.split("+")[1];
			return Integer.parseInt(gapExtension);
		} catch (Exception e) {
			return -1;
		}
	}

	/**
	 * Returns the actual string for the GAPCOSTS parameter which is used to build the URL
	 * 
	 * @return the string representation of the GAPCOSTS parameter formatted for the URL
	 */
	public String getBlastGapCosts() {
		return param.get(GAPCOSTS);
	}

	/**
	 * Sets the GAPCOSTS parameter
	 * 
	 * @param gapCreation integer to use as gap creation value
	 * @param gapExtension integer to use as gap extension value
	 */
	public void setBlastGapCosts(int gapCreation, int gapExtension) {
		String gc = Integer.toString(gapCreation);
		String ge = Integer.toString(gapExtension);
		param.put(GAPCOSTS, gc + "+" + ge);
	}

	/**
	 * Returns the value of the specified substitution matrix
	 * 
	 * @return matrix: the name of the specified substitution matrix
	 */
	public String getBlastMatrix() {
		return param.get(MATRIX_NAME);
	}

	/**
	 * Sets the value for the MATRIX parameter to use for blastall
	 * <p>
	 * Blastall equivalent: -M
	 * 
	 * @param matrix : a String to use as gap creation value
	 * @see BlastMatrix
	 */
	public void setBlastMatrix(BlastMatrix matrix) {
		param.put(MATRIX_NAME, matrix.name());

		boolean gapCostsSet = getBlastGapCreation() != -1 || getBlastGapExtension() != -1;

		if (!gapCostsSet) {
			/*
			 * Setting default values for -G/-E if no other values have been set is necessary because, since BLOSUM62 is
			 * default, the expected values are -G 11 -E 1. If your matrix choice is different, the request will fail,
			 * implicitly expecting GAPCOSTS=11+1
			 */
			switch (matrix) {
			case PAM30:
				setBlastGapCosts(9, 1);
				break;
			case PAM70:
				setBlastGapCosts(10, 1);
				break;
			case PAM250:
				setBlastGapCosts(14, 2);
				break;
			case BLOSUM45:
				setBlastGapCosts(15, 2);
				break;
			case BLOSUM50:
				setBlastGapCosts(13, 2);
				break;
			case BLOSUM80:
			case BLOSUM90:
				setBlastGapCosts(10, 1);
				break;
			}
		}
	}

	/**
	 * Sets the QUERY_FROM and QUERY_TO parameters to be use by blast. Do not use if you want to use the whole sequence
	 * 
	 * @param start QUERY_FROM parameter
	 * @param end QUERY_TO parameter
	 * @see #setBlastFromPosition(int)
	 * @see #setBlastToPosition(int)
	 */
	public void setBlastFromToPosition(int start, int end) {
		if (start >= end) {
			throw new IllegalArgumentException("Start index must be less than end index");
		}
		setBlastFromPosition(start);
		setBlastToPosition(end);
	}

	/**
	 * @return an integer value for the QUERY_FROM parameter
	 * @see #setBlastFromToPosition(int, int)
	 */
	public int getBlastFromPosition() {
//		int a = 0;
//		if (this.param.get("QUERY_FROM") != "-1")
//			a = Integer.parseInt(this.param.get("QUERY_FROM"));
//		else if (this.param.get("QUERY_FROM") == "-1")
//			a = -1;
//		return a;
		return Integer.parseInt(param.get(QUERY_FROM));
	}

	/**
	 * This method set the QUERY_FROM parameter to be use by blast. If you decide to use the whole sequence, do not
	 * use...
	 * <p>
	 * Blastall equivalent: -L
	 * 
	 * @param start : an integer to use for QUERY_FROM
	 */
	public void setBlastFromPosition(int start) {
		param.put("QUERY_FROM", String.valueOf(start));
	}

	/**
	 * @return QUERY_TO parameter
	 * @see #setBlastFromToPosition(int, int)
	 */
	public int getBlastToPosition() {
//		int a = 0;
//		if (this.param.get("QUERY_TO") != "-1")
//			a = Integer.parseInt(this.param.get("QUERY_TO"));
//		else if (this.param.get("QUERY_TO") == "-1")
//			a = -1;
//		return a;
		return Integer.parseInt(param.get(QUERY_TO));
	}

	/**
	 * Sets the QUERY_TO parameter to be use by blast. If you decide to use the whole sequence, do not use...
	 * <p>
	 * Blastall equivalent: -L
	 * 
	 * @param stop : an integer to use as QUERY_TO
	 */
	public void setBlastToPosition(int stop) {
		this.param.put("QUERY_TO", String.valueOf(stop));
	}

	/**
	 * This method is to be used if a request is to use non-default values at submission. Useful for the following
	 * blastall parameters:
	 * <ul>
	 * <li>-r: integer to reward for match. Default = 1</li>
	 * <li>-q: negative integer for penalty to allow mismatch. Default = -3</li>
	 * <li>-y: dropoff for blast extensions in bits, using default if not specified. Default = 20 for blastn, 7 for all
	 * others (except megablast for which it is not applicable).</li>
	 * <li>-X: X dropoff value for gapped alignment, in bits. Default = 30 for blastn/megablast, 15 for all others.</li>
	 * <li>-Z: final X dropoff value for gapped alignement, in bits. Default = 50 for blastn, 25 for all others (except
	 * megablast for which it is not applicable)</li>
	 * <li>-P: equals 0 for multiple hits 1-pass, 1 for single hit 1-pass. Does not apply to blastn ou megablast.</li>
	 * <li>-A: multiple hits window size. Default = 0 (for single hit algorithm)</li>
	 * <li>-I: number of database sequences to save hits for. Default = 500</li>
	 * <li>-Y: effective length of the search space. Default = 0 (0 represents using the whole space)</li>
	 * <li>-z: a real specifying the effective length of the database to use. Default = 0 (0 represents the real size)</li>
	 * <li>-c: an integer representing pseudocount constant for PSI-BLAST. Default = 7</li>
	 * <li>-F: any filtering directive</li>
	 * </ul>
	 * <p>
	 * WARNING!! This method is still very much in flux and might not work as expected...
	 * </p>
	 * <p>
	 * You have to be aware that at no moment is there any error checking on the use of these parameters by this class.
	 * </p>
	 * 
	 * @param advancedOptions : a String with any number of optional parameters with an associated value.
	 */
	public void setBlastAdvancedOptions(String advancedOptions) {
		// Escaping white spaces with + char to comply with QBlast specifications
		param.put(OTHER_ADVANCED, advancedOptions.replaceAll(" ", "+"));
	}

	/**
	 * @return the String with the advanced options
	 */
	public String getBlastAdvancedOptions() {
		return param.get(OTHER_ADVANCED);
	}

	@Override
	public String getAlignmentOption(String key) {
		if (!param.containsKey(key)) {
			throw new IllegalArgumentException("The key named " + key + " is not set for this alignment properties object");
		}
		return param.get(key);
	}

	@Override
	public void setAlignementOption(String key, String val) {
		param.put(key, val);
	}

	@Override
	public Set<String> getAlignmentOptions() {
		return param.keySet();
	}
}
