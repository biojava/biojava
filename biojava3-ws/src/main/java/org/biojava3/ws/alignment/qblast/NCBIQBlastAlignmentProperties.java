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

import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.DATABASE;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.EXPECT;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.GAPCOSTS;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.MATRIX_NAME;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.MEGABLAST;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.OTHER_ADVANCED;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.PROGRAM;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.QUERY_FROM;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.QUERY_TO;
import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.WORD_SIZE;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.biojava3.ws.alignment.RemotePairwiseAlignmentProperties;

/**
 * This class wraps a QBlast search request parameter {@code Map} by adding several convenient parameter addition
 * methods. Other QBlast URL API parameters should be added using
 * {@link #setAlignmentOption(BlastAlignmentParameterEnum, String)}
 * <p/>
 * Required parameters are {@code PROGRAM} and {@code DATABASE}, other parameters are optional
 * 
 * @author Sylvain Foisy, Diploide BioIT
 * @author Gediminas Rimsa
 */
public class NCBIQBlastAlignmentProperties implements RemotePairwiseAlignmentProperties {
	private static final long serialVersionUID = 7158270364392309841L;

	private Map<BlastAlignmentParameterEnum, String> param = new HashMap<BlastAlignmentParameterEnum, String>();

	/**
	 * This method forwards to {@link #getAlignmentOption(BlastAlignmentParameterEnum)}. Consider using it instead.
	 */
	@Override
	public String getAlignmentOption(String key) {
		return getAlignmentOption(BlastAlignmentParameterEnum.valueOf(key));
	}

	/**
	 * This method forwards to {@link #setAlignmentOption(BlastAlignmentParameterEnum, String)}. Consider using it
	 * instead.
	 */
	@Override
	public void setAlignementOption(String key, String val) {
		setAlignmentOption(BlastAlignmentParameterEnum.valueOf(key), val);
	}

	/**
	 * Gets parameters, which are currently set
	 */
	public Set<String> getAlignmentOptions() {
		Set<String> result = new HashSet<String>();
		for (BlastAlignmentParameterEnum parameter : param.keySet()) {
			result.add(parameter.name());
		}
		return result;
	}

	/**
	 * Gets the value of specified parameter or {@code null} if it is not set.
	 */
	public String getAlignmentOption(BlastAlignmentParameterEnum key) {
		return param.get(key);
	}

	/**
	 * Sets the value of specified parameter
	 */
	public void setAlignmentOption(BlastAlignmentParameterEnum key, String val) {
		param.put(key, val);
	}

	/**
	 * Removes given parameter
	 */
	public void removeAlignmentOption(BlastAlignmentParameterEnum key) {
		param.remove(key);
	}

	/**
	 * @return {@linkplain BlastProgramEnum} used for blast run
	 */
	public BlastProgramEnum getBlastProgram() {
		BlastProgramEnum program = BlastProgramEnum.valueOf(getAlignmentOption(PROGRAM));
		boolean isMegablast = BlastProgramEnum.blastn == program && getAlignmentOption(MEGABLAST).equals("on");
		return !isMegablast ? program : BlastProgramEnum.megablast;
	}

	/**
	 * Sets the program to be used with blastall
	 * 
	 * @param program : one of blastall programs
	 */
	public void setBlastProgram(BlastProgramEnum program) {
		if (BlastProgramEnum.megablast != program) {
			setAlignmentOption(PROGRAM, program.name());
			removeAlignmentOption(MEGABLAST);
		} else {
			setAlignmentOption(PROGRAM, BlastProgramEnum.blastn.name());
			setAlignmentOption(MEGABLAST, "on");
		}
	}

	/**
	 * @return name of database used with blastall
	 */
	public String getBlastDatabase() {
		return getAlignmentOption(DATABASE);
	}

	/*
	 * TODO: update comment when URL is available:
	 * A quite exhaustive list of the databases available for QBlast
	 * requests can be found here: <p> http://&lt to_be_completed &gt <p> Blastall equivalent: -d
	 */

	/**
	 * Sets the database to be used with blastall
	 * <p>
	 * A list of available databases can be acquired by calling {@link NCBIQBlastService#printRemoteBlastInfo()}
	 * <p>
	 * Blastall equivalent: -d
	 * 
	 * @param db : a valid name to a NCBI blastable database
	 */
	public void setBlastDatabase(String database) {
		setAlignmentOption(DATABASE, database);
	}

	/**
	 * @return double value of EXPECT parameter used for blast run
	 */
	public double getBlastExpect() {
		if (param.containsKey(EXPECT)) {
			return Double.parseDouble(getAlignmentOption(EXPECT));
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
		setAlignmentOption(EXPECT, Double.toString(expect));
	}

	/**
	 * Returns the value of the WORD_SIZE parameter used for this blast run
	 * 
	 * @return int value of WORD_SIZE used by this search
	 * @throws IllegalArgumentException when program type is not set and program type is not supported
	 */
	public int getBlastWordSize() {
		if (param.containsKey(WORD_SIZE)) {
			return Integer.parseInt(getAlignmentOption(WORD_SIZE));
		}

		// return default word size value
		try {
			BlastProgramEnum programType = getBlastProgram();
			switch (programType) {
			case blastn:
				return 11;
			case megablast:
				return 28;
			case blastp:
			case blastx:
			case tblastn:
			case tblastx:
				return 3;
			default:
				throw new UnsupportedOperationException("Blast program " + programType.name() + " is not supported.");
			}
		} catch (IllegalArgumentException e) {
			throw new IllegalArgumentException("Blast program " + getBlastProgram() + " is not supported.", e);
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
		setAlignmentOption(WORD_SIZE, Integer.toString(word));
	}

	/**
	 * Returns the value for the GAP_CREATION parameter (first half of GAPCOSTS parameter)
	 * 
	 * @return an integer value for gap creation used by this search, -1 if not set or not a number
	 */
	public int getBlastGapCreation() {
		String gapCosts = getAlignmentOption(GAPCOSTS);
		try {
			String gapCreation = gapCosts.split("\\+")[0];
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
		String gapCosts = getAlignmentOption(GAPCOSTS);
		try {
			String gapExtension = gapCosts.split("\\+")[1];
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
		return getAlignmentOption(GAPCOSTS);
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
		setAlignmentOption(GAPCOSTS, gc + "+" + ge);
	}

	/**
	 * Returns the value of the specified substitution matrix
	 * 
	 * @return matrix: the name of the specified substitution matrix
	 */
	public String getBlastMatrix() {
		return getAlignmentOption(MATRIX_NAME);
	}

	/**
	 * Sets the value for the MATRIX parameter to use for blastall
	 * <p>
	 * Blastall equivalent: -M
	 * 
	 * @param matrix : a String to use as gap creation value
	 * @see BlastMatrixEnum
	 */
	public void setBlastMatrix(BlastMatrixEnum matrix) {
		setAlignmentOption(MATRIX_NAME, matrix.name());

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
			case BLOSUM62:
				break;
			default:
				break;
			}
		}
	}

	/**
	 * Sets the QUERY_FROM and QUERY_TO parameters to be use by blast. Do not use if you want to use the whole sequence.<br/>
	 * Blastall equivalent: -L
	 * 
	 * @param start QUERY_FROM parameter
	 * @param end QUERY_TO parameter
	 */
	public void setBlastFromToPosition(int start, int end) {
		if (start >= end) {
			throw new IllegalArgumentException("Start index must be less than end index");
		}
		setAlignmentOption(QUERY_FROM, String.valueOf(start));
		setAlignmentOption(QUERY_TO, String.valueOf(end));
	}

	/**
	 * @return an integer value for the QUERY_FROM parameter
	 * @see #setBlastFromToPosition(int, int)
	 */
	public int getBlastFromPosition() {
		return Integer.parseInt(getAlignmentOption(QUERY_FROM));
	}

	/**
	 * @return QUERY_TO parameter
	 * @see #setBlastFromToPosition(int, int)
	 */
	public int getBlastToPosition() {
		return Integer.parseInt(getAlignmentOption(QUERY_TO));
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
		setAlignmentOption(OTHER_ADVANCED, advancedOptions.replaceAll(" ", "+"));
	}

	/**
	 * @return the String with the advanced options
	 */
	public String getBlastAdvancedOptions() {
		return getAlignmentOption(OTHER_ADVANCED);
	}
}
