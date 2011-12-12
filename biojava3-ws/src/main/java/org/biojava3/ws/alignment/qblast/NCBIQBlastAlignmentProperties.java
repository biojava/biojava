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

import java.lang.Integer;
import java.util.HashMap;
import java.util.Set;
import java.util.Arrays;

import org.biojava3.ws.alignment.RemotePairwiseAlignmentProperties;

/**
 * 
 * This class implements RemotePairwiseAlignmentProperties by specifying several
 * convenient methods used to wrap the addition of Blast alignment parameters.
 * 
 * It is responsible for collecting, doing some basic sanity checks and to
 * create the part of the URL that will hold all alignment parameters to be
 * used.
 * 
 * <p>
 * Many thanks to Matthew Busse for helping in debugging after the migration
 * from BJ1.7 to BJ3.0.
 * </p>
 * 
 * @author Sylvain Foisy, Diploide BioIT
 * @since Biojava 3
 * 
 */
public class NCBIQBlastAlignmentProperties implements
		RemotePairwiseAlignmentProperties {

	private static final long serialVersionUID = 1L;
	private HashMap<String, String> param = new HashMap<String, String>();
	private String cmd;

	/**
	 * 
	 * Initialization of default values
	 * 
	 */
	public NCBIQBlastAlignmentProperties() {
		// Only mandatory parameters
		this.param.put("PROGRAM", "not_set");
		this.param.put("DATABASE", "not_set");
		// Optional parameters to set
		this.param.put("WORD_SIZE", "-1");
		this.param.put("EXPECT", "-1");
		this.param.put("QUERY_FROM", "-1");
		this.param.put("QUERY_TO", "-1");
		this.param.put("MATRIX_NAME", "BLOSUM62");
		this.param.put("GAP_CREATION", "-1");
		this.param.put("GAP_EXTENSION", "-1");
		// Everything else
		this.param.put("OTHER_ADVANCED", "not_set");

		cmd = "CMD=Put";
	}

	/**
	 * This method returns the value of the program used for this particular
	 * blast run.
	 * 
	 * @return program:the name of the blastall program used.
	 * 
	 */
	public String getBlastProgram() {
		return this.param.get("PROGRAM");
	}

	/**
	 * This method set the program to be use with pairwise alignment. This
	 * method does a validation before running on the valid blastall programs:
	 * blastn / megablast / blastp / blastx / tblastn / tblastx
	 * 
	 * @param program
	 *            : one of blastall programs
	 * 
	 * @exception Exception
	 *                if the named program is not a valid blastall options
	 * 
	 */
	public void setBlastProgram(String program) throws Exception {

		boolean isValid = false;
		String[] blastPr = new String[] { "blastn", "blastp", "blastx",
				"megablast", "tblastn", "tblastx" };

		/*
		 * To check if the program called for belongs to the blastPr array
		 */
		if (Arrays.binarySearch(blastPr, program) >= 0) {
			if (program != "megablast")
				this.param.put("PROGRAM", program);
			else
				this.param.put("PROGRAM", "blastn&MEGABLAST=on");
			isValid = true;
		}

		if (!isValid) {
			throw new Exception(
					"Invalid blastall program selection! Use one of valid values: blastn/blastp/blastx/tblastn/tblastx");
		}
	}

	/**
	 * This method returns the value of the database used for this particular
	 * blast run.
	 * 
	 * @return db :the name of the database used
	 */
	public String getBlastDatabase() {
		return this.param.get("DATABASE");
	}

	/**
	 * This method set the database to be used with blastall
	 * 
	 * A quite exhaustive list of the databases available for QBlast requests
	 * can be found here:
	 * 
	 * http://<to_be_completed>
	 * 
	 * Blastall equivalent: -d
	 * 
	 * @param db
	 *            :a valid name to a NCBI blastable database
	 */
	public void setBlastDatabase(String db) {
		this.param.put("DATABASE", db);
	}

	/**
	 * This method returns the value of EXPECT parameter used for this
	 * particular blast run.
	 * 
	 * @return double :the value for EXPECT used by this search
	 */
	public double getBlastExpect() {
		if (this.param.get("EXPECT") != "-1")
			return Double.parseDouble(this.param.get("EXPECT"));
		else
			return 10;
	}

	/**
	 * This method set the EXPECT parameter to be use with blastall
	 * 
	 * Example: if you want a EXPECT of 1e-10, try this:
	 * 
	 * setBlastExpect(Double.parseDouble("1e-10"))
	 * 
	 * Blastall equivalent: -e
	 * 
	 * @param expect
	 *            : a double used to set EXPECT
	 */
	public void setBlastExpect(double expect) {
		String str = Double.toString(expect);
		this.param.put("EXPECT", str);
	}

	/**
	 * This method returns the value of the WORD_SIZE parameter used for this
	 * particular blast run.
	 * 
	 * @return word :the value for WORD_SIZE used by this search
	 */
	public int getBlastWordSize() {

		int word = -1;

		if (this.param.get("WORD_SIZE") != "-1")
			word = Integer.parseInt(this.param.get("WORD_SIZE"));
		else {
			if (this.param.get("PROGRAM") == "blastn")
				word = 11;
			else if (this.param.get("PROGRAM") == "blastp"
					|| this.param.get("PROGRAM") == "blastx"
					|| this.param.get("PROGRAM") == "tblastn"
					|| this.param.get("PROGRAM") == "tblastx")
				word = 3;
			else if (this.param.get("PROGRAM") == "blastn&MEGABLAST=on")
				word = 28;
		}
		return word;
	}

	/**
	 * This method set the WORD_SIZE parameter to be use with blastall.
	 * 
	 * WARNING!! At this point, the method does not verify the validity of your
	 * choice; for example, word size of greater than 5 with blastp returns
	 * error messages from QBlast. Word size range depends on the algorithm
	 * chosen.
	 * 
	 * More at http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node74.html
	 * 
	 * Blastall equivalent: -W
	 * 
	 * @param word
	 *            : an integer used to set WORD_SIZE
	 */
	public void setBlastWordSize(int word) {
		this.param.put("WORD_SIZE", Integer.toString(word));
	}

	/**
	 * 
	 * This method returns the value for the GAP_CREATION parameter
	 * 
	 * @return g :the value for GAP_CREATION used by this search
	 * 
	 */
	public int getBlastGapCreation() {
		return Integer.parseInt(this.param.get("GAP_CREATION"));
	}

	/**
	 * This method set the gap creation value for the GAPCOST parameter to use
	 * for blastall
	 * 
	 * Blastall equivalent : -G
	 * 
	 * @param g
	 *            : an integer to use as gap creation value
	 * 
	 */
	public void setBlastGapCreation(int g) {
		this.param.put("GAP_CREATION", Integer.toString(g));
	}

	/**
	 * 
	 * This method returns the value for the GAP_EXTENSION parameter
	 * 
	 * @return e : an integer for the value for GAP_EXTENSION used by this
	 *         search
	 * 
	 */
	public int getBlastGapExtension() {
		return Integer.parseInt(this.param.get("GAP_EXTENSION"));
	}

	/**
	 * This method set the gap extension value for the GAPCOST parameter to use
	 * for blastall
	 * 
	 * Blastall equivalent: -E
	 * 
	 * @param e
	 *            : an integer to use as gap extension value
	 */
	public void setBlastGapExtension(int e) {
		this.param.put("GAP_EXTENSION", Integer.toString(e));
	}

	/**
	 * This method return the actual string for the GAPCOSTS parameter which is
	 * used to build the URL
	 * 
	 * @return str : the string representation of the GAPCOSTS parameter
	 *         formatted for the URL
	 * 
	 */
	public String getBlastGapCosts() {
		if (this.param.get("GAPCOSTS") != null)
			return this.param.get("GAPCOSTS");
		else
			return "defaults";
	}

	/**
	 * 
	 * This method process the values from GAP_CREATION and GAP_EXTENSION to
	 * generate the values for the actual GAPCOSTS parameter
	 * 
	 */
	private void setBlastGapCosts() {

		String gc = Integer.toString(this.getBlastGapCreation());
		String ge = Integer.toString(this.getBlastGapExtension());

		this.param.put("GAPCOSTS", gc + "+" + ge);
	}

	/**
	 * 
	 * This method returns the value of the specified substitution matrix
	 * 
	 * @return matrix: the name of the specified substitution matrix
	 */
	public String getBlastMatrix() {
		return this.param.get("MATRIX_NAME");
	}

	/**
	 * This method set the value for the MATRIX parameter to use for blastall
	 * 
	 * Allowed matrices:
	 * PAM30,PAM70,PAM90,PAM250,BLOSUM45,BLOSUM50,BLOSUM62,BLOSUM80 (Other
	 * matrices are not useable via QBlast)
	 * 
	 * Blastall equivalent: -M
	 * 
	 * @param mtx
	 *            : a String to use as gap creation value
	 * 
	 * @throws Exception
	 *             if matrix name is not part of allowed BLAST matrices
	 */
	public void setBlastMatrix(String mtx) throws Exception {
		boolean isValid = false;
		String[] blastMat = new String[] { "BLOSUM45", "BLOSUM50", "BLOSUM62",
				"BLOSUM80", "BLOSUM90", "PAM250", "PAM30", "PAM70" };

		/*
		 * To check if the matrix called for belongs to the blastMat array
		 */
		if (Arrays.binarySearch(blastMat, mtx) >= 0) {
			this.param.put("MATRIX_NAME", mtx);
			isValid = true;

			/*
			 * This step is necessary because, since BLOSUM62 is default, the
			 * expected values are -G 11 -E 1. If your matrix choice is
			 * different, the request will fail, implicitly expecting
			 * &GAPCOSTS=11+1
			 */
			if (mtx != "BLOSUM62") {
				/*
				 * Setting default values for -G/-E if no other values have been
				 * set via setBlastGapCreation/setBlastGapExtension
				 */
				if (mtx == "PAM30") {
					if (this.getBlastGapCreation() == -1
							&& this.getBlastGapExtension() == -1) {
						this.setBlastGapCreation(9);
						this.setBlastGapExtension(1);
					}
				} else if (mtx == "PAM70") {
					if (this.getBlastGapCreation() == -1
							&& this.getBlastGapExtension() == -1) {
						this.setBlastGapCreation(10);
						this.setBlastGapExtension(1);
					}

				} else if (mtx == "PAM250") {
					if (this.getBlastGapCreation() == -1
							&& this.getBlastGapExtension() == -1) {
						this.setBlastGapCreation(14);
						this.setBlastGapExtension(2);
					}

				} else if (mtx == "BLOSUM45") {
					if (this.getBlastGapCreation() == -1
							&& this.getBlastGapExtension() == -1) {
						this.setBlastGapCreation(15);
						this.setBlastGapExtension(2);
					}
				} else if (mtx == "BLOSUM50") {
					if (this.getBlastGapCreation() == -1
							&& this.getBlastGapExtension() == -1) {
						this.setBlastGapCreation(13);
						this.setBlastGapExtension(2);
					}
				} else if (mtx == "BLOSUM80") {
					if (this.getBlastGapCreation() == -1
							&& this.getBlastGapExtension() == -1) {
						this.setBlastGapCreation(10);
						this.setBlastGapExtension(1);
					}
				} else if (mtx == "BLOSUM90") {
					if (this.getBlastGapCreation() == -1
							&& this.getBlastGapExtension() == -1) {
						this.setBlastGapCreation(10);
						this.setBlastGapExtension(1);
					}
				}
			}
		}
		
		if (!isValid)
			throw new Exception(
					"Invalid blastp substitution matrix selection! Use one of valid values: PAM30,PAM70,PAM250,BLOSUM45,BLOSUM50,BLOSUM62,BLOSUM80\n");
	}

	/**
	 * This method returns the value for the QUERY_FROM parameter
	 * 
	 * @return an integer value for the QUERY_FROM parameter
	 * 
	 */
	public int getBlastFromPosition() {
		int a = 0;
		if (this.param.get("QUERY_FROM") != "-1")
			a = Integer.parseInt(this.param.get("QUERY_FROM"));
		else if (this.param.get("QUERY_FROM") == "-1")
			a = -1;
		return a;
	}

	/**
	 * 
	 * This method set the QUERY_FROM parameter to be use by blast. It needs the
	 * corresponding setBlastToPosition() to work. If you decide to use the
	 * whole sequence, do not use...
	 * 
	 * Blastall equivalent: -L
	 * 
	 * @param start
	 *            : an integer to use for QUERY_FROM
	 * 
	 */
	public void setBlastFromPosition(int start) {
		this.param.put("QUERY_FROM", String.valueOf(start));
	}

	/**
	 * 
	 * This method returns the value for the QUERY_TO parameter
	 * 
	 * @return an integer value for the QUERY_TO parameter
	 * 
	 */
	public int getBlastToPosition() {
		int a = 0;
		if (this.param.get("QUERY_TO") != "-1")
			a = Integer.parseInt(this.param.get("QUERY_TO"));
		else if (this.param.get("QUERY_TO") == "-1")
			a = -1;
		return a;
	}

	/**
	 * 
	 * This method set the QUERY_TO parameter to be use by blast. It needs the
	 * corresponding setBlastFromPosition(). If you decide to use the whole
	 * sequence, do not use...
	 * 
	 * Blastall equivalent: -L
	 * 
	 * @param stop
	 *            : an integer to use as QUERY_TO
	 */
	public void setBlastToPosition(int stop) {
		this.param.put("QUERY_TO", String.valueOf(stop));
	}

	/**
	 * 
	 * This method is to be used if a request is to use non-default values at
	 * submission. According to QBlast info, the accepted parameters for PUT
	 * requests are:
	 * 
	 * 
	 * Useful for the following blastall parameters:
	 * 
	 * <ul>
	 * <li>-r: integer to reward for match. Default = 1</li>
	 * <li>-q: negative integer for penalty to allow mismatch. Default = -3</li>
	 * <li>-y: dropoff for blast extensions in bits, using default if not
	 * specified. Default = 20 for blastn, 7 for all others (except megablast
	 * for which it is not applicable).</li>
	 * <li>-X: X dropoff value for gapped alignment, in bits. Default = 30 for
	 * blastn/megablast, 15 for all others.</li>
	 * <li>-Z: final X dropoff value for gapped alignement, in bits. Default =
	 * 50 for blastn, 25 for all others (except megablast for which it is not
	 * applicable)</li>
	 * <li>-P: equals 0 for multiple hits 1-pass, 1 for single hit 1-pass. Does
	 * not apply to blastn ou megablast.</li>
	 * <li>-A: multiple hits window size. Default = 0 (for single hit algorithm)
	 * </li>
	 * <li>-I: number of database sequences to save hits for. Default = 500</li>
	 * <li>-Y: effective length of the search space. Default = 0 (0 represents
	 * using the whole space)</li>
	 * <li>-z: a real specifying the effective length of the database to use.
	 * Default = 0 (0 represents the real size)</li>
	 * <li>-c: an integer representing pseudocount constant for PSI-BLAST.
	 * Default = 7</li>
	 * <li>-F: any filtering directive</li>
	 * </ul>
	 * 
	 * <p>
	 * WARNING!! This method is still very much in flux and might not work as
	 * expected...
	 * </p>
	 * <p>
	 * You have to be aware that at no moment is there any error checking on the
	 * use of these parameters by this class.
	 * </p>
	 * 
	 * @param aStr
	 *            : a String with any number of optional parameters with an
	 *            associated value.
	 * 
	 */
	public void setBlastAdvancedOptions(String aStr) {

		// Escaping white spaces with + char to
		// comply with QBlast specifications
		String tmp = aStr.replaceAll(" ", "+");

		this.param.put("OTHER_ADVANCED", tmp);
	}

	/**
	 * 
	 * Simply return the string given as argument via setBlastAdvancedOptions
	 * 
	 * @return advanced :the string with the advanced options
	 */
	public String getBlastAdvancedOptions() {
		return this.param.get("OTHER_ADVANCED");
	}

	/**
	 * 
	 * This method will return the part of the URL submitted to QBlast that has
	 * all the alignment parameters defined for this alignment.
	 * 
	 * Warning!! It does not contain any sequence, gid or Genbank Identifier
	 * These are added in the workings of the NCBIQBlastService class.
	 * 
	 * @return cmd: part of the URL used for this alignment request
	 * 
	 */
	public String getBlastCommandsToQBlast() {
		return this.cmd;
	}

	/**
	 * 
	 * This method is responsible for building the String with all the alignment
	 * parameters to use with a given request via your program
	 * 
	 * It does basic sanity checks on the values but nothing else at this
	 * point...
	 * 
	 */
	public void setBlastCommandsToQBlast() throws Exception {

		/*
		 * blastall program has to be set...
		 */
		if (this.getBlastProgram() == "not_set") {
			throw new Exception(
					"Impossible to execute QBlast request. Your program has not been set correctly.\n");
		} else {
			this.cmd = this.cmd + "&PROGRAM=" + this.getBlastProgram();
		}

		/*
		 * A database has to be specified...
		 */
		if (this.getBlastDatabase() == "not_set") {
			throw new Exception(
					"Impossible to execute QBlast request. Your database has not been set correctly.\n");
		} else {
			this.cmd = this.cmd + "&DATABASE=" + this.getBlastDatabase();
		}

		/*
		 * This code block deals with what to do with the various non-mandatory
		 * parameters
		 */
		if (this.getBlastExpect() != -1) {
			this.cmd = this.cmd + "&EXPECT=" + this.getBlastExpect();
		}

		if (this.getBlastWordSize() != -1) {
			this.cmd = this.cmd + "&WORD_SIZE="
					+ this.getAlignmentOption("WORD_SIZE");
		}

		if (this.getBlastFromPosition() != -1
				&& this.getBlastToPosition() != -1) {
			this.cmd = this.cmd + "&QUERY_FROM=" + this.getBlastFromPosition()
					+ "&QUERY_TO=" + this.getBlastToPosition();
		}

		if (this.getBlastProgram() != "blastn") {
			if (this.getBlastMatrix() != "BLOSUM62") {
				this.cmd = this.cmd + "&MATRIX_NAME=" + this.getBlastMatrix();
				if (this.getBlastGapCreation() != -1
						&& this.getBlastGapCreation() != -1) {
					this.setBlastGapCosts();
					cmd = cmd + "&GAPCOSTS=" + this.getBlastGapCosts();
				}
			}
		}
		// if (this.getBlastAdvancedOptions()!="not_set") {
		// cmd = cmd + "&OTHER_ADVANCED=" +
		// this.getAlignmentOption("OTHER_ADVANCED");
		// }
	}

	/**
	 * 
	 * A way to start a new URL from scratch after analyzing one sequence before
	 * going to the next if new parameters are necessary
	 * 
	 */
	public void reinitializeBlastCommandsToQBlast() {
		this.cmd = "CMD=Put";
	}

	/*
	 * 
	 * These three methods are necessary to comply with Interface definition
	 * 
	 * Could be useful
	 */
	public String getAlignmentOption(String key) throws Exception {
		if (param.containsKey(key)) {
			return this.param.get(key);
		} else {
			throw new Exception("The key named " + key
					+ " is not set in this RemoteQBlastOutputProperties object");
		}
	}

	public void setAlignementOption(String key, String val) {
		this.param.put(key, val);
	}

	public Set<String> getAlignmentOptions() {
		return param.keySet();
	}
}
