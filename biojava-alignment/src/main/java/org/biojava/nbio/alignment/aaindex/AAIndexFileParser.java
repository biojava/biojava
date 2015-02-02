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
package org.biojava.nbio.alignment.aaindex;

import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;



public class AAIndexFileParser {


	int scale = -1;

	Map<String,SubstitutionMatrix<AminoAcidCompound>> matrices;

	ScaledSubstitutionMatrix currentMatrix;
	String currentRows;
	String currentCols;
	int currentRowPos;
	private short[][] matrix;
	short max;
	short min;
	private List<AminoAcidCompound> rows, cols;
	boolean inMatrix;
	boolean symmetricMatrix ;


	public AAIndexFileParser(){
		matrices  = new HashMap<String, SubstitutionMatrix<AminoAcidCompound>>();
	}

	/** parse an inputStream that points to an AAINDEX database file
	 * 
	 * @param inputStream
	 * @throws IOException 
	 */
	public void parse(InputStream inputStream) throws IOException {

		currentMatrix = null;
		currentRows = "";
		currentCols = "";
		max = Short.MIN_VALUE;
		min = Short.MAX_VALUE;
		inMatrix = false;
		
		BufferedReader buf = new BufferedReader (new InputStreamReader (inputStream));
		String line = null;
		line = buf.readLine();

		while (  line != null ) {

			if ( line.startsWith("//")) {
				finalizeMatrix();
				inMatrix = false;

			} else if ( line.startsWith("H ")){
				// a new matric!
				newMatrix(line);
			} else if ( line.startsWith("D ")) {
				currentMatrix.setDescription(line.substring(2));
			} else if ( line.startsWith("M ")) {
				initMatrix(line);
				inMatrix = true;
			} else if ( line.startsWith("  ")){
				if ( inMatrix)
					processScores(line);
			}

			line = buf.readLine();
		}

	}


	//  process a line such as >    -0.3     1.6     0.7     0.8    -2.6     3.0<
	private void processScores(String line) {

		String[] values = line.trim().split(" +");

		// increment the current row we are talking about
		currentRowPos++;

		
		
		for ( int i =0 ; i < values.length ; i++){
		
			if ( values[i].endsWith(".")) {
				values[i] = values[i] + "0";
			}

			// special case: MEHP950101
			if (values[i].equals("-")) {
				values[i] = "0";
			}
		
			if ( scale == -1 ) {
				scale = determineScale(values[0]);				
			}
			
			
			Float score = Float.parseFloat(values[i]);
			score = scale * score;

			Short s = (short) Math.round(score);
			
			matrix[currentRowPos][i] = s;

			if ( values.length < cols.size() || ( symmetricMatrix)){
				//System.out.println(values.length + " " + cols.size() + " " + currentRowPos + " " + i + " " +  line);
				
				matrix[i][currentRowPos] = s;
				
				symmetricMatrix = true;
				
			}

			if ( score > max)
				max = s;
			if ( score < min)
				min = s;


		}
	}

	private int determineScale(String value) {
		
		String[] spl = value.split("\\.");
		
		if (spl.length <= 1)
			return 1;
		
		String digits = spl[1];
		
		return (int)Math.round(Math.pow(10, digits.length()));
		
	}

	// process a line of type >M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV<
	private void initMatrix(String line) {
		String[] spl = line.split(" ");

		// trim off the final , character 
		currentRows = spl[3].substring(0, spl[3].length()-1);
		currentCols = spl[6];
		currentRowPos = -1;

		int nrRows = currentRows.length();
		int nrCols = currentCols.length();

		matrix = new short[nrRows][nrCols];

		rows = new ArrayList<AminoAcidCompound>();
		cols = new ArrayList<AminoAcidCompound>();


		//System.out.println(">" + currentRows+"<");
		AminoAcidCompoundSet compoundSet = AminoAcidCompoundSet.getAminoAcidCompoundSet();
		for ( int i = 0 ; i < currentRows.length() ; i ++){
			char c = currentRows.charAt(i);
			AminoAcidCompound aa = compoundSet.getCompoundForString(c+""); 

			rows.add(aa);
		}

		for ( int i = 0 ; i < currentCols.length() ; i ++){
			char c = currentRows.charAt(i);
			AminoAcidCompound aa = compoundSet.getCompoundForString(c+""); 

			cols.add(aa);
		}



	

		currentMatrix.setScale(scale);
	}

	
	private void newMatrix(String line) {
		symmetricMatrix = false;
		scale = -1;
		
		currentMatrix = new ScaledSubstitutionMatrix();
		currentMatrix.setName(line.substring(2));
		
		
		//System.out.println("new Matrix " + currentMatrix.getName());
	}

	// 
	private SubstitutionMatrix<AminoAcidCompound> finalizeMatrix() {

		currentMatrix.setMatrix(matrix);
		currentMatrix.setMax(max);
		currentMatrix.setMin(min);
		currentMatrix.setCols(cols);
		currentMatrix.setRows(rows);
		currentMatrix.setScale(scale);
		matrices.put(currentMatrix.getName(), currentMatrix);

		return currentMatrix;

	}

	public Map<String, SubstitutionMatrix<AminoAcidCompound>> getMatrices() {
		return matrices;
	}
}
