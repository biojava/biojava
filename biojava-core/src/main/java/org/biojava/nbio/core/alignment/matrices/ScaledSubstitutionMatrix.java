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
package org.biojava.nbio.core.alignment.matrices;

import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.template.CompoundSet;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;


/**
 * The biojava-alignment module represents substitution matrices with short
 * values. This is for performance reasons. Some substitution matrices, however,
 * are provided as float values with up to 2 decimal places.
 * <p>
 * In order to be able to use them in the alignment module these are scaled in
 * order to be able to represent as short values.
 * The method {@link #getScale()} provides access to the scaling factor.
 *  
 *  
 * @author Andreas Prlic
 *
 */
public class ScaledSubstitutionMatrix implements
		SubstitutionMatrix<AminoAcidCompound> {
	
    private static final String comment = "#";
    
    private String description, name;
    private short[][] matrix;
    private short max, min;
    private AminoAcidCompoundSet compoundSet;
    
    private List<AminoAcidCompound> rows, cols;
	
    private int scale;
    
    public ScaledSubstitutionMatrix(){
    	compoundSet = AminoAcidCompoundSet.getAminoAcidCompoundSet();
    }
    
    
    
    public int getScale() {
		return scale;
	}

	public void setScale(int scale) {
		this.scale = scale;
	}

	@Override
	public String getDescription() {
		return description;
	}
	@Override
	public void setDescription(String description) {
		this.description = description;
	}
	@Override
	public String getName() {
		return name;
	}
	@Override
	public void setName(String name) {
		this.name = name;
	}
	@Override
	public short[][] getMatrix() {
		return matrix;
	}
	public void setMatrix(short[][] matrix) {
		this.matrix = matrix;
	}
	public short getMax() {
		return max;
	}
	public void setMax(short max) {
		this.max = max;
	}
	public short getMin() {
		return min;
	}
	public void setMin(short min) {
		this.min = min;
	}
	public List<AminoAcidCompound> getRows() {
		return rows;
	}
	public void setRows(List<AminoAcidCompound> rows) {
		this.rows = rows;
	}
	public List<AminoAcidCompound> getCols() {
		return cols;
	}
	public void setCols(List<AminoAcidCompound> cols) {
		this.cols = cols;
	}
	public static String getComment() {
		return comment;
	}
	
	  /**
     * Returns in a format similar to the standard NCBI files.
     */
    @Override
    public String toString() {
    	
    	String newline = System.getProperty("line.separator");
        StringBuilder s = new StringBuilder();
        
        
               
        StringTokenizer st = new StringTokenizer(description, newline);
        while (st.hasMoreTokens()) {
            String line = st.nextToken();
            if (!line.startsWith(comment)) {
                s.append(comment);
            }
            s.append(String.format("%s%n", line));
        }
        
        if ( scale != 1)
        	s.append("# Matrix scaled by a factor of " + scale + newline );
        s.append(getMatrixAsString());
        return s.toString();
    }
	
	
	
	@Override
	public CompoundSet<AminoAcidCompound> getCompoundSet() {
		return compoundSet;
	}
	  @Override
	    public String getMatrixAsString() {
	        StringBuilder s = new StringBuilder();
	        
	        
	        
	        int lengthCompound = compoundSet.getMaxSingleCompoundStringLength(), lengthRest =
	                Math.max(Math.max(Short.toString(min).length(), Short.toString(max).length()), lengthCompound) + 1;
	        
	        String padCompound = "%" + Integer.toString(lengthCompound) + "s",
	                padRest = "%" + Integer.toString(lengthRest);
	        
	        for (int i = 0; i < lengthCompound; i++) {
	            s.append(" ");
	        }
	        for (AminoAcidCompound col : cols) {
	            s.append(String.format(padRest + "s", compoundSet.getStringForCompound(col)));
	        }
	        s.append(String.format("%n"));
	        for (AminoAcidCompound row : rows) {
	            s.append(String.format(padCompound, compoundSet.getStringForCompound(row)));
	            for (AminoAcidCompound col : cols) {
	                s.append(String.format(padRest + "d", getValue(row, col)));
	            }
	            s.append(String.format("%n"));
	        }
	        return s.toString();
	    }
	@Override
	public short getMaxValue() {
		return max;
	}
	@Override
	public short getMinValue() {
		return min;
	}
	@Override
	public short getValue(AminoAcidCompound from, AminoAcidCompound to) {
		 int row = rows.indexOf(from), col = cols.indexOf(to);
	        if (row == -1 || col == -1) {
	            row = cols.indexOf(from);
	            col = rows.indexOf(to);
	            if (row == -1 || col == -1) {
	                return min;
	            }
	        }
	        return matrix[row][col];

		
	}
	
	
	@Override
	public SubstitutionMatrix<AminoAcidCompound> normalizeMatrix(short scale) {
		return null;
	}
    

	@Override
	public Map<AminoAcidCompound, Short> getRow(AminoAcidCompound row) {
		int rowIndex = rows.indexOf(row);
		Map<AminoAcidCompound, Short> map = new HashMap<AminoAcidCompound, Short>();
		for (int colIndex = 0; colIndex < matrix[rowIndex].length; colIndex++) {
			map.put(cols.get(colIndex), matrix[rowIndex][colIndex]);
		}
		return map;
	}

	@Override
	public Map<AminoAcidCompound, Short> getColumn(AminoAcidCompound column) {
		int colIndex = cols.indexOf(column);
		Map<AminoAcidCompound, Short> map = new HashMap<AminoAcidCompound, Short>();
		for (int i = 0; i < matrix.length; i++) {
			map.put(rows.get(i), matrix[i][colIndex]);
		}
		return map;
	}


}
