package org.biojava3.alignment.aaindex;

import java.util.List;
import java.util.StringTokenizer;

import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.template.CompoundSet;


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

	public String getDescription() {
		return description;
	}
	public void setDescription(String description) {
		this.description = description;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
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
    
    

}
