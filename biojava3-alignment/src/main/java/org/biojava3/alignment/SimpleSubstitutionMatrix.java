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
 * Created on June 9, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import java.util.StringTokenizer;

import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;

/**
 * Implements a data structure which holds the score (penalty or bonus) given during alignment for the exchange of one
 * {@link Compound} in a sequence for another.
 *
 * @author Mark Chapman
 * @param <C> each element of the matrix corresponds to a pair of {@link Compound}s of type C
 */
public class SimpleSubstitutionMatrix<C extends Compound> implements SubstitutionMatrix<C> {

    private static final String comment = "#";

    private CompoundSet<C> compoundSet;
    private String description, name;
    private short[][] matrix;
    private short max, min;
    private List<C> rows, cols;

    /**
     * Creates a substitution matrix using the defaults (BLOSUM 62).
     *
     * @throws ClassCastException if type parameter C is not {@link AminoAcidCompound}
     */
    @SuppressWarnings("unchecked") // TODO proper generic type checking instead of possible ClassCastException
    public SimpleSubstitutionMatrix() {
        this((CompoundSet<C>) AminoAcidCompoundSet.getAminoAcidCompoundSet(), new InputStreamReader(
                SimpleSubstitutionMatrix.class.getResourceAsStream("/blosum62.txt")), "blosum62");
    }

    /**
     * Creates a substitution matrix by reading in a file.
     *
     * @param compoundSet the {@link CompoundSet} on which the matrix is defined
     * @param fileInput file parsed for a substitution matrix
     * @throws FileNotFoundException if fileInput parameter cannot be read
     */
    public SimpleSubstitutionMatrix(CompoundSet<C> compoundSet, File fileInput) throws FileNotFoundException {
        this(compoundSet, new BufferedReader(new FileReader(fileInput)), fileInput.getName());
    }

    /**
     * Creates a substitution matrix by parsing some input.
     *
     * @param compoundSet the {@link CompoundSet} on which the matrix is defined
     * @param input input parsed for a substitution matrix
     * @param name the name (short description) of this matrix
     */
    public SimpleSubstitutionMatrix(CompoundSet<C> compoundSet, Reader input, String name) {
        this(compoundSet, new Scanner(input), name);
    }

    /**
     * Creates a substitution matrix by parsing a String.
     *
     * @param compoundSet the {@link CompoundSet} on which the matrix is defined
     * @param matrixInput String parsed for a substitution matrix
     * @param name the name (short description) of this matrix
     */
    public SimpleSubstitutionMatrix(CompoundSet<C> compoundSet, String matrixInput, String name) {
        this(compoundSet, new Scanner(matrixInput), name);
    }

    /**
     * Creates an identity substitution matrix from match and replace values.
     *
     * @param compoundSet the {@link CompoundSet} on which the matrix is defined
     * @param match matrix value used for equivalent {@link Compound}s
     * @param replace matrix value used for differing {@link Compound}s
     */
    public SimpleSubstitutionMatrix(CompoundSet<C> compoundSet, short match, short replace) {
        this.compoundSet = compoundSet;
        description = "Identity matrix. All replaces and all matches are treated equally.";
        name = "IDENTITY_" + match + "_" + replace;
        max = (match > replace) ? match : replace;
        min = (match < replace) ? match : replace;
        rows = cols = compoundSet.getAllCompounds();
        matrix = new short[rows.size()][cols.size()];
        for (int r = 0; r < rows.size(); r++) {
            for (int c = 0; c < cols.size(); c++) {
                try {
                    matrix[r][c] = (compoundSet.compoundsEquivalent(rows.get(r), cols.get(c))) ? match : replace;
                } catch (UnsupportedOperationException e) {
                    matrix[r][c] = (r == c) ? match : replace;
                }
            }
        }
    }

    // helper constructor that creates a substitution matrix by parsing input
    private SimpleSubstitutionMatrix(CompoundSet<C> compoundSet, Scanner input, String name) {
        this.compoundSet = compoundSet;
        this.name = name;
        max = Short.MIN_VALUE;
        min = Short.MAX_VALUE;
        rows = new ArrayList<C>();
        cols = new ArrayList<C>();
        StringBuilder descriptionIn = new StringBuilder();
        List<short[]> matrixIn = new ArrayList<short[]>();
        while(input.hasNextLine()) {
            String line = input.nextLine();
            if (line.startsWith(comment)) {
                descriptionIn.append(String.format("%s%n", line));
            } else if (!line.trim().isEmpty()) {
                StringTokenizer st = new StringTokenizer(line);
                if (cols.isEmpty()) {
                    while (st.hasMoreTokens()) {
                        cols.add(compoundSet.getCompoundForString(st.nextToken()));
                    }
                } else {
                    rows.add(compoundSet.getCompoundForString(st.nextToken()));
                    short[] row = new short[cols.size()];
                    for (int i = 0; i < row.length && st.hasMoreTokens(); i++) {
                        row[i] = Short.parseShort(st.nextToken());
                        max = (max > row[i]) ? max : row[i];
                        min = (min < row[i]) ? min : row[i];
                    }
                    matrixIn.add(row);
                }
            }
        }
        input.close();
        description = descriptionIn.toString();
        matrix = new short[rows.size()][cols.size()];
        for (int i = 0; i < rows.size(); i++) {
            matrix[i] = matrixIn.get(i);
        }
    }

    @Override
    public CompoundSet<C> getCompoundSet() {
        return compoundSet;
    }

    @Override
    public String getDescription() {
        return description;
    }

    @Override
    public short[][] getMatrix() {
        short[][] copy = new short[matrix.length][matrix[0].length];
        for (int i = 0; i < copy.length; i++) {
            copy[i] = Arrays.copyOf(matrix[i], matrix[i].length);
        }
        return copy;
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
        for (C col : cols) {
            s.append(String.format(padRest + "s", compoundSet.getStringForCompound(col)));
        }
        s.append(String.format("%n"));
        for (C row : rows) {
            s.append(String.format(padCompound, compoundSet.getStringForCompound(row)));
            for (C col : cols) {
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
    public String getName() {
        return name;
    }

    @Override
    public short getValue(C from, C to) {
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
    public SubstitutionMatrix<C> normalizeMatrix(short scale) {
        // TODO SubstitutionMatrix<C> normalizeMatrix(short)
        return null;
    }

    @Override
    public void setDescription(String description) {
        this.description = description;
    }

    @Override
    public void setName(String name) {
        this.name = name;
    }

    /**
     * Returns in a format similar to the standard NCBI files.
     */
    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        StringTokenizer st = new StringTokenizer(description, "\n\r");
        while (st.hasMoreTokens()) {
            String line = st.nextToken();
            if (!line.startsWith(comment)) {
                s.append(comment);
            }
            s.append(String.format("%s%n", line));
        }
        s.append(getMatrixAsString());
        return s.toString();
    }

}
