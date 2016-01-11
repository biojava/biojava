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
 * Created on June 7, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.core.alignment.template;

import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;

import java.util.Map;

/**
 * Defines a data structure which holds the score (penalty or bonus) given during alignment for the exchange of one
 * {@link Compound} in a sequence for another.
 *
 * @author Mark Chapman
 * @author Paolo Pavan
 * @param <C> each element of the matrix corresponds to a pair of {@link Compound}s of type C
 */
public interface SubstitutionMatrix<C extends Compound> {

    /**
     * Returns the {@link CompoundSet} on which the matrix is defined.
     *
     * @return the {@link CompoundSet} on which the matrix is defined
     */
    CompoundSet<C> getCompoundSet();

    /**
     * Returns the description of this matrix.
     *
     * @return description
     */
    String getDescription();

    /**
     * Returns entire matrix.
     *
     * @return matrix
     */
    short[][] getMatrix();

    /**
     * Returns this matrix as a formatted String with {@link Compound} labels along the axes.
     *
     * @return this matrix as a formatted String
     */
    String getMatrixAsString();

    /**
     * Returns the maximum value in this matrix.
     *
     * @return the maximum value in this matrix
     */
    short getMaxValue();

    /**
     * Returns the minimum value in this matrix.
     *
     * @return the minimum value in this matrix
     */
    short getMinValue();

    /**
     * Returns the name (short description) of this matrix.
     *
     * @return name
     */
    String getName();

    /**
     * Returns value in matrix for conversion from first {@link Compound} to the second.  If an argument does not
     * belong to the {@link CompoundSet}, this could either throw an {@link IllegalArgumentException} or it could
     * return {@link #getMinValue()}.
     *
     * @param from original {@link Compound}
     * @param to replacement {@link Compound}
     * @return value in matrix for conversion from first {@link Compound} to the second
     * @throws IllegalArgumentException possibly, if an argument does not belong to the {@link CompoundSet}
     */
    short getValue(C from, C to);

    /**
     * Rescales the matrix so that to {@link #getMaxValue()} - {@link #getMinValue()} = scale.
     *
     * @param scale new normalization scale of this matrix
     * @throws IllegalArgumentException if scale < 1
     */
    SubstitutionMatrix<C> normalizeMatrix(short scale);

    /**
     * Sets the description of this matrix.
     *
     * @param description new description
     */
    void setDescription(String description);

    /**
     * Sets the name (short description) of this matrix.
     *
     * @param name new name
     */
    void setName(String name);

    Map<C, Short> getRow(C row);
    
    Map<C, Short> getColumn(C column);
    
}
