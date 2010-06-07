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

package org.biojava3.alignment.template;

import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;

public interface SubstitutionMatrix<S extends CompoundSet<C>, C extends Compound> {

    S getCompoundSet();

    String getDescription();

    short[][] getMatrix();

    String getMatrixAsString();

    short getMaxValue();

    short getMinValue();

    String getName();

    short getValue(C from, C to);

    void normalizeMatrix(short scale);

    void setDescription(String description);

    void setName(String name);

}
