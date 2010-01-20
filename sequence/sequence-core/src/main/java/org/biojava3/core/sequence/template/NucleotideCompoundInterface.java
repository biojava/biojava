/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.template;

/**
 *
 * @author Scooter
 */
public interface NucleotideCompoundInterface<C extends Compound> extends Compound{

    public void setComplement(NucleotideCompoundInterface nucleotideCompound);
    public C getComplement();
}
