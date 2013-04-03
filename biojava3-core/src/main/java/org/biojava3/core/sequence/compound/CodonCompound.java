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
 * Created on 01-21-2010
 */
package org.biojava3.core.sequence.compound;

import org.biojava3.core.sequence.template.AbstractCompound;
import org.biojava3.core.sequence.template.Compound;
/**
 * Define a codon
 * @author Andy Yates
 */
public class CodonCompound extends AbstractCompound {

  private final NucleotideCompound one;
  private final NucleotideCompound two;
  private final NucleotideCompound three;
  private final boolean start;

  public CodonCompound(NucleotideCompound one, NucleotideCompound two,
      NucleotideCompound three, boolean start) {
    super(one.toString()+two.toString()+three.toString());
    this.one = one;
    this.two = two;
    this.three = three;
    this.start = start;
  }

  public boolean equalsIgnoreCase(Compound compound) {
    if (compound == null) {
      return false;
    }
    if (!(compound instanceof CodonCompound)) {
      return false;
    }
    CodonCompound them = (CodonCompound) compound;
    return toString().equalsIgnoreCase(them.toString());
  }

  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (!(obj instanceof CodonCompound)) {
      return false;
    }
    CodonCompound them = (CodonCompound)obj;
    return toString().equals(them.toString());
  }

  public int hashCode() {
    return toString().hashCode();
  }

  public NucleotideCompound getOne() {
    return one;
  }

  public NucleotideCompound getTwo() {
    return two;
  }

  public NucleotideCompound getThree() {
    return three;
  }

  public boolean isStart() {
    return start;
  }

  public String getDescription() {
    // TODO Auto-generated method stub
    return null;
  }

  public String getLongName() {
    // TODO Auto-generated method stub
    return null;
  }

  public Float getMolecularWeight() {
    // TODO Auto-generated method stub
    return null;
  }

  public String getShortName() {
    // TODO Auto-generated method stub
    return null;
  }

  public void setDescription(String description) {
    // TODO Auto-generated method stub

  }

  public void setLongName(String longName) {
    // TODO Auto-generated method stub

  }

  public void setMolecularWeight(Float molecularWeight) {
    // TODO Auto-generated method stub

  }

  public void setShortName(String shortName) {
    // TODO Auto-generated method stub

  }

}
