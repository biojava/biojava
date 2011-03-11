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
package org.biojava3.core.sequence.template;

/**
 * The details of a Compound
 * @author Andy Yates
 */
public abstract class AbstractCompound implements Compound {

  private final String base;
  private final String upperedBase;
  private String shortName = null;
  private String longName = null;
  private String description = null;
  private Float molecularWeight = null;

  public AbstractCompound(String base) {
    this.base = base;
    this.upperedBase = base.toUpperCase();
  }

  public String getBase() {
    return base;
  }

  public String getUpperedBase() {
    return upperedBase;
  }

  public String getDescription() {
    return description;
  }

  public void setDescription(String description) {
    this.description = description;
  }

  public String getShortName() {
    return shortName;
  }

  public void setShortName(String shortName) {
    this.shortName = shortName;
  }

  public String getLongName() {
    return longName;
  }

  public void setLongName(String longName) {
    this.longName = longName;
  }

  public Float getMolecularWeight() {
    return molecularWeight;
  }

  public void setMolecularWeight(Float molecularWeight) {
    this.molecularWeight = molecularWeight;
  }

  @Override
  public String toString() {
    return base;
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (!(obj instanceof AbstractCompound)) {
      return false;
    }
    AbstractCompound them = (AbstractCompound) obj;
    return this.base.equals(them.base);
  }

  @Override
  public int hashCode() {
    return this.base.hashCode();
  }

  @Override
  public boolean equalsIgnoreCase(Compound compound) {
    if (compound == null) {
      return false;
    }
    if (!(compound instanceof AbstractCompound)) {
      return false;
    }
    AbstractCompound them = (AbstractCompound) compound;
    return this.base.toString().equalsIgnoreCase(them.base.toString());
  }
}
