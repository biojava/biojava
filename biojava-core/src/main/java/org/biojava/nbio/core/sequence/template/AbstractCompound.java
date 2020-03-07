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
package org.biojava.nbio.core.sequence.template;


/**
 * The details of a Compound
 *
 * @author Andy Yates
 */
public abstract class AbstractCompound implements Compound {

    private final String base;
    private final String upperedBase;
    private String shortName = null;
    private String longName = null;
    private String description = null;
    protected float molecularWeight = Float.NaN;
    final int hash;

    protected AbstractCompound(String base) {
        this.base = base;
        this.upperedBase = base.toUpperCase();
        this.hash = base.hashCode();
    }

    public String getBase() {
        return base;
    }

    public String getUpperedBase() {
        return upperedBase;
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
    public String getShortName() {
        return shortName;
    }

    @Override
    public void setShortName(String shortName) {
        this.shortName = shortName;
    }

    @Override
    public String getLongName() {
        return longName;
    }

    @Override
    public void setLongName(String longName) {
        this.longName = longName;
    }

    @Override
    public final float getMolecularWeight() {
        return molecularWeight;
    }

    @Override
    public final String toString() {
        return base;
    }

    @Override
    public final boolean equals(Object obj) {
        return obj == this || (obj instanceof AbstractCompound && this.base.equals(((AbstractCompound) obj).base));
    }

    @Override
    public final int hashCode() {
        return hash;
    }

    @Override
    public final boolean equalsIgnoreCase(Compound compound) {
        if (this == compound) return true;
        if (!(compound instanceof AbstractCompound)) {
            return false;
        }
        AbstractCompound them = (AbstractCompound) compound;
        return this.base.equalsIgnoreCase(them.base);
    }
}
