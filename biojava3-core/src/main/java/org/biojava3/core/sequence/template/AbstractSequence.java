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
 *
 * @author Richard Holland
 * @auther Scooter Willis
 *
 */
package org.biojava3.core.sequence.template;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.TaxonomyID;

import org.biojava3.core.sequence.storage.ArrayListSequenceReader;

public abstract class AbstractSequence<C extends Compound> implements Sequence<C> {

    private TaxonomyID taxonomy;
    private AccessionID accession;
    private SequenceReader<C> sequenceStorage = null;
    private CompoundSet<C> compoundSet;
    private AnnotationType annotationType = AnnotationType.UNKNOWN;
    private String description;
    private String originalHeader;
    private Collection<Object> userCollection;
    private Integer bioBegin = null;
    private Integer bioEnd = null;
    private AbstractSequence<C> parentSequence = null;

    public AbstractSequence() {
    }

    public AbstractSequence(String seqString, CompoundSet<C> compoundSet) {
        setCompoundSet(compoundSet);
        sequenceStorage = new ArrayListSequenceReader<C>();
        sequenceStorage.setCompoundSet(this.getCompoundSet());
        sequenceStorage.setContents(seqString);
    }

    public AbstractSequence(ProxySequenceReader<C> proxyLoader, CompoundSet<C> compoundSet) {
        setCompoundSet(compoundSet);
        this.sequenceStorage = proxyLoader;
    }

    /**
     * @return the bioBegin
     */
    public Integer getBioBegin() {
        if (bioBegin == null) {
            return 1;
        } else {
            return bioBegin;
        }
    }

    /**
     * @param bioBegin the bioBegin to set
     */
    public void setBioBegin(Integer begin) {
        this.bioBegin = begin;
    }

    /**
     * @return the bioEnd
     */
    public Integer getBioEnd() {
        if (bioEnd == null) {
            return this.getLength();
        } else {
            return bioEnd;
        }
    }

    /**
     * @param bioEnd the bioEnd to set
     */
    public void setBioEnd(Integer end) {
        this.bioEnd = end;
    }

    /**
     * Provided for convience if the developer needs to associate data with a sequence
     *
     * @return
     */
    public Collection<Object> getUserCollection() {

        return userCollection;
    }

    /**
     *
     * @param userCollection
     */
    public void setUserCollection(Collection<Object> userCollection) {
        this.userCollection = userCollection;
    }

    /**
     * @return the annotation
     */
    public AnnotationType getAnnotationType() {
        return annotationType;
    }

    /**
     * @param annotation the annotation to set
     */
    public void setAnnotationType(AnnotationType annotationType) {
        this.annotationType = annotationType;
    }

    /**
     * @return the description
     */
    public String getDescription() {
        return description;
    }

    /**
     * @param description the description to set
     */
    public void setDescription(String description) {
        this.description = description;
    }

    /**
     * @return the originalHeader
     */
    public String getOriginalHeader() {
        return originalHeader;
    }

    /**
     * @param originalHeader the originalHeader to set
     */
    public void setOriginalHeader(String originalHeader) {
        this.originalHeader = originalHeader;
    }

    /**
     * @return the parentSequence
     */
    public AbstractSequence<C> getParentSequence() {
        return parentSequence;
    }

    /**
     * @param parentSequence the parentSequence to set
     */
    public void setParentSequence(AbstractSequence<C> parentSequence) {
        this.parentSequence = parentSequence;
    }

    public enum AnnotationType {

        CURATED, PREDICTED, UNKNOWN;
    }

    /**
     * @return the accession
     */
    public AccessionID getAccession() {
        return accession;
    }

    /**
     * @param accession the accession to set
     */
    public void setAccession(AccessionID accession) {
        this.accession = accession;
    }

    /**
     * @return the species
     */
    public TaxonomyID getTaxonomy() {
        return taxonomy;
    }

    /**
     * @param species the species to set
     */
    public void setTaxonomy(TaxonomyID taxonomy) {
        this.taxonomy = taxonomy;
    }

    public CompoundSet<C> getCompoundSet() {
        if (compoundSet != null) {
            return compoundSet;
        }
        if (parentSequence != null) {
            return parentSequence.getCompoundSet();
        }
        return null;


    }

    public void setCompoundSet(CompoundSet<C> compoundSet) {
        this.compoundSet = compoundSet;
    }

    @Override
    public String toString() {
        return getSequenceAsString();
    }

    private SequenceReader<C> getSequenceStorage() {
        if (sequenceStorage != null) {
            return sequenceStorage;
        }
        if (parentSequence != null) {
            return parentSequence.getSequenceStorage();
        }
        return null;
    }

    /**
     *
     * @param begin
     * @param end
     * @param strand 
     * @return
     */
    public String getSequenceAsString(Integer begin, Integer end, Strand strand) {
        SequenceReader<C> ss = getSequenceStorage();
        return ss.getSequenceAsString(begin, end, strand);
    }

    /**
     * Default case is to assume strand is positive because only CDSSequence can be either positive or negative Strand.
     * @return
     */
    public String getSequenceAsString() {
        return getSequenceAsString(this.getBioBegin(), this.getBioEnd(), Strand.POSITIVE);
    }

    public List<C> getAsList() {
        return getSequenceStorage().getAsList();
    }

    public C getCompoundAt(int position) {
        return getSequenceStorage().getCompoundAt(position);
    }

    public int getIndexOf(C compound) {
        return getSequenceStorage().getIndexOf(compound);
    }

    public int getLastIndexOf(C compound) {
        return getSequenceStorage().getLastIndexOf(compound);
    }

    public int getLength() {
        return getSequenceStorage().getLength();
    }

    public SequenceView<C> getSubSequence(final Integer bioStart, final Integer bioEnd) {
        return new SequenceProxyView<C>(AbstractSequence.this, bioStart, bioEnd);
    }

    public Iterator<C> iterator() {
        return getSequenceStorage().iterator();
    }

    public int countCompounds(C... compounds) {
      return SequenceMixin.countCompounds(this, compounds);
    }
}
