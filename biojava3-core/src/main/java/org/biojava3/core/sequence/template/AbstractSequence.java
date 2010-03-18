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
import org.biojava3.core.sequence.TaxonomyID;
import org.biojava3.core.sequence.storage.ArrayListSequenceBackingStore;

public abstract class AbstractSequence<C extends Compound> implements Sequence<C> {

    private TaxonomyID taxonomy;
    private AccessionID accession;
    private SequenceBackingStore<C> backingStore = null;
    private CompoundSet<C> compoundSet;
    private AnnotationType annotationType = AnnotationType.UNKNOWN;
    private String description;
    private String originalHeader;
    private Collection<Object> userCollection;


    public AbstractSequence() {
    }

    public AbstractSequence(String seqString, CompoundSet<C> compoundSet) {
        setCompoundSet(compoundSet);
        backingStore = new ArrayListSequenceBackingStore<C>();
        backingStore.setCompoundSet(this.getCompoundSet());
        backingStore.setContents(seqString);
    }

    public AbstractSequence(SequenceProxyLoader<C> proxyLoader, CompoundSet<C> compoundSet) {
        setCompoundSet(compoundSet);
        this.backingStore = proxyLoader;
    }

    /**
     * Provided for convience if the developer needs to associate data with a sequence
     *
     * @return
     */

    public Collection<Object> getUserCollection(){

        return userCollection;
    }

    /**
     *
     * @param userCollection
     */

    public void setUserCollection(Collection<Object> userCollection){
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
        return compoundSet;
    }

    public void setCompoundSet(CompoundSet<C> compoundSet) {
        this.compoundSet = compoundSet;
    }

    @Override
    public String toString() {
        return getSequenceAsString();
    }

    public String getSequenceAsString() {
        return backingStore.getSequenceAsString();
    }

    public List<C> getAsList() {
        return backingStore.getAsList();
    }

    public C getCompoundAt(int position) {
        return backingStore.getCompoundAt(position);
    }

    public int getIndexOf(C compound) {
        return backingStore.getIndexOf(compound);
    }

    public int getLastIndexOf(C compound) {
        return backingStore.getLastIndexOf(compound);
    }

    public int getLength() {
        return backingStore.getLength();
    }

    public SequenceView<C> getSubSequence(final int start, final int end) {
        return new AbstractSequenceView<C>() {

            public int getEnd() {
                return end;
            }

            public int getStart() {
                return start;
            }

            public Sequence<C> getViewedSequence() {
                return AbstractSequence.this;
            }

            public String getSequenceAsString() {
              return getViewedSequence().getSequenceAsString().substring(start-1, end);
            }

            @Override
            public AccessionID getAccession() {
                return getViewedSequence().getAccession();
            }


        };
    }

    public Iterator<C> iterator() {
        return backingStore.iterator();
    }
}
