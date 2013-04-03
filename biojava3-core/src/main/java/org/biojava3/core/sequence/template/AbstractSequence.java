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

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.TaxonomyID;
import org.biojava3.core.sequence.features.AbstractFeature;
import org.biojava3.core.sequence.features.DatabaseReferenceInterface;
import org.biojava3.core.sequence.features.FeatureInterface;
import org.biojava3.core.sequence.features.FeaturesKeyWordInterface;
import org.biojava3.core.sequence.loader.UniprotProxySequenceReader;
import org.biojava3.core.sequence.location.SequenceLocation;
import org.biojava3.core.sequence.location.SimpleLocation;
import org.biojava3.core.sequence.location.template.Location;
import org.biojava3.core.sequence.storage.ArrayListSequenceReader;

/**
 *
 * The base class for DNA, RNA and Protein sequences.
 * @param <C>
 */
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
    private String source = null;
    private ArrayList<String> notesList = new ArrayList<String>();
    private Double sequenceScore = null;
    private FeaturesKeyWordInterface featuresKeyWord = null;
    private DatabaseReferenceInterface databaseReferences = null;
    private ArrayList<FeatureInterface<AbstractSequence<C>, C>> features =
            new ArrayList<FeatureInterface<AbstractSequence<C>, C>>();
    private LinkedHashMap<String, ArrayList<FeatureInterface<AbstractSequence<C>, C>>> groupedFeatures =
            new LinkedHashMap<String, ArrayList<FeatureInterface<AbstractSequence<C>, C>>>();

    public AbstractSequence() {
    }

    /**
     * Create a Sequence from a simple string where the values should be found in compoundSet
     * @param seqString
     * @param compoundSet
     */
    public AbstractSequence(String seqString, CompoundSet<C> compoundSet) {
        setCompoundSet(compoundSet);
        sequenceStorage = new ArrayListSequenceReader<C>();
        sequenceStorage.setCompoundSet(this.getCompoundSet());
        sequenceStorage.setContents(seqString);
    }

    /**
     * A ProxySequenceReader allows abstraction of both the storage of the sequence data and the location
     * of the sequence data. A variety of use cases are possible. A ProxySequenceReader that knows the offset and of the sequence in
     * a large fasta file. A ProxySequenceReader that can pull Sequence data from Uniprot, NCBI or a custom database.
     * If the ProxySequecneReader implements various interfaces then the sequence will set those interfaces so that calls to
     * various methods will be valid.
     *
     * @param proxyLoader
     * @param compoundSet
     */
    public AbstractSequence(SequenceReader<C> proxyLoader, CompoundSet<C> compoundSet) {
        setCompoundSet(compoundSet);
        setProxySequenceReader(proxyLoader);
    }

    /**
     * Very important method that allows external mappings of sequence data and features. This method
     * will gain additional interface inspection that allows external data sources with knowledge
     * of features for a sequence to be supported. 
     *  
     * @param proxyLoader
     */
    public void setProxySequenceReader(SequenceReader<C> proxyLoader) {
        this.sequenceStorage = proxyLoader;
        if (proxyLoader instanceof FeaturesKeyWordInterface) {
            this.setFeaturesKeyWord((FeaturesKeyWordInterface) sequenceStorage);
        }
        if (proxyLoader instanceof DatabaseReferenceInterface) {
            this.setDatabaseReferences((DatabaseReferenceInterface) sequenceStorage);
        }
        if(getAccession() == null && proxyLoader instanceof UniprotProxySequenceReader){ // we have lots of unsupported operations for this call so quick fix to allow this tow rork
            this.setAccession(proxyLoader.getAccession());
        }
    }

    public SequenceReader<C> getProxySequenceReader() {
        return (SequenceReader<C>) sequenceStorage;
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
     * Provided for convince if the developer needs to associate data with a sequence
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

    /**
     * Added support for the source of this sequence for GFF3 export
     * If a sub sequence doesn't have  source then check for parent source
     * @return the source
     */
    public String getSource() {
        if (source != null) {
            return source;
        }
        if (parentSequence != null) {
            return parentSequence.getSource();
        }
        return null;
    }

    /**
     * Added support for the source of this sequence for GFF3 export
     * @param source the source to set
     */
    public void setSource(String source) {

        this.source = source;
    }

    /**
     * Add notes about this sequence that will get exported for GFF3
     * @param note
     */
    public void addNote(String note) {
        notesList.add(note);
    }

    public void removeNote(String note) {
        notesList.remove(note);
    }

    /**
     * @return the notesList
     */
    public ArrayList<String> getNotesList() {
        return notesList;
    }

    /**
     * @param notesList the notesList to set
     */
    public void setNotesList(ArrayList<String> notesList) {
        this.notesList = notesList;
    }

    /**
     * Provide place holder for a metric that indicate a score associated with the sequence
     * @return the sequenceScore
     */
    public Double getSequenceScore() {
        return sequenceScore;
    }

    /**
     * @param sequenceScore the sequenceScore to set
     */
    public void setSequenceScore(Double sequenceScore) {
        this.sequenceScore = sequenceScore;
    }

    /**
     * Return features at a sequence position by type
     * @param featureType
     * @param bioSequencePosition
     * @return
     */
    public List<FeatureInterface<AbstractSequence<C>, C>> getFeatures(String featureType, int bioSequencePosition) {
        ArrayList<FeatureInterface<AbstractSequence<C>, C>> featureHits =
                new ArrayList<FeatureInterface<AbstractSequence<C>, C>>();
        List<FeatureInterface<AbstractSequence<C>, C>> features = getFeaturesByType(featureType);
        if (features != null) {
            for (FeatureInterface<AbstractSequence<C>, C> feature : features) {
                if (bioSequencePosition >= feature.getLocations().getStart().getPosition() && bioSequencePosition <= feature.getLocations().getEnd().getPosition()) {
                    featureHits.add(feature);
                }
            }
        }
        return featureHits;
    }

    /**
     * Return features at a sequence position
     * @param featureType
     * @param bioSequencePosition
     * @return
     */
    public List<FeatureInterface<AbstractSequence<C>, C>> getFeatures(int bioSequencePosition) {
        ArrayList<FeatureInterface<AbstractSequence<C>, C>> featureHits =
                new ArrayList<FeatureInterface<AbstractSequence<C>, C>>();
        if (features != null) {
            for (FeatureInterface<AbstractSequence<C>, C> feature : features) {
                if (bioSequencePosition >= feature.getLocations().getStart().getPosition() && bioSequencePosition <= feature.getLocations().getEnd().getPosition()) {
                    featureHits.add(feature);
                }
            }
        }
        return featureHits;
    }

    /**
     *
     * @return
     */
    public List<FeatureInterface<AbstractSequence<C>, C>> getFeatures() {
        return features;
    }

    /**
     * Method to help set the proper details for a feature as it relates to a sequence
     * where the feature needs to have a location on the sequence
     * @param bioStart
     * @param bioEnd
     * @param feature
     */
    public void addFeature(int bioStart, int bioEnd, FeatureInterface<AbstractSequence<C>, C> feature) {
        SequenceLocation<AbstractSequence<C>, C> sequenceLocation =
                new SequenceLocation<AbstractSequence<C>, C>(bioStart, bioEnd, this);
        feature.setLocation(sequenceLocation);
        addFeature(feature);
    }

    /**
     * Add a feature to this sequence. The feature will be added to the collection where the order is start position and if more than
     * one feature at the same start position then longest is added first. This helps on doing feature layout for displaying features
     * in SequenceFeaturePanel
     * @param feature
     */
    public void addFeature(FeatureInterface<AbstractSequence<C>, C> feature) {
        features.add(feature);
        ArrayList<FeatureInterface<AbstractSequence<C>, C>> featureList = groupedFeatures.get(feature.getType());
        if (featureList == null) {
            featureList = new ArrayList<FeatureInterface<AbstractSequence<C>, C>>();
            groupedFeatures.put(feature.getType(), featureList);
        }
        featureList.add(feature);
        Collections.sort(features, AbstractFeature.LOCATION_LENGTH);
        Collections.sort(featureList, AbstractFeature.LOCATION_LENGTH);
    }

    /**
     * Remove a feature from the sequence
     * @param feature
     */
    public void removeFeature(FeatureInterface<AbstractSequence<C>, C> feature) {
        features.remove(feature);
        ArrayList<FeatureInterface<AbstractSequence<C>, C>> featureList = groupedFeatures.get(feature.getType());
        if (featureList != null) {
            featureList.remove(feature);
            if (featureList.isEmpty()) {
                groupedFeatures.remove(feature.getType());
            }
        }
    }

    /**
     *
     * @param type
     * @return
     */
    public List<FeatureInterface<AbstractSequence<C>, C>> getFeaturesByType(String type) {
        List<FeatureInterface<AbstractSequence<C>, C>> features = groupedFeatures.get(type);
        if (features == null) {
            features = new ArrayList<FeatureInterface<AbstractSequence<C>, C>>();
        }
        return features;
    }

    /**
     * @return the featuresKeyWord
     */
    public FeaturesKeyWordInterface getFeaturesKeyWord() {
        return featuresKeyWord;
    }

    /**
     * @param featuresKeyWord the featuresKeyWord to set
     */
    public void setFeaturesKeyWord(FeaturesKeyWordInterface featuresKeyWord) {
        this.featuresKeyWord = featuresKeyWord;
    }

    /**
     * @return the databaseReferences
     */
    public DatabaseReferenceInterface getDatabaseReferences() {
        return databaseReferences;
    }

    /**
     * @param databaseReferences the databaseReferences to set
     */
    public void setDatabaseReferences(DatabaseReferenceInterface databaseReferences) {
        this.databaseReferences = databaseReferences;
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
    public String getSequenceAsString(Integer bioStart, Integer bioEnd, Strand strand) {
        Location loc = new SimpleLocation(bioStart, bioEnd, strand);
        return loc.getSubSequence(this).getSequenceAsString();
    }

    /**
     * Default case is to assume strand is positive because only CDSSequence can be either positive or negative Strand.
     * @return
     */
    public String getSequenceAsString() {
        return SequenceMixin.toString(this);

    }

    /**
     *
     * @return
     */
    public List<C> getAsList() {
        return SequenceMixin.toList(this);
    }

    /**
     *
     * @param position The 1-indexed position of the amino acid
     * @return
     */
    public C getCompoundAt(int position) {
        return getSequenceStorage().getCompoundAt(position);
    }

    /**
     *
     * @param compound
     * @return The first index of compound in this sequence (1-based)
     */
    public int getIndexOf(C compound) {
        return getSequenceStorage().getIndexOf(compound);
    }

    /**
     *
     * @param compound
     * @return The last index of compound in this sequence (1-based)
     */
    public int getLastIndexOf(C compound) {
        return getSequenceStorage().getLastIndexOf(compound);
    }

    /**
     *
     * @return
     */
    public int getLength() {
        return getSequenceStorage().getLength();
    }

    /**
     *
     * @param bioStart
     * @param bioEnd
     * @return
     */
    public SequenceView<C> getSubSequence(final Integer bioStart, final Integer bioEnd) {
        return new SequenceProxyView<C>(this, bioStart, bioEnd);
    }

    /**
     *
     * @return
     */
    public Iterator<C> iterator() {
        return getSequenceStorage().iterator();
    }

    /**
     *
     * @param compounds
     * @return
     */
    public int countCompounds(C... compounds) {
        return SequenceMixin.countCompounds(this, compounds);
    }

    /**
     *
     * @return
     */
    @Override
    public SequenceView<C> getInverse() {
        return SequenceMixin.inverse(this);
    }
    
    //TODO needs equals and hashcode
}
