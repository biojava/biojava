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
 */
package org.biojava.bio.seq.io.agave;
import org.biojava.bio.Annotation;

/**
 * <p>This interface defines mapping from BioJava into AGAVE format.
 * As data from different sources is stored differently in BioJava
 * it is impossible to define universal mapping from BioJava to Agave.
 * Currently I implemented two mappings:</p>
 *
 * <p>( embl  ->  )biojava -> agave</p>
 * <p>( agave ->  )biojava -> agave</p>
 *
 * @author Hanning Ni     Doubletwist Inc
 */
public interface AGAVEAnnotFilter
{
        public static final int FORWARD = 1 ;
        public static final int COMPLEMENT = -1 ;
        public static final int BOTH_FORWARD_COMPLEMENT = 0 ;

        public String getAccession(Annotation annot);
        public String getLabel(Annotation annot);
        public String getElementId(Annotation annot);
        public String getSequenceId(Annotation annot);
        public String getKeyword(Annotation annot);
        public String getOrganism(Annotation annot);
        public String getDescription(Annotation annot);
        public String getNote (Annotation annot);
        public String getVersion(Annotation annot);
        public String getOS(Annotation annot);
        public String getMolType(Annotation annot);
        public String getTaxonId(Annotation annot);
        public String getCloneId(Annotation annot);
        public String getCloneLibrary(Annotation annot);
        public String getChromosome(Annotation annot);
        public String getMapPosition(Annotation annot);
        public String getEcNumber(Annotation annot);
        public String getCreateDate(Annotation annot);
        public String getUpdateDate(Annotation annot);
        public AGAVEDbId[]  getAltIds(Annotation annot);
        public AGAVEXrefs[]  getXrefs(Annotation annot);
        public AGAVEMapLocation[]  getMapLocation(Annotation annot);
        public AGAVERelatedAnnot[]  getRelatedAnnot(Annotation annot);
        public String[] getElementIds(Annotation annot) ;
        public String[] getExonIds(Annotation annot) ;
        public String getChromNum(Annotation annot);
        public AGAVEProperty[] getProperty(Annotation annot, String type);
        public AGAVEDbId  getDbId(Annotation annot) ;
        public String getGroupOrder(Annotation annot) ;
        public String getFeatureType(Annotation annot) ;
        public String getResultType(Annotation annot) ;
        public String getConfidence(Annotation annot) ;
        public String getAlignLength(Annotation annot) ;
        public String getAlignUnits(Annotation annot) ;
        public String getMatchDesc(Annotation annot) ;
        public String getMatchAlign(Annotation annot) ;
        public AGAVEQueryRegion getQueryRegion(Annotation annot );
        public AGAVEMatchRegion getMatchRegion(Annotation annot );
        public AGAVEIdAlias[] getIdAlias(Annotation annot);
        public String getClassifySystem(Annotation annot);
        public String getClassifyId(Annotation annot);
        public String getClassifyType(Annotation annot);
}
