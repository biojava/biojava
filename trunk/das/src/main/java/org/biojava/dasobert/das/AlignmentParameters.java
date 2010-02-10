/*
 *                  BioJava development code
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
 * Created on Dec 10, 2005
 *
 */
package org.biojava.dasobert.das;

import org.biojava.dasobert.dasregistry.Das1Source;
import org.biojava.dasobert.dasregistry.DasCoordinateSystem;

/** a class that stores the arguments that can be sent to the AlignmentThread class
 * 
 * @author Andreas Prlic
 *
 */
public class AlignmentParameters {

    String query;
    String subject;
    String queryPDBChainId;
    String subjectPDBChainId;
    
    
    DasCoordinateSystem queryCoordinateSystem;
    DasCoordinateSystem subjectCoordinateSystem;
    Das1Source[] dasSources;
    
    
    public static String DEFAULT_PDBCOORDSYS     = "PDBresnum,Protein Structure";
    public static String DEFAULT_UNIPROTCOORDSYS = "UniProt,Protein Sequence";
    public static String DEFAULT_ENSPCOORDSYS    = "Ensembl,Protein Sequence";
    
    
    public AlignmentParameters() {
        super();
        dasSources = new Das1Source[0];

    }

    public DasCoordinateSystem getDefaultPDBCoordSys(){
        return DasCoordinateSystem.fromString(DEFAULT_PDBCOORDSYS);
    }
    public DasCoordinateSystem getDefaultUniProtCoordSys(){
        return DasCoordinateSystem.fromString(DEFAULT_UNIPROTCOORDSYS);
    }
    public DasCoordinateSystem getDefaultEnspCoordSys(){
        return DasCoordinateSystem.fromString(DEFAULT_ENSPCOORDSYS);
    }
    

    public Das1Source[] getDasSources() {
        return dasSources;
    }



    public void setDasSources(Das1Source[] dasSources) {
        this.dasSources = dasSources;
    }
    
    public void setDasSource(Das1Source dasSource){
    	
    	dasSources = new Das1Source[] {dasSource};
    	
    }



    public String getQuery() {
        return query;
    }



    public void setQuery(String query) {
        this.query = query;
    }



    public DasCoordinateSystem getQueryCoordinateSystem() {
        return queryCoordinateSystem;
    }



    public void setQueryCoordinateSystem(DasCoordinateSystem queryCoordinateSystem) {
        this.queryCoordinateSystem = queryCoordinateSystem;
    }



    public String getQueryPDBChainId() {
        return queryPDBChainId;
    }



    public void setQueryPDBChainId(String queryPDBChainId) {
        this.queryPDBChainId = queryPDBChainId;
    }



    public String getSubject() {
        return subject;
    }



    public void setSubject(String subject) {
        this.subject = subject;
    }



    public DasCoordinateSystem getSubjectCoordinateSystem() {
        return subjectCoordinateSystem;
    }



    public void setSubjectCoordinateSystem(
            DasCoordinateSystem subjectCoordinateSystem) {
        this.subjectCoordinateSystem = subjectCoordinateSystem;
    }



    public String getSubjectPDBChainId() {
        return subjectPDBChainId;
    }



    public void setSubjectPDBChainId(String subjectPDBChainId) {
        this.subjectPDBChainId = subjectPDBChainId;
    }
    
    
    

}
