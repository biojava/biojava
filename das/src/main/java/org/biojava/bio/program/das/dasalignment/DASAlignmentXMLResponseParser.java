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
 * Created on 13.5.2004
 * @author Andreas Prlic
 *
 */

package org.biojava.bio.program.das.dasalignment ;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.biojava.bio.Annotation;
import org.biojava.bio.program.ssbind.AnnotationFactory;
import org.biojava.bio.structure.AtomImpl;
import org.xml.sax.Attributes;
import org.xml.sax.helpers.DefaultHandler;

/** A class to Parse the XML response of a DAS Alignment service. 
 * returns an Alignment object.
 *
 * @author Andreas Prlic
 * @since 1.4
 */
public class DASAlignmentXMLResponseParser  extends DefaultHandler{
    ArrayList alignments ;
    Alignment alignment      ;
    HashMap current_object   ;
    String  current_position ;
    HashMap current_block  ;
    HashMap  current_segment ;
    ArrayList segments ;
    
    ArrayList vectors ;
    ArrayList matrices ;
    String current_detailProperty;
    String intObjectId ;
    String cigar ;
    String txt;
    
    public DASAlignmentXMLResponseParser() {
        super() ;
        //System.out.println("in init DASAlignmentXMLResponseParser");
        
        current_position = "start"   ;
        cigar                  = ""  ;
        txt                    = ""  ;
        current_detailProperty = ""  ;
        alignment  = new Alignment() ;
        
        alignments = new ArrayList() ;
        segments   = new ArrayList() ;
        vectors    = new ArrayList() ;
        matrices   = new ArrayList() ;
        
        
    }
    
    /**
     * Returns the alignments.
     *
     * @return an array of Alignment objects 
     
     */
    public Alignment[] getAlignments() {
        
        return (Alignment[])alignments.toArray(new Alignment[alignments.size()]);
    }
    
    /**
     * Returns Alignment at position ...
     *
     * @param position  an int
     * @return an Alignment object
     */
    public Alignment getAlignment(int position) {
        Alignment ra = (Alignment) alignments.get(position) ; 
        return ra ;
    }
    
    public void startElement (String uri, String name, String qName, Attributes atts){
        //System.out.println("startElement " + qName) ;
        if (qName.equals("alignObject")) OBJECThandler     (atts);
        if (qName.equals("sequence")   ) SEQUENCEhandler   (atts);
        if (qName.equals("score")      ) SCOREhandler      (atts);
        if (qName.equals("block")      ) BLOCKhandler      (atts);
        if (qName.equals("segment")    ) SEGMENThandler    (atts);
        if (qName.equals("cigar")      ) CIGARhandler      (atts);
        if (qName.equals("geo3D")      ) GEO3Dhandler     (atts);
        if (qName.equals("vector")     ) VECTORhandler     (atts);
        if (qName.equals("matrix")     ) MATRIXhandler     (atts);
        if (qName.equals("alignObjectDetail")) DETAILhandler(atts);
        
        
    }
    
    public void endElement (String uri, String name, String qName){
        //System.out.println("endElement >" + qName + "< >" + name + "<") ;
        if (qName.equals("geo3D")) {	
            intObjectId = null ;
        }
        
        if (qName.equals("alignObject")) {
            try {
                Annotation oba = AnnotationFactory.makeAnnotation(current_object) ;
                alignment.addObject(oba);
            } catch ( DASException  e) {
                e.printStackTrace() ;
            }
            current_object = new HashMap() ;
            ArrayList details = new ArrayList();
            current_object.put("details",details);
        }	
        if (qName.equals("segment")) {
            Annotation sega = AnnotationFactory.makeAnnotation(current_segment);
            //current_block.add(current_segment);
            segments.add(sega) ;
            current_segment = new HashMap() ;
            
        }
        if (qName.equals("block")) {
            try {
                current_block.put("segments",segments);
                Annotation bloa = AnnotationFactory.makeAnnotation(current_block);
                alignment.addBlock(bloa);
            } catch ( DASException  e) {
                e.printStackTrace() ;
            }
            current_block = new HashMap() ;
            segments = new ArrayList();
        }
        if (qName.equals("alignment")){
            alignments.add(alignment) ;
            
            alignment = new Alignment() ;
        }
        if (qName.equals("cigar")) {
            cigar = cigar.trim();
            current_segment.put("cigar",cigar);
            cigar = "" ;
        }
        if (qName.equals("alignObjectDetail")) {
            String detailValue = txt;
            HashMap detail = new HashMap();
            detail.put("property", current_detailProperty);
            detail.put("detail",  detailValue);
            
            Annotation a = AnnotationFactory.makeAnnotation(detail);
            List details = (List) current_object.get("details");
            details.add(a);
            current_object.put("details",details);
            //System.out.println("adding new detail " + detail);
        }
    }
    
    private void DETAILhandler(Attributes atts){
        //. TODO: deal with dbSource property
        //System.out.println("SUPPORTING DETAIL!!!!!!!!!!!!");
        current_detailProperty = atts.getValue("property");
    }
    
    private void SEGMENThandler(Attributes atts) {
        current_position = "segment";
        current_segment  = new HashMap() ;
        
        String id     = atts.getValue("intObjectId");
        String start  = atts.getValue("start");
        String end    = atts.getValue("end");
        // orientation, not implemented yet ...
        current_segment.put("intObjectId",id);
        if ( start != null ) {
            current_segment.put("start",start);
        }
        if ( end != null ) {
            current_segment.put("end",end) ;
        }
   	}
    
    
    private void CIGARhandler(Attributes atts) {
        current_position = "cigar" ;
        cigar = "" ;
    }
    private void BLOCKhandler(Attributes atts) {
        current_block = new HashMap();
        String blockOrder = atts.getValue("blockOrder");
        
        current_block.put("blockOrder",blockOrder);
        
        try {
            String blockScore = atts.getValue("blockScore");
            if ( blockScore != null ) {
                current_block.put("blockScore",blockScore);
            }
        } catch (Exception e) {} ;
        
    }
    
    
    private void SEQUENCEhandler(Attributes atts) {
        //System.out.println("sequence");
        current_position = "sequence" ;
        String start    = atts.getValue("start");
        String end      = atts.getValue("end");
        if ( start != null ) {
            current_object.put("seqStart",start);
        }
        if ( end != null ) {
            current_object.put("seqEnd",end);
        }
        
    }
    
    
    private void SCOREhandler(Attributes atts) {
        System.out.println("SCOREhandler not implemented,yet...");
        
        HashMap score = new HashMap() ;
        
        String metnam = atts.getValue("methodName");
        String value  = atts.getValue("value");
        score.put("methodName",metnam);
        score.put("value",value);
        Annotation s = AnnotationFactory.makeAnnotation(score);
        try {
            alignment.addScore(s);
        } catch (DASException e ) {
            e.printStackTrace();
            return ;
        }
    }
    
    private void GEO3Dhandler(Attributes atts) {
        intObjectId      = atts.getValue("intObjectId");
    }
    
    private void VECTORhandler(Attributes atts) {
        //System.out.println("VECTORhandler");
        //String intObjectId      = atts.getValue("intObjectId");
        String xs = atts.getValue("x");
        String ys = atts.getValue("y");
        String zs = atts.getValue("z");
        
        AtomImpl atom = new AtomImpl();
        atom.setX(Double.parseDouble(xs));
        atom.setY(Double.parseDouble(ys));
        atom.setZ(Double.parseDouble(zs));
        HashMap vector = new HashMap() ;
        vector.put("intObjectId"   ,intObjectId);
        vector.put("vector",atom);
        
        try {
            Annotation oba = AnnotationFactory.makeAnnotation(vector) ;
            alignment.addVector(oba);
        } catch ( DASException  e) {
            e.printStackTrace() ;
        }	
    }
    
    private void MATRIXhandler(Attributes atts) {
        HashMap matrix = new HashMap() ;
        //String intObjectId      = atts.getValue("intObjectId");
        matrix.put("intObjectId"  ,intObjectId);
        
        for ( int x=1; x<=3; x++){
            for ( int y=1; y<=3; y++){
                String mat = "mat"+x+y;
                String val = atts.getValue(mat);
                //System.out.println(mat+" " + val);
                matrix.put(mat,val) ;
            }
        }
        
        try {
            Annotation oba = AnnotationFactory.makeAnnotation(matrix) ;
            alignment.addMatrix(oba);
        } catch ( DASException  e) {
            e.printStackTrace() ;
        }
    }
    
    
    private void OBJECThandler(Attributes atts) {
        // found a new object
        String dbAccessionId    = atts.getValue("dbAccessionId");
        String objectVersion    = atts.getValue("objectVersion");
        String intObjectId      = atts.getValue("intObjectId");
        String type             = "" ;
        
        try { type = atts.getValue("type");} catch (Exception e) {} 
        
        
        String dbSource         = atts.getValue("dbSource");
        String dbVersion        = atts.getValue("dbVersion");
        //System.out.println("here" + dbAccessionId + " | " + objectVersion + " | " + intObjectId + " | " + dbSource + " | " + dbVersion + " | " + type);
        String dbCoordSys       = atts.getValue("dbCoordSys");
        //System.out.println("there" + dbCoordSys);
        
        HashMap object = new HashMap() ;
        object.put("dbAccessionId" ,dbAccessionId);
        object.put("objectVersion" ,objectVersion);
        object.put("intObjectId"   ,intObjectId);
        ArrayList details = new ArrayList();
        object.put("details",details);
        
        object.put("dbSource"      ,dbSource) ;
        //System.out.println("daga");
        if ( dbCoordSys != null ) {
            object.put("dbCoordSys"    ,dbCoordSys);
        } 
        //System.out.println("dong");
        object.put("dbVersion"     ,dbVersion) ;
        
        if ( type != null ){
            object.put("type",type); 
        } 
        
        
        current_object = object ;
        //System.out.println("done");
        
    }
    
    
    public void startDocument() {
        //System.out.println("start document");
    }
    
    public void endDocument ()	{
        
    }
    
    public void characters (char ch[], int start, int length){
        txt = "";
        for (int i = start; i < start + length; i++) {
            txt += ch[i] ;
        }
        if ( current_position == "cigar"){
            cigar += txt ;
            
        }
        if (current_position == "sequence"){
            //System.out.println(txt);
            current_object.put("sequence",txt);
        }
        
    }
    
}
