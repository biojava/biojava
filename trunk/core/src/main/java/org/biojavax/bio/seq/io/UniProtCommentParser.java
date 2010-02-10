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

package org.biojavax.bio.seq.io;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.seq.io.ParseException;
import org.biojavax.Comment;


/**
 *
 * @author Richard Holland
 * @since 1.5
 */
public class UniProtCommentParser {
    
    /**
     * Creates a new instance of UniProtCommentParser.
     */
    public UniProtCommentParser() {
        this.interactions = new ArrayList();
        this.isoforms = new ArrayList();
        this.events = new ArrayList();
        this.KMs = new ArrayList();
        this.VMaxes = new ArrayList();
        this.seqCautions = new ArrayList();
    }
    
    // the prefix for comments
    private static final String PREFIX = "-!-";
    
    /**
     * A name for a comment type.
     */
    public static final String BIOPHYSICOCHEMICAL_PROPERTIES = "BIOPHYSICOCHEMICAL PROPERTIES";
    
    /**
     * A name for a comment type.
     */
    public static final String DATABASE = "DATABASE";
    
    /**
     * A name for a comment type.
     */
    public static final String MASS_SPECTROMETRY = "MASS SPECTROMETRY";
    
    /**
     * A name for a comment type.
     */
    public static final String ALTERNATIVE_PRODUCTS = "ALTERNATIVE PRODUCTS";
    
    /**
     * A name for a comment type.
     */
    public static final String INTERACTION = "INTERACTION";
    
    /**
     * A name for a comment type.
     */
    public static final String PTM = "PTM";
    
    /**
     * A name for a comment type.
     */
    public static final String SEQUENCE_CAUTION = "SEQUENCE CAUTION";
    
    /**
     * Parses the comment string from the given comment and populates
     * the internal fields appropriately.  If the comment is not a
     * UniProt comment (does not start with -!-) then an exception is
     * thrown.
     * @param c the comment to parse.
     * @throws ParseException if the comment was not parseable.
     */
    public void parseComment(Comment c) throws ParseException {
        this.parseComment(c.getComment());
    }
    
    /**
     * Parses the comment string from the given comment and populates
     * the internal fields appropriately. If the comment is not a
     * UniProt comment (does not start with -!-) then an exception is
     * thrown.
     * @param c the comment to parse.
     * @throws ParseException if the comment was not parseable.
     */
    public void parseComment(String c) throws ParseException {
        if (!isParseable(c)) throw new ParseException("Comment is not a UniProt structured comment. Comment was "+c);
        
        String comment = new String(c); //keep the original just in case...
        // do the parsing here.
        try{
            c = c.replaceAll("\\s+", " ").trim(); // replace all multi-spaces and newlines with single spaces
            // our comment is now one long string, -!- TYPE: [prefix: key=value; | key=value; | text]
            c = c.substring(PREFIX.length()+1); // chomp "-!- "
            String type = c.substring(0,c.indexOf(':')); // find type
            this.setCommentType(type); // remember type
            c = c.substring(c.indexOf(':')+1); // chomp type and colon
            if (c.endsWith(".")) c=c.substring(0,c.length()-1); // chomp trailing full stop
            
            // what we have left is the [prefix: key=value; | key=value; | text.] section
            if (this.getCommentType().equalsIgnoreCase(BIOPHYSICOCHEMICAL_PROPERTIES)) {
            /*
CC   -!- BIOPHYSICOCHEMICAL PROPERTIES:
CC       Absorption:
CC         Abs(max)=xx nm;
CC         Note=free_text;
CC       Kinetic parameters:
CC         KM=xx unit for substrate [(free_text)];
CC         Vmax=xx unit enzyme [free_text];
CC         Note=free_text;
CC       pH dependence:
CC         free_text;
CC       Redox potential:
CC         free_text;
CC       Temperature dependence:
CC         free_text;
             */
                do {
                    String[] parts = c.split(";");
                    if (parts.length==1) {
                        // we are one of the last three options on the list
                        int firstColon = parts[0].indexOf(':');
                        String key = parts[0].substring(0,firstColon).trim();
                        String value = parts[0].substring(firstColon+1).trim();
                        if (key.equalsIgnoreCase("pH dependence")) this.setPHDependence(value);
                        else if (key.equalsIgnoreCase("Redox potential")) this.setRedoxPotential(value);
                        else if (key.equalsIgnoreCase("Temperature dependence")) this.setTemperatureDependence(value);
                        // skip to next chunk
                        c = c.substring(c.indexOf(";")+1);
                    } else {
                        // we are one of the first two options on the list
                        int skippos = -1;
                        String key = parts[0].split(":")[0].trim();
                        if (key.equalsIgnoreCase("Absorption")) {
                            String[] subparts = parts[0].split(":")[1].split("=");
                            this.setAbsorptionMax(subparts[1].trim());
                            subparts = parts[1].split("=");
                            this.setAbsorptionNote(subparts[1].trim());
                            skippos = 2;
                        } else if (key.equalsIgnoreCase("Kinetic parameters")) {
                            int partCount = 0;
                            String[] subparts = parts[partCount].split(":")[1].split("=");
                            key = subparts[0].trim();
                            String value = subparts[1].trim();
                            while (!key.equalsIgnoreCase("Note")) {
                                if (key.equalsIgnoreCase("KM")) this.getKMs().add(value);
                                else if (key.equalsIgnoreCase("VMax")) this.getVMaxes().add(value);
                                subparts = parts[++partCount].split("=");
                                key = subparts[0].trim();
                                value = subparts[1].trim();
                            }
                            this.setKineticsNote(value);
                        }
                        // skip to next chunk
                        int chunkpos = c.indexOf(parts[skippos]);
                        c = c.substring(chunkpos);
                    }
                    c = c.trim();
                } while (c.length()>0);
            } else if (this.getCommentType().equalsIgnoreCase(DATABASE)) {
            /*
CC   -!- DATABASE: NAME=Text[; NOTE=Text][; WWW="Address"][; FTP="Address"].
             */
                c = c.substring(0,c.length()-1); // chomp trailing dot
                String[] parts = c.split(";");
                for (int i = 0; i < parts.length; i++) {
                    String[] subparts = parts[i].split("=");
                    String key = subparts[0].trim();
                    String value = subparts[1].trim();
                    if (key.equalsIgnoreCase("NAME")) this.setDatabaseName(value);
                    else if (key.equalsIgnoreCase("NOTE")) this.setNote(value);
                    else if (key.equalsIgnoreCase("WWW") || key.equalsIgnoreCase("FTP")) this.setUri(value);
                }
            } else if (this.getCommentType().equalsIgnoreCase(MASS_SPECTROMETRY)) {
            /*
CC   -!- MASS SPECTROMETRY: MW=XXX[; MW_ERR=XX]; METHOD=XX; RANGE=XX-XX[ (Name)]; NOTE={Free text (Ref.n)|Ref.n}.
             */
                c = c.substring(0,c.length()-1); // chomp trailing dot
                String[] parts = c.split(";");
                for (int i = 0; i < parts.length; i++) {
                    String[] subparts = parts[i].split("=");
                    String key = subparts[0].trim();
                    String value = subparts[1].trim();
                    if (key.equalsIgnoreCase("MW")) this.setMolecularWeight(Integer.parseInt(value));
                    else if (key.equalsIgnoreCase("MW_ERR")) this.setMolWeightError(new Integer(value));
                    else if (key.equalsIgnoreCase("METHOD")) this.setMolWeightMethod(value);
                    else if (key.equalsIgnoreCase("RANGE")) {
                        if (value.indexOf(' ')>-1) value = value.substring(0, value.indexOf(' ')); // drop name
                        String[] locs = value.split("-");
                        this.setMolWeightRangeStart(Integer.parseInt(locs[0]));
                        this.setMolWeightRangeEnd(Integer.parseInt(locs[1]));
                    } else if (key.equalsIgnoreCase("NOTE")) this.setNote(value);
                }
            } else if (this.getCommentType().equalsIgnoreCase(INTERACTION)) {
            /*
CC   -!- INTERACTION:
CC       {{SP_Ac:identifier[ (xeno)]}|Self}; NbExp=n; IntAct=IntAct_Protein_Ac, IntAct_Protein_Ac;
             */
                String[] parts = c.split(";");
                Interaction interact = null;
                for (int i = 0; i < parts.length; i++) {
                    String[] subparts = parts[i].split("=");
                    String key = subparts[0].trim();
                    String value = null;
                    if (key.equalsIgnoreCase("Self")) {
                        // start new self-self interaction
                        interact = new Interaction();
                        interact.setID("Self");
                        interact.setOrganismsDiffer(false);
                        this.getInteractions().add(interact);
                    } else if (subparts.length==1) {
                        // start new protein-protein interaction
                        subparts = key.split(":");
                        boolean differ = false;
                        if (subparts[1].indexOf("(xeno)")>-1) {
                            differ = true;
                            subparts[1] = subparts[1].substring(0,subparts[1].indexOf("(xeno)"));
                        }
                        interact = new Interaction();
                        interact.setID(subparts[0].trim());
                        interact.setLabel(subparts[1].trim());
                        interact.setOrganismsDiffer(differ);
                        this.getInteractions().add(interact);
                    } else {
                        value = subparts[1].trim();
                        // continue existing interaction
                        if (key.equalsIgnoreCase("NbExp")) interact.setNumberExperiments(Integer.parseInt(value));
                        else if (key.equalsIgnoreCase("IntAct")) {
                            subparts = value.split(",");
                            interact.setFirstIntActID(subparts[0].trim());
                            interact.setSecondIntActID(subparts[1].trim());
                        }
                    }
                }
            } else if (this.getCommentType().equalsIgnoreCase(ALTERNATIVE_PRODUCTS)) {
            /*
CC   -!- ALTERNATIVE PRODUCTS:
CC       Event=Alternative promoter;
CC         Comment=Free text;
CC       Event=Alternative splicing; Named isoforms=n;
CC         Comment=Optional free text;
CC       Name=Isoform_1; Synonyms=Synonym_1[, Synonym_n];
CC         IsoId=Isoform_identifier_1[, Isoform_identifier_n]; Sequence=Displayed;
CC         Note=Free text;
CC       Name=Isoform_n; Synonyms=Synonym_1[, Synonym_n];
CC         IsoId=Isoform_identifier_1[, Isoform_identifier_n]; Sequence=VSP_identifier_1 [, VSP_identifier_n];
CC         Note=Free text;
CC       Event=Alternative initiation;
CC         Comment=Free text;
             */
                Event event = null;
                Isoform isoform = null;
                String[] parts = c.split(";");
                for (int i = 0; i < parts.length; i++) {
                    String[] subparts = parts[i].split("=");
                    if (subparts.length<2)
                        continue;
                    String key = subparts[0].trim();
                    String value = subparts[1].trim();
                    if (key.equalsIgnoreCase("Event")) {
                        // new event
                        event = new Event();
                        this.getEvents().add(event);
                        event.setType(value);
                    } else if (key.equalsIgnoreCase("Name")) {
                        // new isoform
                        isoform = new Isoform();
                        this.getIsoforms().add(isoform);
                        isoform.getNames().add(value);
                    } else if (key.equalsIgnoreCase("Synonyms")) {
                        subparts = value.split(",");
                        for (int j = 0; j < subparts.length; j++) isoform.getNames().add(subparts[j].trim());
                    } else if (key.equalsIgnoreCase("IsoId")) {
                        subparts = value.split(",");
                        for (int j = 0; j < subparts.length; j++) isoform.getIsoIDs().add(subparts[j].trim());
                    } else if (key.equalsIgnoreCase("Sequence")) {
                        if (value.equalsIgnoreCase("Displayed")) isoform.setSequenceType("Displayed");
                        else if (value.equalsIgnoreCase("Not described")) isoform.setSequenceType("Not described");
                        else if (value.equalsIgnoreCase("External")) isoform.setSequenceType("External");
                        else {
                            isoform.setSequenceType("Described");
                            isoform.setSequenceRef(value);
                        }
                    } else if (key.equalsIgnoreCase("Note")) {
                        isoform.setNote(value);
                    } else if (key.equalsIgnoreCase("Named isoforms")) {
                        event.setNamedIsoforms(Integer.parseInt(value));
                    } else if (key.equalsIgnoreCase("Comment")) {
                        event.setComment(value);
                    }
                }
            } else if (this.getCommentType().equalsIgnoreCase(SEQUENCE_CAUTION)) {
            /*
CC   -!- SEQUENCE_CAUTION: Sequence=Sequence; Type=Type;[ Positions=Positions;][ Note=Note;]
             */
                SeqCaution seqc = null;
                c = c.substring(0,c.length()-1); // chomp trailing dot
                String[] parts = c.split(";");
                for (int i = 0; i < parts.length; i++) {
                    String[] subparts = parts[i].split("=");
                    String key = subparts[0].trim();
                    String value = subparts[1].trim();
                    if (key.equalsIgnoreCase("SEQUENCE")) {
                        seqc = new SeqCaution();
                        this.getSeqCautions().add(seqc);
                        seqc.setSequence(value);
                    } else if (key.equalsIgnoreCase("TYPE")) seqc.setType(value);
                    else if (key.equalsIgnoreCase("POSITIONS")) seqc.setPositions(value);
                    else if (key.equalsIgnoreCase("NOTE")) seqc.setNote(value);
                }
            } else {
                // all others are just free text.
                this.setText(c);
            }
        }catch(RuntimeException ex){
            throw new ParseException(ex, "Cannot parse the comment: "+comment);
        }
        // all done
    }
    
    /**
     * Returns true if the comment may be parseable (starts with -!-).
     * @param c the comment to check.
     * @return true if it starts with -!-, false otherwise.
     */
    public static boolean isParseable(Comment c) {
        return isParseable(c.getComment());
    }
    
    /**
     * Returns true if the comment may be parseable (starts with -!-).
     * @param c the comment to check.
     * @return true if it starts with -!-, false otherwise.
     */
    public static boolean isParseable(String c) {
        return c.trim().startsWith(PREFIX);
    }
    
    /**
     * Generates a comment string based on the current values of the
     * internal fields.
     * @return the comment string representing the current settings.
     * @throws ParseException if the current settings do not allow the
     * creation of a correct comment string.
     */
    public String generate() throws ParseException {
        StringBuffer sb = new StringBuffer();
        sb.append(PREFIX);
        sb.append(" ");
        sb.append(this.getCommentType());
        sb.append(": ");
        
        // output the specifics
        if (this.getCommentType().equalsIgnoreCase(BIOPHYSICOCHEMICAL_PROPERTIES)) {
            /*
CC   -!- BIOPHYSICOCHEMICAL PROPERTIES:
CC       Absorption:
CC         Abs(max)=xx nm;
CC         Note=free_text;
CC       Kinetic parameters:
CC         KM=xx unit for substrate [(free_text)];
CC         Vmax=xx unit enzyme [free_text];
CC         Note=free_text;
CC       pH dependence:
CC         free_text;
CC       Redox potential:
CC         free_text;
CC       Temperature dependence:
CC         free_text;
             */
            if (this.getAbsorptionNote()!=null) {
                // we have an absorption line!
                sb.append("\nAbsorption:\n  Abs(max)=");
                sb.append(this.getAbsorptionMax());
                sb.append(";\n  Note=");
                sb.append(this.getAbsorptionNote());
                sb.append(";");
            }
            if (this.getKineticsNote()!=null) {
                // we have a kinetics note!
                sb.append("\nKinetic parameters:\n");
                for (Iterator j = this.getKMs().iterator(); j.hasNext(); ) {
                    sb.append("  KM=");
                    sb.append((String)j.next());
                    sb.append(";\n");
                }
                for (Iterator j = this.getVMaxes().iterator(); j.hasNext(); ) {
                    sb.append("  VMax=");
                    sb.append((String)j.next());
                    sb.append(";\n");
                }
                sb.append("  Note=");
                sb.append(this.getKineticsNote());
                sb.append(";");
            }
            if (this.getPHDependence()!=null) {
                sb.append("\npH dependence:\n  ");
                sb.append(this.getPHDependence());
                sb.append(";");
            }
            if (this.getRedoxPotential()!=null) {
                sb.append("\nRedox potential:\n  ");
                sb.append(this.getRedoxPotential());
                sb.append(";");
            }
            if (this.getTemperatureDependence()!=null) {
                sb.append("\nTemperature dependence:\n  ");
                sb.append(this.getTemperatureDependence());
                sb.append(";");
            }
            
        } else if (this.getCommentType().equalsIgnoreCase(DATABASE)) {
            if (this.getDatabaseName()==null) throw new ParseException("Database name is missing");
            /*
CC   -!- DATABASE: NAME=Text[; NOTE=Text][; WWW="Address"][; FTP="Address"].
             */
            sb.append("NAME=");
            sb.append(this.getDatabaseName());
            if (this.getNote()!=null) {
                sb.append("; NOTE=");
                sb.append(this.getNote());
            }
            if (this.getUri()!=null) {
                sb.append("; ");
                if (this.getUri().startsWith("ftp")) sb.append(" FTP=");
                else sb.append(" WWW=");
                sb.append(this.getUri());
            }
            sb.append(".");
            
        } else if (this.getCommentType().equalsIgnoreCase(MASS_SPECTROMETRY)) {
            /*
CC   -!- MASS SPECTROMETRY: MW=XXX[; MW_ERR=XX]; METHOD=XX; RANGE=XX-XX[ (Name)]; NOTE={Free text (Ref.n)|Ref.n}.
             */
            sb.append("MW=");
            sb.append(""+this.getMolecularWeight());
            if (this.getMolWeightError()!=null) {
                sb.append("; MW_ERR=");
                sb.append(""+this.getMolWeightError());
            }
            sb.append("; METHOD=");
            sb.append(this.getMolWeightMethod());
            sb.append("; RANGE=");
            sb.append(""+this.getMolWeightRangeStart());
            sb.append("-");
            sb.append(""+this.getMolWeightRangeEnd());
            sb.append("; NOTE=");
            sb.append(this.getNote());
            sb.append(".");
            
        } else if (this.getCommentType().equalsIgnoreCase(INTERACTION)) {
            /*
CC   -!- INTERACTION:
CC       {{SP_Ac:identifier[ (xeno)]}|Self}; NbExp=n; IntAct=IntAct_Protein_Ac, IntAct_Protein_Ac;
             */
            for (Iterator i = this.getInteractions().iterator(); i.hasNext(); ) {
                Interaction interact = (Interaction)i.next();
                sb.append("\n"); // each interaction starts on a new line
                if (interact.getID().equals("Self")) {
                    sb.append("Self; ");
                } else {
                    sb.append(interact.getID());
                    sb.append(":");
                    sb.append(interact.getLabel());
                    if (interact.isOrganismsDiffer()) sb.append(" (xeno)");
                    sb.append("; ");
                }
                sb.append("NbExp=");
                sb.append(""+interact.getNumberExperiments());
                sb.append("; ");
                sb.append("IntAct=");
                sb.append(interact.getFirstIntActID());
                sb.append(", ");
                sb.append(interact.getSecondIntActID());
                sb.append(";");
            }
            
        } else if (this.getCommentType().equalsIgnoreCase(ALTERNATIVE_PRODUCTS)) {
            /*
CC   -!- ALTERNATIVE PRODUCTS:
CC       Event=Alternative promoter;
CC         Comment=Free text;
CC       Event=Alternative splicing; Named isoforms=n;
CC         Comment=Optional free text;
CC       Name=Isoform_1; Synonyms=Synonym_1[, Synonym_n];
CC         IsoId=Isoform_identifier_1[, Isoform_identifier_n]; Sequence=Displayed;
CC         Note=Free text;
CC       Name=Isoform_n; Synonyms=Synonym_1[, Synonym_n];
CC         IsoId=Isoform_identifier_1[, Isoform_identifier_n]; Sequence=VSP_identifier_1 [, VSP_identifier_n];
CC         Note=Free text;
CC       Event=Alternative initiation;
CC         Comment=Free text;
             */
            for (Iterator i = this.getEvents().iterator(); i.hasNext(); ) {
                Event event = (Event)i.next();
                sb.append("\n"); // each one starts on a new line
                sb.append("Event=");
                sb.append(event.getType());
                if (event.getType().equals("Alternative splicing")) {
                    sb.append("; Named isoforms=");
                    sb.append(""+event.getNamedIsoforms());
                }
                sb.append(";\n  Comment="); // comment is indented two on next line
                sb.append(event.getComment());
                sb.append(";");
                if (event.getType().equals("Alternative splicing")) {
                    for (Iterator j = this.getIsoforms().iterator(); j.hasNext(); ) {
                        Isoform isoform = (Isoform)j.next();
                        sb.append("\nName="); // each isoform on a new line
                        sb.append(isoform.getNames().get(0));
                        sb.append("; Synonyms=");
                        for (int k =1 ; k < isoform.getNames().size(); k++) {
                            sb.append(isoform.getNames().get(k));
                            if (k<isoform.getNames().size()-1) sb.append(", ");
                        }
                        sb.append(";\n  IsoId="); // isoid on new line indented by two
                        sb.append(isoform.getIsoIDs().get(0));
                        for (int k =1 ; k < isoform.getIsoIDs().size(); k++) {
                            sb.append(isoform.getIsoIDs().get(k));
                            if (k<isoform.getIsoIDs().size()-1) sb.append(", ");
                        }
                        sb.append("; Sequence=");
                        if (isoform.getSequenceRef()!=null) sb.append(isoform.getSequenceRef());
                        else sb.append(isoform.getSequenceType());
                        sb.append(";\n  Note="); // note is indented by two as well
                        sb.append(isoform.getNote());
                        sb.append(";");
                    }
                }
            }
        } else if (this.getCommentType().equalsIgnoreCase(SEQUENCE_CAUTION)) {
                /*
CC   -!- SEQUENCE_CAUTION: Sequence=Sequence; Type=Type;[ Positions=Positions;][ Note=Note;]
                 */
            for (Iterator i = this.getSeqCautions().iterator(); i.hasNext(); ) {
                SeqCaution seqc = (SeqCaution)i.next();
                sb.append("\n"); // each one starts on a new line
                sb.append("Sequence=");
                sb.append(seqc.getSequence());
                sb.append("; Type=");
                sb.append(seqc.getType());
                if (seqc.getPositions()!=null) {
                    sb.append("; Positions=");
                    sb.append(seqc.getPositions());
                }
                if (this.getNote()!=null) {
                    sb.append("; Note=");
                    sb.append(seqc.getNote());
                }
                sb.append(";");
            }
        } else {
            // just append free text for all others.
            sb.append(this.getText());
            if (!this.getText().endsWith(".")) sb.append(".");
        }
        
        // return it
        return sb.toString();
    }
    
    /**
     * Holds value of property commentType.
     */
    private String commentType;
    
    /**
     * Getter for property commentType.
     * @return Value of property commentType.
     */
    public String getCommentType() {
        
        return this.commentType;
    }
    
    /**
     * Setter for property commentType.
     * @param commentType New value of property commentType.
     */
    public void setCommentType(String commentType) {
        
        this.commentType = commentType;
    }
    
    /**
     * Holds value of property text.
     */
    private String text;
    
    /**
     * Getter for property text.
     * @return Value of property text.
     */
    public String getText() {
        
        return this.text;
    }
    
    /**
     * Setter for property text.
     * @param text New value of property text.
     */
    public void setText(String text) {
        
        this.text = text;
    }
    
    /**
     * Holds value of property databaseName.
     */
    private String databaseName;
    
    /**
     * Getter for property databaseName.
     * @return Value of property databaseName.
     */
    public String getDatabaseName() {
        
        return this.databaseName;
    }
    
    /**
     * Setter for property databaseName.
     * @param databaseName New value of property databaseName.
     */
    public void setDatabaseName(String databaseName) {
        
        this.databaseName = databaseName;
    }
    
    /**
     * Holds value of property note.
     */
    private String note;
    
    /**
     * Getter for property note.
     * @return Value of property note.
     */
    public String getNote() {
        
        return this.note;
    }
    
    /**
     * Setter for property note.
     * @param note New value of property note.
     */
    public void setNote(String note) {
        
        this.note = note;
    }
    
    /**
     * Holds value of property uri.
     */
    private String uri;
    
    /**
     * Getter for property uri.
     * @return Value of property uri.
     */
    public String getUri() {
        
        return this.uri;
    }
    
    /**
     * Setter for property uri.
     * @param uri New value of property uri.
     */
    public void setUri(String uri) {
        
        this.uri = uri;
    }
    
    /**
     * Holds value of property molecularWeight.
     */
    private int molecularWeight;
    
    /**
     * Getter for property molecularWeight.
     * @return Value of property molecularWeight.
     */
    public int getMolecularWeight() {
        
        return this.molecularWeight;
    }
    
    /**
     * Setter for property molecularWeight.
     * @param molecularWeight New value of property molecularWeight.
     */
    public void setMolecularWeight(int molecularWeight) {
        
        this.molecularWeight = molecularWeight;
    }
    
    /**
     * Holds value of property molWeightError.
     */
    private Integer molWeightError;
    
    /**
     * Getter for property molWeightError.
     * @return Value of property molWeightError.
     */
    public Integer getMolWeightError() {
        
        return this.molWeightError;
    }
    
    /**
     * Setter for property molWeightError.
     * @param molWeightError New value of property molWeightError.
     */
    public void setMolWeightError(Integer molWeightError) {
        
        this.molWeightError = molWeightError;
    }
    
    /**
     * Holds value of property molWeightRangeStart.
     */
    private int molWeightRangeStart;
    
    /**
     * Getter for property molWeightRangeStart.
     * @return Value of property molWeightRangeStart.
     */
    public int getMolWeightRangeStart() {
        
        return this.molWeightRangeStart;
    }
    
    /**
     * Setter for property molWeightRangeStart.
     * @param molWeightRangeStart New value of property molWeightRangeStart.
     */
    public void setMolWeightRangeStart(int molWeightRangeStart) {
        
        this.molWeightRangeStart = molWeightRangeStart;
    }
    
    /**
     * Holds value of property molWeightRangeEnd.
     */
    private int molWeightRangeEnd;
    
    /**
     * Getter for property molWeightRangeEnd.
     * @return Value of property molWeightRangeEnd.
     */
    public int getMolWeightRangeEnd() {
        
        return this.molWeightRangeEnd;
    }
    
    /**
     * Setter for property molWeightRangeEnd.
     * @param molWeightRangeEnd New value of property molWeightRangeEnd.
     */
    public void setMolWeightRangeEnd(int molWeightRangeEnd) {
        
        this.molWeightRangeEnd = molWeightRangeEnd;
    }
    
    /**
     * Holds value of property molWeightMethod.
     */
    private String molWeightMethod;
    
    /**
     * Getter for property molWeightMethod.
     * @return Value of property molWeightMethod.
     */
    public String getMolWeightMethod() {
        
        return this.molWeightMethod;
    }
    
    /**
     * Setter for property molWeightMethod.
     * @param molWeightMethod New value of property molWeightMethod.
     */
    public void setMolWeightMethod(String molWeightMethod) {
        
        this.molWeightMethod = molWeightMethod;
    }
    
    /**
     * Holds value of property interactions.
     */
    private List interactions;
    
    /**
     * Getter for property interactions.
     * @return Value of property interactions.
     */
    public List getInteractions() {
        
        return this.interactions;
    }
    
    /**
     * Setter for property interactions.
     * @param interactions New value of property interactions.
     */
    public void setInteractions(List interactions) {
        
        this.interactions = interactions;
    }
    
    /**
     * Holds value of property seqCautions.
     */
    private List seqCautions;
    
    /**
     * Getter for property seqCautions.
     * @return Value of property seqCautions.
     */
    public List getSeqCautions() {
        
        return this.seqCautions;
    }
    
    /**
     * Setter for property seqCautions.
     * @param seqCautions New value of property seqCautions.
     */
    public void setSeqCautions(List seqCautions) {
        
        this.seqCautions = seqCautions;
    }
    
    /**
     * Holds value of property events.
     */
    private List events;
    
    /**
     * Getter for property events.
     * @return Value of property events.
     */
    public List getEvents() {
        
        return this.events;
    }
    
    /**
     * Setter for property events.
     * @param events New value of property events.
     */
    public void setEvents(List events) {
        
        this.events = events;
    }
    
    /**
     * Holds value of property isoforms.
     */
    private List isoforms;
    
    /**
     * Getter for property isoforms.
     * @return Value of property isoforms.
     */
    public List getIsoforms() {
        
        return this.isoforms;
    }
    
    /**
     * Setter for property isoforms.
     * @param isoforms New value of property isoforms.
     */
    public void setIsoforms(List isoforms) {
        
        this.isoforms = isoforms;
    }
    
    /**
     * Holds value of property absorptionMax.
     */
    private String absorptionMax;
    
    /**
     * Getter for property absorptionMax.
     * @return Value of property absorptionMax.
     */
    public String getAbsorptionMax() {
        
        return this.absorptionMax;
    }
    
    /**
     * Setter for property absorptionMax.
     * @param absorptionMax New value of property absorptionMax.
     */
    public void setAbsorptionMax(String absorptionMax) {
        
        this.absorptionMax = absorptionMax;
    }
    
    /**
     * Holds value of property absorptionNote.
     */
    private String absorptionNote;
    
    /**
     * Getter for property absorptionNote.
     * @return Value of property absorptionNote.
     */
    public String getAbsorptionNote() {
        
        return this.absorptionNote;
    }
    
    /**
     * Setter for property absorptionNote.
     * @param absorptionNote New value of property absorptionNote.
     */
    public void setAbsorptionNote(String absorptionNote) {
        
        this.absorptionNote = absorptionNote;
    }
    
    /**
     * Holds value of property KMs.
     */
    private List KMs;
    
    /**
     * Getter for property KMs.
     * @return Value of property KMs.
     */
    public List getKMs() {
        
        return this.KMs;
    }
    
    /**
     * Setter for property KMs.
     * @param KMs New value of property KMs.
     */
    public void setKMs(List KMs) {
        
        this.KMs = KMs;
    }
    
    /**
     * Holds value of property VMaxes.
     */
    private List VMaxes;
    
    /**
     * Getter for property VMaxes.
     * @return Value of property VMaxes.
     */
    public List getVMaxes() {
        
        return this.VMaxes;
    }
    
    /**
     * Setter for property VMaxes.
     * @param VMaxes New value of property VMaxes.
     */
    public void setVMaxes(List VMaxes) {
        
        this.VMaxes = VMaxes;
    }
    
    /**
     * Holds value of property kineticsNote.
     */
    private String kineticsNote;
    
    /**
     * Getter for property kineticsNote.
     * @return Value of property kineticsNote.
     */
    public String getKineticsNote() {
        
        return this.kineticsNote;
    }
    
    /**
     * Setter for property kineticsNote.
     * @param kineticsNote New value of property kineticsNote.
     */
    public void setKineticsNote(String kineticsNote) {
        
        this.kineticsNote = kineticsNote;
    }
    
    /**
     * Holds value of property PHDependence.
     */
    private String PHDependence;
    
    /**
     * Getter for property PHDependence.
     * @return Value of property PHDependence.
     */
    public String getPHDependence() {
        
        return this.PHDependence;
    }
    
    /**
     * Setter for property PHDependence.
     * @param PHDependence New value of property PHDependence.
     */
    public void setPHDependence(String PHDependence) {
        
        this.PHDependence = PHDependence;
    }
    
    /**
     * Holds value of property redoxPotential.
     */
    private String redoxPotential;
    
    /**
     * Getter for property redoxPotential.
     * @return Value of property redoxPotential.
     */
    public String getRedoxPotential() {
        
        return this.redoxPotential;
    }
    
    /**
     * Setter for property redoxPotential.
     * @param redoxPotential New value of property redoxPotential.
     */
    public void setRedoxPotential(String redoxPotential) {
        
        this.redoxPotential = redoxPotential;
    }
    
    /**
     * Holds value of property temperatureDependence.
     */
    private String temperatureDependence;
    
    /**
     * Getter for property temperatureDependence.
     * @return Value of property temperatureDependence.
     */
    public String getTemperatureDependence() {
        
        return this.temperatureDependence;
    }
    
    /**
     * Setter for property temperatureDependence.
     * @param temperatureDependence New value of property temperatureDependence.
     */
    public void setTemperatureDependence(String temperatureDependence) {
        
        this.temperatureDependence = temperatureDependence;
    }
    
    /**
     * A class to describe protein-protein interactions.
     */
    public static class Interaction {
        /**
         * Holds value of property ID.
         */
        private String ID;
        
        /**
         * Getter for property ID.
         * @return Value of property ID.
         */
        public String getID() {
            
            return this.ID;
        }
        
        /**
         * Setter for property ID.
         * @param ID New value of property ID.
         */
        public void setID(String ID) {
            
            this.ID = ID;
        }
        
        /**
         * Holds value of property label.
         */
        private String label;
        
        /**
         * Getter for property label.
         * @return Value of property label.
         */
        public String getLabel() {
            
            return this.label;
        }
        
        /**
         * Setter for property label.
         * @param label New value of property label.
         */
        public void setLabel(String label) {
            
            this.label = label;
        }
        
        /**
         * Holds value of property organismsDiffer.
         */
        private boolean organismsDiffer;
        
        /**
         * Getter for property organismsDiffer.
         * @return Value of property organismsDiffer.
         */
        public boolean isOrganismsDiffer() {
            
            return this.organismsDiffer;
        }
        
        /**
         * Setter for property organismsDiffer.
         * @param organismsDiffer New value of property organismsDiffer.
         */
        public void setOrganismsDiffer(boolean organismsDiffer) {
            
            this.organismsDiffer = organismsDiffer;
        }
        
        /**
         * Holds value of property firstIntActID.
         */
        private String firstIntActID;
        
        /**
         * Getter for property firstIntActID.
         * @return Value of property firstIntActID.
         */
        public String getFirstIntActID() {
            
            return this.firstIntActID;
        }
        
        /**
         * Setter for property firstIntActID.
         * @param firstIntActID New value of property firstIntActID.
         */
        public void setFirstIntActID(String firstIntActID) {
            
            this.firstIntActID = firstIntActID;
        }
        
        /**
         * Holds value of property secondIntActID.
         */
        private String secondIntActID;
        
        /**
         * Getter for property secondIntActID.
         * @return Value of property secondIntActID.
         */
        public String getSecondIntActID() {
            
            return this.secondIntActID;
        }
        
        /**
         * Setter for property secondIntActID.
         * @param secondIntActID New value of property secondIntActID.
         */
        public void setSecondIntActID(String secondIntActID) {
            
            this.secondIntActID = secondIntActID;
        }
        
        /**
         * Holds value of property numberExperiments.
         */
        private int numberExperiments;
        
        /**
         * Getter for property numberExperiments.
         * @return Value of property numberExperiments.
         */
        public int getNumberExperiments() {
            
            return this.numberExperiments;
        }
        
        /**
         * Setter for property numberExperiments.
         * @param numberExperiments New value of property numberExperiments.
         */
        public void setNumberExperiments(int numberExperiments) {
            
            this.numberExperiments = numberExperiments;
        }
    }
    
    /**
     * A class to describe events for alternative product comments.
     */
    public static class Event {
        /**
         * Holds value of property type.
         */
        private String type;
        
        /**
         * Getter for property type.
         * @return Value of property type.
         */
        public String getType() {
            
            return this.type;
        }
        
        /**
         * Setter for property type.
         * @param type New value of property type.
         */
        public void setType(String type) {
            
            this.type = type;
        }
        
        /**
         * Holds value of property comment.
         */
        private String comment;
        
        /**
         * Getter for property comment.
         * @return Value of property comment.
         */
        public String getComment() {
            
            return this.comment;
        }
        
        /**
         * Setter for property comment.
         * @param comment New value of property comment.
         */
        public void setComment(String comment) {
            
            this.comment = comment;
        }
        
        /**
         * Holds value of property namedIsoforms.
         */
        private int namedIsoforms;
        
        /**
         * Getter for property namedIsoforms.
         * @return Value of property namedIsoforms.
         */
        public int getNamedIsoforms() {
            
            return this.namedIsoforms;
        }
        
        /**
         * Setter for property namedIsoforms.
         * @param namedIsoforms New value of property namedIsoforms.
         */
        public void setNamedIsoforms(int namedIsoforms) {
            
            this.namedIsoforms = namedIsoforms;
        }
        
    }
    
    /**
     * A class to describe isoforms for alternative product comments.
     */
    public static class Isoform {
        
        /**
         * Creates a new instance.
         */
        public Isoform() {
            this.names = new ArrayList();
            this.isoIDs = new ArrayList();
        }
        
        /**
         * Holds value of property names.
         */
        private List names;
        
        /**
         * Getter for property names.
         * @return Value of property names.
         */
        public List getNames() {
            
            return this.names;
        }
        
        /**
         * Setter for property names.
         * @param names New value of property names.
         */
        public void setNames(List names) {
            
            this.names = names;
        }
        
        /**
         * Holds value of property isoIDs.
         */
        private List isoIDs;
        
        /**
         * Getter for property isoIDs.
         * @return Value of property isoIDs.
         */
        public List getIsoIDs() {
            
            return this.isoIDs;
        }
        
        /**
         * Setter for property isoIDs.
         * @param isoIDs New value of property isoIDs.
         */
        public void setIsoIDs(List isoIDs) {
            
            this.isoIDs = isoIDs;
        }
        
        /**
         * Holds value of property sequenceType.
         */
        private String sequenceType;
        
        /**
         * Getter for property sequenceType.
         * @return Value of property sequenceType.
         */
        public String getSequenceType() {
            
            return this.sequenceType;
        }
        
        /**
         * Setter for property sequenceType.
         * @param sequenceType New value of property sequenceType.
         */
        public void setSequenceType(String sequenceType) {
            
            this.sequenceType = sequenceType;
        }
        
        /**
         * Holds value of property sequenceRef.
         */
        private String sequenceRef;
        
        /**
         * Getter for property sequenceRef.
         * @return Value of property sequenceRef.
         */
        public String getSequenceRef() {
            
            return this.sequenceRef;
        }
        
        /**
         * Setter for property sequenceRef.
         * @param sequenceRef New value of property sequenceRef.
         */
        public void setSequenceRef(String sequenceRef) {
            
            this.sequenceRef = sequenceRef;
        }
        
        /**
         * Holds value of property note.
         */
        private String note;
        
        /**
         * Getter for property note.
         * @return Value of property note.
         */
        public String getNote() {
            
            return this.note;
        }
        
        /**
         * Setter for property note.
         * @param note New value of property note.
         */
        public void setNote(String note) {
            
            this.note = note;
        }
        
    }
    
    /**
     * A class to describe seq caution entries.
     */
    public static class SeqCaution {
        
        private String sequence;
        
        private String type;
        
        private String positions;
        
        private String note;
        
        /**
         * Creates a new instance.
         */
        public SeqCaution() {
        }
        
        /**
         * @return the note
         */
        public String getNote() {
            return note;
        }
        
        /**
         * @param note the note to set
         */
        public void setNote(String note) {
            this.note = note;
        }
        
        /**
         * @return the positions
         */
        public String getPositions() {
            return positions;
        }
        
        /**
         * @param positions the positions to set
         */
        public void setPositions(String positions) {
            this.positions = positions;
        }
        
        /**
         * @return the sequence
         */
        public String getSequence() {
            return sequence;
        }
        
        /**
         * @param sequence the sequence to set
         */
        public void setSequence(String sequence) {
            this.sequence = sequence;
        }
        
        /**
         * @return the type
         */
        public String getType() {
            return type;
        }
        
        /**
         * @param type the type to set
         */
        public void setType(String type) {
            this.type = type;
        }
        
    }
}
