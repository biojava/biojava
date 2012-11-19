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
 * Created on Sep 3, 2007
 *
 */
package org.biojava.bio.structure;

import java.io.Serializable;
import java.lang.reflect.Method;
import java.util.Formatter;
import java.util.Locale;


/** A class to represent database cross references. This is just a simple bean that contains the infor from one
 * DBREF line
 *
 * @author Andreas Prlic
 * @since 4:56:14 PM
 * @version %I% %G%
 */
public class DBRef implements PDBRecord, Serializable{

	/**
    *
    */
   private static final long serialVersionUID = -1050178577542224379L;
   Structure parent;
	String idCode;
    Character chainId;
    int seqbegin;
    char insertBegin;
    int seqEnd;
    char insertEnd;

    String database;
    String dbAccession;
    String dbIdCode;

    int dbSeqBegin;
    char idbnsBegin;
    int dbSeqEnd;
    char idbnsEnd;

    private Long id;

    public DBRef() {
    	insertBegin = ' ';
    	insertEnd   = ' ';
    	idbnsBegin  = ' ';
    	idbnsEnd    = ' ';
    }

    /** Get the ID used by Hibernate.
     *
     * @return the ID used by Hibernate
     * @see #setId(Long)
     */
    public Long getId() {
        return id;
    }

    /** Set the ID used by Hibernate.
     *
     * @param id the id assigned by Hibernate
     * @see #getId()
     */
    public void setId(Long id) {
        this.id = id;
    }

    /** Set the structure object that this DBRef relates to.
     *
     * @param s a structure object
     * @see #getParent()
     */
    public void setParent(Structure s){
    	parent = s;

    }

    /** Get the structure object that this DBRef relates to.
     *
     * @return s a structure object
     * @see #setParent(Structure)
     */
    public Structure getParent(){
    	return parent;
    }

    /** Convert the DBRef object to a DBREF record as it is used in PDB files
     *
     * @return a PDB - DBREF formatted line
     */
    public String toPDB(){



         StringBuffer buf = new StringBuffer();
         toPDB(buf);
         return buf.toString();

    }

    /** Append the PDB representation of this DBRef to the provided StringBuffer
     *
     * @param buf the StringBuffer to write to.
     */
    public void toPDB(StringBuffer buf){
    	 Formatter formatter = new Formatter(new StringBuilder(),Locale.UK);
//        DBREF  3ETA A  990  1295  UNP    P06213   INSR_HUMAN    1017   1322
//        DBREF  3EH2 A    2   767  UNP    P53992   SC24C_HUMAN    329   1094
//        DBREF 3EH2 A    2   767     UNP   P53992  SC24C_HUMAN   329   1094
//        DBREF  3ETA A  990  1295  UNP    P06213   INSR_HUMAN    1017   1322
        formatter.format("DBREF  %4s %1s %4d%1s %4d%1s %-6s %-8s %-12s%6d%1c%6d%1c",
                idCode, chainId,seqbegin,insertBegin,seqEnd,insertEnd,
                database,dbAccession,dbIdCode,
                dbSeqBegin,idbnsBegin,dbSeqEnd,idbnsEnd
                );

        buf.append(formatter.toString().trim());
        formatter.close();

    }
    /** String representation of a DBRef.
     * @return a String
     */
    public String toString(){
        StringBuffer buf = new StringBuffer();

        try {

            @SuppressWarnings("rawtypes")
			Class c = Class.forName("org.biojava.bio.structure.DBRef");
            Method[] methods  = c.getMethods();

            for (int i = 0; i < methods.length; i++) {
                Method m = methods[i];

                String name = m.getName();

                if ( name.substring(0,3).equals("get")) {
                    if (name.equals("getClass"))
                            continue;
                    Object o  = m.invoke(this, new Object[]{});
                    if ( o != null){
                        buf.append(name.substring(3,name.length()));
                        buf.append(": " + o + " ");
                    }
                }
            }
        } catch (Exception e){
            e.printStackTrace();
        }

        return buf.toString();
    }


    /** get the idCode for this entry
     *
     * @return the idCode
     * @see #setIdCode(String)
     */
    public String getIdCode() {
		return idCode;
	}

    /** Set the idCode for this entry.
     *
     * @param idCode the idCode for this entry
     * @see #getIdCode()
     */
	public void setIdCode(String idCode) {
		this.idCode = idCode;
	}

	/** The chain ID of the corresponding chain.
     *
     * @return chainId the ID of the corresponding chain.
     * @see #setChainId(Character)
     */
    public Character getChainId() {
        return chainId;
    }


    /** The chain ID of the corresponding chain.
     *
     * @param chainId the ID of the corresponding chain
     * @see #getChainId()
     */
    public void setChainId(Character chainId) {
        this.chainId = chainId;
    }


    /** The database of the db-ref.
     * uses the abbreviation as provided in the PDB files:
     *
     *<pre>   Database name                         database
                                     (code in columns 27 - 32)
    ----------------------------------------------------------
    GenBank                               GB
    Protein Data Bank                     PDB
    Protein Identification Resource       PIR
    SWISS-PROT                            SWS
    TREMBL                                TREMBL
    UNIPROT                               UNP
    </pre>
     * @return name of database of this DBRef
     * @see #setDatabase(String)
     */
    public String getDatabase() {
        return database;
    }

    /** Specifies the database value.
     *
     * @param database the database
     * @see #getDatabase()
     */
    public void setDatabase(String database) {
        this.database = database;
    }

    /** Sequence database accession code.
     * @return the dbAccession
     * @see #setDbAccession(String)
     * */
    public String getDbAccession() {
        return dbAccession;
    }

    /** Sequence database accession code.
     * @param dbAccession the dbAccession
     * @see #getDbAccession()
     * */
    public void setDbAccession(String dbAccession) {
        this.dbAccession = dbAccession;
    }


    /** Sequence database          identification code.
     *
     * @return the dbIdCode
     * @see #setDbIdCode(String)
     */
    public String getDbIdCode() {
        return dbIdCode;
    }

    /** Sequence database          identification code.
     *
     * @param dbIdCode identification code
     * @see #getDbIdCode()
     */
    public void setDbIdCode(String dbIdCode) {
        this.dbIdCode = dbIdCode;
    }

    /** Initial sequence number of the
    database seqment.
     * @return position
     * @see #setDbSeqBegin(int)
     */
    public int getDbSeqBegin() {
        return dbSeqBegin;
    }


    /** Initial sequence number of the
    database seqment.
     * @param dbSeqBegin a sequence position
     * @see #getDbSeqBegin()
     *
     */
    public void setDbSeqBegin(int dbSeqBegin) {
        this.dbSeqBegin = dbSeqBegin;
    }


    /** Ending sequence position  of the database segment.
     * @return dbSeqEnd
     * @see #setDbSeqEnd(int)
     */
    public int getDbSeqEnd() {
        return dbSeqEnd;
    }


    /** The begin of the sequence position in the database
     *
     * @param dbSeqEnd sequence position
     * @see #getDbSeqEnd()
     */
    public void setDbSeqEnd(int dbSeqEnd) {
        this.dbSeqEnd = dbSeqEnd;
    }

    /** Insertion code of initial residue of the segment, if PDB is the
    reference.
     * @return idbnsBegin isnertion code
     * @see #setIdbnsBegin(char)
     * */
    public char getIdbnsBegin() {
        return idbnsBegin;
    }

    /** Insertion code of initial residue of the segment, if PDB is the
    reference.
     * @param idbnsBegin insertion code
     * @see #getIdbnsBegin()
     * */
    public void setIdbnsBegin(char idbnsBegin) {
        this.idbnsBegin = idbnsBegin;
    }

    /** Insertion code of the ending
    residue of the segment, if PDB is
    the reference.
     * @return idbnsEnd insertion code
     * @see #setIdbnsEnd(char)
     */
    public char getIdbnsEnd() {
        return idbnsEnd;
    }


    /** Insertion code of the ending
    residue of the segment, if PDB is
    the reference.
     * @param idbnsEnd the insertion code
     * @see #setIdbnsEnd(char)
     */
    public void setIdbnsEnd(char idbnsEnd) {
        this.idbnsEnd = idbnsEnd;
    }

    /** Initial insertion code of the PDB sequence segment.
     *
     * @return insertBegin
     * @see #setInsertBegin(char)
     */

    public char getInsertBegin() {
        return insertBegin;
    }

    /** Initial insertion code of the PDB sequence segment.
     *
     * @param insertBegin
     * @see #getInsertBegin()
     */

    public void setInsertBegin(char insertBegin) {
        this.insertBegin = insertBegin;
    }

    /** Ending insertion code of the PDB sequence segment.
     *
     * @return insertEnd insertion code
     * @see #setInsertEnd(char)
     */
    public char getInsertEnd() {
        return insertEnd;
    }

    /** Ending insertion code of the PDB sequence segment.
     *
     * @param insertEnd insertEnd
     * @see #getInsertEnd()
     *
     */
    public void setInsertEnd(char insertEnd) {
        this.insertEnd = insertEnd;
    }

    /**   Initial sequence number of the PDB sequence segment.
     *
     * @return start seq. position
     * @see #setSeqBegin
     */
    public int getSeqBegin() {
        return seqbegin;
    }

    /**   Initial sequence number of the PDB sequence segment.
     *
     * @param seqbegin start seq. position
     * @see #getSeqBegin()
     */
    public void setSeqBegin(int seqbegin) {
        this.seqbegin = seqbegin;
    }

    /**Ending sequence number   of the PDB sequence segment.
     *
     * @return sequence end position
     * @see #getSeqEnd()
     */
    public int getSeqEnd() {
        return seqEnd;
    }

    /**Ending sequence number   of the PDB sequence segment.
     *
     * @param seqEnd sequence end position
     * @see #setSeqEnd(int)
     *
     */
    public void setSeqEnd(int seqEnd) {
        this.seqEnd = seqEnd;
    }

}
