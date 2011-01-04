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
 * Created on 28.04.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.bio.structure;

import java.io.Serializable;

import org.biojava.bio.structure.io.FileConvert;


/**
 * Implementation of an Atom of a PDB file.
 * currently the coordinates of an atom are represented by a double[3] array.
 * @author Andreas Prlic
 * @since 1.4
 * @version %I% %G%
 */
public class AtomImpl implements Atom,Serializable, PDBRecord {

    /**
    *
    */
   private static final long serialVersionUID = -2258364127420562883L;
   String name     ;
    String fullName ;
    Element element;
    double[] coords ;
    String pdbline  ;
    int pdbserial   ;

    double occupancy ;
    double tempfactor;

    Character altLoc ;
    Group parent;
    long id;

    public AtomImpl () {
        name     = null        ;
        element = Element.R;
        fullName = null        ;
        coords   = new double[3];
        pdbline  = ""          ;
        occupancy  = 0.0       ;
        tempfactor = 0.0       ;
        altLoc = new Character(' ');
        parent = null;
    }
    /** Get the Hibernate database ID.
     *
     * @return the id
     * @see #setId(long)
     */
	public long getId() {
		return id;
	}

	/** Set the Hibernate database ID.
     *
     * @param id the hibernate id
     * @see #getId()
     */
	public void setId(long id) {
		this.id = id;
	}

    /** trimmed version of atom name, e.g. "CA"
     * @see #getName
     */
    public void   setName(String s) { name = s ;}

    /**
     * Gets this object's name.
     * @return a String representing the name value
     * @see #setName
     */
    public String getName()         { return name ;}

    /** set full name of atom e.g. " CA " .
     * @see #getFullName
     */
    public void   setFullName( String s ) { fullName =s ; }

    /** get full name of atom e.g. " CA ".
     * @return a String representing the full name value
     * @see #setFullName
     */
    public String getFullName()           { return fullName; }

    /** set PDB atom number.
     * @see #getPDBserial
     */
    public void setPDBserial(int i) { pdbserial = i    ; }

    /** get PDB atom number.
     * @see #setPDBserial
     */
    public int  getPDBserial()      { return pdbserial ; }

    /** the coordinates.
     * @see #getCoords
     */
    public void     setCoords( double[] c ) { coords = c   ; }

    /** get the coordinates as a double[3] array .
     * @return an array of doubles representing the coords value
     * @see #setCoords
     */
    public double[] getCoords()            { return coords ; }

    public void setX(double x) {
        coords[0] = x ;
    }
    public void setY(double y) {
        coords[1] = y ;
    }
    public void setZ(double z) {
        coords[2] = z ;
    }

    /** Get the X coordinate.
     * @see #setX(double)
     * */
    public double getX() { return coords[0]; }

    /** Get the Y coordinate.
     * @see #setY(double)
     * */
    public double getY() { return coords[1]; }

    /** Get the Z coordinate.
     * @see #setZ(double)
     *
     * */
    public double getZ() { return coords[2]; }

    /** set alternate Location.
     * @see #getAltLoc
     */
    public void setAltLoc(Character c) {
        altLoc = c ;
    }
    /** get alternate Location.
     * @return a Character object representing the alt loc value
     * @see #setAltLoc
     */
    public Character getAltLoc() {
        return altLoc ;
    }


    /** store the whole line.
     * @see #getPDBline
     * @deprecated
     */
    public void   setPDBline(String s) { pdbline = s;}

    /** get the whole line .
     * @return a String representing the PDBline value
     * @see #setPDBline
     * @deprecated
     */
    public String getPDBline() { return pdbline ;}

    /** string representation. */
    public String toString() {
        String str = fullName +" (" + name +") " + element + " " + pdbserial + " " + coords[0] + " " + coords[1] + " " + coords[2];
        return str ;
    }

    public void   setOccupancy(double occu){ occupancy = occu ;} ;
    public double getOccupancy(){ return occupancy; } ;

    public void   setTempFactor(double temp){ tempfactor = temp ;} ;
    public double getTempFactor(){ return tempfactor; } ;

    /** returns and identical copy of this  object .
     * @return  and identical copy of this  object
     */
    public Object clone() {
        AtomImpl n = new AtomImpl();
        n.setOccupancy(getOccupancy());
        n.setTempFactor(getTempFactor());
        n.setAltLoc(getAltLoc());
        double[] coords = getCoords();
        n.setX(coords[0]);
        n.setY(coords[1]);
        n.setZ(coords[2]);
        n.setPDBserial(getPDBserial());
        n.setFullName(getFullName());
        n.setName(getName());
        n.setElement(getElement());
        return n ;
    }

    public void setParent(Group parent) {

       setGroup(parent);
    }

    public Group getParent() {
        return getGroup();
    }
    
    
    public void setGroup(Group parent){
    	this.parent = parent;
    }
    
    public Group getGroup(){
    	return parent;
    }
	public Element getElement() {
		return element;
	}
	
	public void setElement(Element e) {
		this.element = e;
	
	}
	
	public String toPDB() {
		
		return FileConvert.toPDB(this);
	}
	public void toPDB(StringBuffer buf) {
		FileConvert.toPDB(this,buf);
		
	}


}
