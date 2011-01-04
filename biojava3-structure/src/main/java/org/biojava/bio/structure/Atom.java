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

/**
 * A simple interface for an Atom.
 * The coordinates can be accessed via the
 * {@link #getCoords()}, or the 
 * {@link #getX()}, {@link #getY()}, {@link #getZ()} methods.
 * There are a few additional methods here to provide some PDB specific information.
 * 
 
 * @author Andreas Prlic
 * @author Horvath Tamas
 * @version %I% %G%
 * @since 1.4
 * 
 */
public interface Atom extends Cloneable, PDBRecord {
    
    /** set trimmed version of atom name, e.g. "CA". 
     * @param s  a String specifying the name value
     * @see #getName
     */
    public void   setName(String s);
    
    /** get trimmed version of atom name, e.g. "CA".
     * @return a String representing the name value 
     * @see #setName
     */
    public String getName();
    
    /** set full name of atom e.g. " CA ". 
     * @param s  a String specifying the full name value
     * @see #getFullName
     */
    public void   setFullName(String s) ;
    
    /** get full name of atom e.g. " CA ".
     * @return a String representing the full name value 
     * @see #setFullName
     */
    public String getFullName();
    
    /** set element of the atom name, e.g. Element.Fe
     * @param e  an Element enumeration
     * @see #getElement
     */
    public void   setElement(Element e);
    
    /** get element of the atom, e.g. Element.Ca
     * @return an Element enumeration 
     * @see #setElement
     */
    public Element getElement();
    
    /** set PDB atom number. 
     * @param i  an int specifying the PDBserial value 
     * @see #getPDBserial
     */
    public void setPDBserial(int i) ;
    
    /** get PDB atom number.
     * @return an int representing the PDBserial value 
     * @see #setPDBserial
     */
    public int  getPDBserial() ;
    
    /** set the coordinates.  
     * @param c  an array of doubles specifying the coords value
     * @see #getCoords
     */    
    public void    setCoords(double[] c);
    
    /** get the coordinates. 
     * @return an array of doubles representing the coords value
     * @see #setCoords
     */    
    public double[] getCoords() ;
    
    /** Set the X coordinate.
     * @param x  a double
     * @see #getX()
     */
    public void setX(double x);
    
    /** Set the Y coordinate.
     * @param y  a double
     * @see #getY()
     */
    public void setY(double y);
    
    /** Set the Z coordinate.
     * @param z  a double
     * @see #getZ()     
     */
    public void setZ(double z);
    
    /** Get coordinate X. 
     * @return a double
     * @see #setX(double)
     */    
    public double getX() ;
    
    /** Get coordinate Y. 
     * @return a double
     * @see #setY(double)
     */    
    public double getY() ;
    
    /** Get coordinate Z. 
     * @return a double
     * @see #setZ(double)
     */    
    public double getZ() ;
    
    /** get set alternate Location.
     * @param c  a Character object specifying the alt loc value 
     * @see #getAltLoc
     */
    public void setAltLoc(Character c);
    
    /** get set alternate Location. 
     * @return a Character object representing the alt loc value
     * @see #setAltLoc
     */
    public Character getAltLoc();
    
    /** store the whole line.  
     * @param s  a String specifying the PDBline value
     * @see #getPDBline
     * @deprecated replaced by {@link #toPDB()}
     */
    public void   setPDBline(String s) ;
    
    /** store the whole line.  
     * @return a String representing the PDBline value
     * @see #setPDBline
     * @deprecated @deprecated replaced by {@link #toPDB()}
     */
    
    public String getPDBline() ;
    
    /** set occupancy. 
     * @param occupancy  a double specifying the occupancy value
     * @see #getOccupancy
     */
    public void setOccupancy(double occupancy) ;
    
    /** get occupancy. 
     * @return a double representing the occupancy value
     * @see #setOccupancy
     */
    public double getOccupancy();
    
    /** get set temp factor .
     * @param temp  a double specifying the temp factor value
     * @see #getTempFactor
     */
    public void   setTempFactor(double temp) ;
    
    /** get set temp factor. 
     * @return a double representing the temp factor value
     * @see #setTempFactor
     */
    public double getTempFactor() ;
    
    /** returns and identical copy of this  object .
     * @return  and identical copy of this  object 
     */
    public Object clone();
    
    /** Sets the back-reference to its parent Group.
     * @param parent the parent Group
     * @see #getParent()
     * @deprecated replaced by {@link #setGroup(Group)}
     */
    public void setParent(Group parent) ; 
    
    /** Sets the back-reference to its parent Group.
     * @param parent the parent Group
     * @see #getGroup()
     */
    
    public void setGroup(Group parent);
    
     /** Returns the parent Group of the Atom.
     * returns null if the referenced object is not Group 
     * @return Group the parent Group of the Atom, or null
     * @see #setParent(Group)   
     * @deprecated replaced by {@link #getGroup()}  
     */
    public Group getParent();
    
    /** Returns the parent Group of the Atom.
     * returns null if the referenced object is not Group 
     * @return Group the parent Group of the Atom, or null
     * @see #setParent(Group)     
     */
    public Group getGroup();
    
    
}
