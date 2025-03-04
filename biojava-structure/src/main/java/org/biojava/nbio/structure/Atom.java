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
package org.biojava.nbio.structure;

import java.util.List;

import javax.vecmath.Point3d;

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

	/**
	 * Set atom name, e.g. "CA".
	 * @param s  a trimmed String specifying the name value
	 * @see #getName
	 */
	void   setName(String s);

	/**
	 * Get atom name, e.g. "CA".
	 * Beware that some PDB atom names are ambiguous (e.g. CA, which means C-alpha or Calcium),
	 * the ambiguity can simply be resolved by also checking the element with {@link #getElement()}
	 * @return a trimmed String representing the name value
	 * @see #setName
	 */
	String getName();

	/**
	 * Set element of the atom name, e.g. {@link Element#Fe}
	 * @param e  an Element enumeration
	 * @see #getElement
	 */
	void   setElement(Element e);

	/**
	 * Get element of the atom, e.g. {@link Element#Ca}
	 * @return an Element enumeration
	 * @see #setElement
	 */
	Element getElement();

	/**
	 * Set PDB atom number.
	 * @param i  an int specifying the PDBserial value
	 * @see #getPDBserial
	 */
	void setPDBserial(int i) ;

	/**
	 * Get PDB atom number.
	 * @return an int representing the PDBserial value
	 * @see #setPDBserial
	 */
	int  getPDBserial() ;

	/**
	 * Set the coordinates.
	 * @param c  an array of doubles specifying the coords value
	 * @see #getCoords
	 */
	void    setCoords(double[] c);

	/**
	 * Get the coordinates.
	 * @return a new array of doubles representing the coords value
	 * @see #setCoords
	 * @see #getCoordsAsPoint3d()
	 */
	double[] getCoords() ;

	/**
	 * Get the coordinates.
	 * <p>
	 * Internally the coordinates are represented as Point3d so this
	 * is recommended over {@link #getCoords()}
	 * @return a reference to the Point3d coordinates
	 * @see #getCoords()
	 */
	Point3d getCoordsAsPoint3d();

	/**
	 * Set the X coordinate.
	 * @param x  a double
	 * @see #getX()
	 */
	void setX(double x);

	/**
	 * Set the Y coordinate.
	 * @param y  a double
	 * @see #getY()
	 */
	void setY(double y);

	/**
	 * Set the Z coordinate.
	 * @param z  a double
	 * @see #getZ()
	 */
	void setZ(double z);

	/**
	 * Get coordinate X.
	 * @return a double
	 * @see #setX(double)
	 */
	double getX() ;

	/**
	 * Get coordinate Y.
	 * @return a double
	 * @see #setY(double)
	 */
	double getY() ;

	/**
	 * Get coordinate Z.
	 * @return a double
	 * @see #setZ(double)
	 */
	double getZ() ;

	/**
	 * Set alternate Location.
	 * @param c  a Character object specifying the alt loc value
	 * @see #getAltLoc
	 */
	void setAltLoc(Character c);

	/**
	 * Get alternate Location.
	 * @return a Character object representing the alt loc value. Default altLoc ('.' in mmCIF files)
	 * is represented by ' ' (space character, ascii 32).
	 * @see #setAltLoc
	 */
	Character getAltLoc();

	/**
	 * Set occupancy.
	 * @param occupancy  a float specifying the occupancy value
	 * @see #getOccupancy
	 */
	void setOccupancy(float occupancy) ;

	/**
	 * Get occupancy.
	 * @return a float representing the occupancy value
	 * @see #setOccupancy
	 */
	float getOccupancy();

	/**
	 * Set temp factor .
	 * @param temp  a float specifying the temp factor value
	 * @see #getTempFactor
	 */
	void   setTempFactor(float temp) ;

	/**
	 * Get temp factor.
	 * @return a float representing the temp factor value
	 * @see #setTempFactor
	 */
	float getTempFactor() ;

	/**
	 * Return an identical copy of this  object .
	 * @return  an identical copy of this  object
	 */
	Object clone();

	/**
	 * Set the back-reference to its parent Group.
	 * @param parent the parent Group
	 * @see #getGroup()
	 */

	void setGroup(Group parent);

	/**
	 * Return the parent Group of the Atom.
	 * returns null if the referenced object is not Group
	 * @return Group the parent Group of the Atom, or null
	 * @see #setGroup(Group)
	 */
	Group getGroup();

	/**
	 * Add a bond
	 * @param bond to be added
	 * @see #getBonds()
	 */
	void addBond(Bond bond);

	/**
	 * Get all {@link Bond}s this atom is part of.
	 *
	 * @return a list of {@link Bond}s or null if no bonds exist for this Atom
	 */
	List<Bond> getBonds();

	/**
	 * Sets the bonds
	 * @param bonds
	 */
	void setBonds(List<Bond> bonds);

	/**
	 * Test if another atom has a bond to this atom
	 *
	 * @param other
	 * @return
	 */
	boolean hasBond(Atom other);

	/**
	 * Get the charge of this atom
	 *
	 * @return a the integer charge.
	 */
	short getCharge();

	/**
	 * Set the charge of this atom
	 *
	 */
	void setCharge(short charge);
}
