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
package org.biojava.nbio.structure.contact;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;


/**
 * A grid to be used for calculating atom contacts through a spatial hashing algorithm.
 * <p>
 * The grid is composed of cells of size of the cutoff so that the distances that need to be calculated
 * are reduced to those within each cell and to the neighbouring cells.
 * <p>
 * Usage, for generic 3D points:
 * <pre>
 *  Point3d[] points = ...;
 *  Grid grid = new Grid(8.0);
 *  grid.addCoords(points);
 *  List&lt;Contact&gt; contacts = getIndicesContacts();
 * </pre>
 * Usage, for atoms:
 * <pre>
 *  Atom[] atoms = ...;
 *  Grid grid = new Grid(8.0);
 *  grid.addCoords(atoms);
 *  AtomContactSet contacts = getAtomContacts();
 * </pre> 
 *
 * @author Jose Duarte
 *
 */
public class Grid {

	/**
	 * The scale: we use units of hundredths of Angstroms (thus cutoffs can be specified with a maximum precision of 0.01A)
	 */
	private static final int SCALE=100;

	private GridCell[][][] cells;

	private double cutoff;
	private int cellSize;

	private Point3d[] iAtoms;
	private Point3d[] jAtoms;
	
	private Atom[] iAtomObjects;
	private Atom[] jAtomObjects;

	// the bounds in int grid coordinates
	private int[] bounds;

	// the i and j bounding boxes in original double coordinates
	private BoundingBox ibounds;
	private BoundingBox jbounds;

	private boolean noOverlap; // if the 2 sets of atoms are found not to overlap then this is set to true

	/**
	 * Creates a <code>Grid</code>, the cutoff is in the same units as the coordinates
	 * (Angstroms if they are atom coordinates) and can
	 * be specified to a precision of 0.01.
	 * @param cutoff
	 */
	public Grid(double cutoff) {
		this.cutoff = cutoff;
		this.cellSize = (int) Math.floor(cutoff*SCALE);
		this.noOverlap = false;
	}

	private int getFloor(double number) {
		return (cellSize*((int)Math.floor(number*SCALE/cellSize)));
	}

	private int xintgrid2xgridindex(int xgridDim) {
		return (xgridDim-bounds[0])/cellSize;
	}

	private int yintgrid2ygridindex(int ygridDim) {
		return (ygridDim-bounds[1])/cellSize;
	}

	private int zintgrid2zgridindex(int zgridDim) {
		return (zgridDim-bounds[2])/cellSize;
	}

	/**
	 * Adds the i and j atoms and fills the grid. Their bounds will be computed.
	 * Subsequent call to {@link #getIndicesContacts()} or {@link #getAtomContacts()} will produce the interatomic contacts.
	 * @param iAtoms
	 * @param jAtoms
	 */
	public void addAtoms(Atom[] iAtoms, Atom[] jAtoms) {
		addAtoms(iAtoms, null, jAtoms, null);
	}

	/**
	 * Adds the i and j atoms and fills the grid, passing their bounds (array of size 6 with x,y,z minima and x,y,z maxima)
	 * This way the bounds don't need to be recomputed.
	 * Subsequent call to {@link #getIndicesContacts()} or {@link #getAtomContacts()} will produce the interatomic contacts.
	 * @param iAtoms
	 * @param icoordbounds
	 * @param jAtoms
	 * @param jcoordbounds
	 */
	public void addAtoms(Atom[] iAtoms, BoundingBox icoordbounds, Atom[] jAtoms, BoundingBox jcoordbounds) {
		this.iAtoms = Calc.atomsToPoints(iAtoms);
		this.iAtomObjects = iAtoms;

		if (icoordbounds!=null) {
			this.ibounds = icoordbounds;
		} else {
			this.ibounds = new BoundingBox(this.iAtoms);
		}

		this.jAtoms = Calc.atomsToPoints(jAtoms);
		this.jAtomObjects = jAtoms;

		if (jAtoms==iAtoms) {
			this.jbounds=ibounds;
		} else {
			if (jcoordbounds!=null) {
				this.jbounds = jcoordbounds;
			} else {
				this.jbounds = new BoundingBox(this.jAtoms);

			}
		}

		fillGrid();
	}

	/**
	 * Adds a set of atoms, subsequent call to {@link #getIndicesContacts()} or {@link #getAtomContacts()} will produce the interatomic contacts.
	 * The bounding box of the atoms will be computed based on input array
	 * @param atoms
	 */
	public void addAtoms(Atom[] atoms) {
		addAtoms(atoms, (BoundingBox) null);
	}

	/**
	 * Adds a set of atoms, subsequent call to {@link #getIndicesContacts()} or {@link #getAtomContacts()} will produce the interatomic contacts.
	 * The bounds calculated elsewhere can be passed, or if null they are computed.
	 * @param atoms
	 * @param bounds
	 */
	public void addAtoms(Atom[] atoms, BoundingBox bounds) {
		this.iAtoms = Calc.atomsToPoints(atoms);
		this.iAtomObjects = atoms;

		if (bounds!=null) {
			this.ibounds = bounds;
		} else {
			this.ibounds = new BoundingBox(iAtoms);
		}

		this.jAtoms = null;
		this.jAtomObjects = null;
		this.jbounds = null;

		fillGrid();
	}
	
	/**
	 * Adds the i and j coordinates and fills the grid. Their bounds will be computed.
	 * Subsequent call to {@link #getIndicesContacts()} will produce the 
	 * contacts, i.e. the set of points within distance cutoff.
	 * 
	 * Subsequent calls to method {@link #getAtomContacts()} will produce a NullPointerException 
	 * since this only adds coordinates and no atom information.
	 * @param iAtoms
	 * @param jAtoms
	 */
	public void addCoords(Point3d[] iAtoms, Point3d[] jAtoms) {
		addCoords(iAtoms, null, jAtoms, null);
	}

	/**
	 * Adds the i and j coordinates and fills the grid, passing their bounds (array of size 6 with x,y,z minima and x,y,z maxima)
	 * This way the bounds don't need to be recomputed.
	 * Subsequent call to {@link #getIndicesContacts()} will produce the 
	 * contacts, i.e. the set of points within distance cutoff.
	 * 
	 * Subsequent calls to method {@link #getAtomContacts()} will produce a NullPointerException 
	 * since this only adds coordinates and no atom information.
	 * @param iAtoms
	 * @param icoordbounds
	 * @param jAtoms
	 * @param jcoordbounds
	 */
	public void addCoords(Point3d[] iAtoms, BoundingBox icoordbounds, Point3d[] jAtoms, BoundingBox jcoordbounds) {
		this.iAtoms = iAtoms;
		this.iAtomObjects = null;

		if (icoordbounds!=null) {
			this.ibounds = icoordbounds;
		} else {
			this.ibounds = new BoundingBox(this.iAtoms);
		}

		this.jAtoms = jAtoms;
		this.jAtomObjects = null;

		if (jAtoms==iAtoms) {
			this.jbounds=ibounds;
		} else {
			if (jcoordbounds!=null) {
				this.jbounds = jcoordbounds;
			} else {
				this.jbounds = new BoundingBox(this.jAtoms);

			}
		}

		fillGrid();
	}

	/**
	 * Adds a set of coordinates, subsequent call to {@link #getIndicesContacts()} will produce the 
	 * contacts, i.e. the set of points within distance cutoff.
	 * The bounding box of the atoms will be computed based on input array.
	 * Subsequent calls to method {@link #getAtomContacts()} will produce a NullPointerException 
	 * since this only adds coordinates and no atom information.
	 * @param atoms
	 */
	public void addCoords(Point3d[] atoms) {
		addCoords(atoms, (BoundingBox) null);
	}

	/**
	 * Adds a set of coordinates, subsequent call to {@link #getIndicesContacts()} will produce the 
	 * contacts, i.e. the set of points within distance cutoff.
	 * The bounds calculated elsewhere can be passed, or if null they are computed.
	 * Subsequent calls to method {@link #getAtomContacts()} will produce a NullPointerException 
	 * since this only adds coordinates and no atom information.
	 * @param atoms
	 * @param bounds
	 */
	public void addCoords(Point3d[] atoms, BoundingBox bounds) {
		this.iAtoms = atoms;
		this.iAtomObjects = null;

		if (bounds!=null) {
			this.ibounds = bounds;
		} else {
			this.ibounds = new BoundingBox(iAtoms);
		}

		this.jAtoms = null;
		this.jAtomObjects = null;
		this.jbounds = null;

		fillGrid();
	}

	/**
	 * Creates the grid based on the boundaries defined by all atoms given (iAtoms and jAtoms)
	 * and places the atoms in their corresponding grid cells.
	 * Checks also if the i and j grid overlap, i.e. the enclosing bounds of
	 * the 2 grids (i and j) are no more than one cell size apart. If they don't
	 * overlap then they are too far apart so there's nothing to calculate, we set
	 * the noOverlap flag and then {@link #getIndicesContacts()} will do no calculation at all.
	 */
	private void fillGrid() {

		if (jbounds!=null && !ibounds.overlaps(jbounds, cutoff)) {
			//System.out.print("-");
			noOverlap = true;
			return;
		}

		findFullGridIntBounds();

		cells = new GridCell[1+(bounds[3]-bounds[0])/cellSize]
		                    [1+(bounds[4]-bounds[1])/cellSize]
		                    [1+(bounds[5]-bounds[2])/cellSize];

		int i = 0;
		for (Point3d atom:iAtoms) {

			int xind = xintgrid2xgridindex(getFloor(atom.x));
			int yind = yintgrid2ygridindex(getFloor(atom.y));
			int zind = zintgrid2zgridindex(getFloor(atom.z));
			if (cells[xind][yind][zind]==null) {
				cells[xind][yind][zind] = new GridCell(this);
			}
			cells[xind][yind][zind].addIindex(i);
			i++;
		}

		if (jAtoms==null) return;

		int j = 0;
		for (Point3d atom:jAtoms) {

			int xind = xintgrid2xgridindex(getFloor(atom.x));
			int yind = yintgrid2ygridindex(getFloor(atom.y));
			int zind = zintgrid2zgridindex(getFloor(atom.z));
			if (cells[xind][yind][zind]==null) {
				cells[xind][yind][zind] = new GridCell(this);
			}
			cells[xind][yind][zind].addJindex(j);
			j++;
		}

	}

	/**
	 * Calculates an int array of size 6 into member variable bounds:
	 * - elements 0,1,2: minimum x,y,z of the iAtoms and jAtoms
	 * - elements 3,4,5: maximum x,y,z of the iAtoms and jAtoms
	 */
	private void findFullGridIntBounds() {
		int[] iIntBounds = getIntBounds(ibounds);

		bounds = new int[6];
		if (jbounds==null) {
			bounds = iIntBounds;
		} else {
			int[] jIntBounds = getIntBounds(jbounds);
			bounds[0] = Math.min(iIntBounds[0],jIntBounds[0]);
			bounds[1] = Math.min(iIntBounds[1],jIntBounds[1]);
			bounds[2] = Math.min(iIntBounds[2],jIntBounds[2]);
			bounds[3] = Math.max(iIntBounds[3],jIntBounds[3]);
			bounds[4] = Math.max(iIntBounds[4],jIntBounds[4]);
			bounds[5] = Math.max(iIntBounds[5],jIntBounds[5]);
		}
	}

	/**
	 * Returns an int array of size 6 :
	 * - elements 0,1,2: minimum x,y,z (in grid int coordinates) of the given atoms
	 * - elements 3,4,5: maximum x,y,z (in grid int coordinates) of the given atoms
	 * @return
	 */
	private int[] getIntBounds(BoundingBox coordbounds) {
		int[] bs = new int[6];
		bs[0] = getFloor(coordbounds.xmin);
		bs[1] = getFloor(coordbounds.ymin);
		bs[2] = getFloor(coordbounds.zmin);
		bs[3] = getFloor(coordbounds.xmax);
		bs[4] = getFloor(coordbounds.ymax);
		bs[5] = getFloor(coordbounds.zmax);
		return bs;
	}

	/**
	 * Returns all contacts, i.e. all atoms that are within the cutoff distance.
	 * If both iAtoms and jAtoms are defined then contacts are between iAtoms and jAtoms,
	 * if jAtoms is null, then contacts are within the iAtoms.
	 * @return
	 */
	public AtomContactSet getAtomContacts() {

		AtomContactSet contacts = new AtomContactSet(cutoff);

		List<Contact> list = getIndicesContacts();

		if (jAtomObjects == null) {
			for (Contact cont : list) {
				contacts.add(new AtomContact(new Pair<Atom>(iAtomObjects[cont.getI()],iAtomObjects[cont.getJ()]),cont.getDistance()));
			}

		} else {
			for (Contact cont : list) {
				contacts.add(new AtomContact(new Pair<Atom>(iAtomObjects[cont.getI()],jAtomObjects[cont.getJ()]),cont.getDistance()));
			}
		}

		return contacts;
	}
	
	/**
	 * Returns all contacts, i.e. all atoms that are within the cutoff distance.
	 * If both iAtoms and jAtoms are defined then contacts are between iAtoms and jAtoms,
	 * if jAtoms is null, then contacts are within the iAtoms.
	 * @return
	 * @deprecated use {@link #getAtomContacts()} instead
	 */	
	@Deprecated
	public AtomContactSet getContacts() {
		return getAtomContacts();
	}
	
	/**
	 * Returns all contacts, i.e. all atoms that are within the cutoff distance, as simple Contact objects containing the atom indices pairs and the distance.
	 * If both iAtoms and jAtoms are defined then contacts are between iAtoms and jAtoms,
	 * if jAtoms is null, then contacts are within the iAtoms.
	 * @return
	 */
	public List<Contact> getIndicesContacts() {

		List<Contact> list = new ArrayList<>();

		// if the 2 sets of atoms are not overlapping they are too far away and no need to calculate anything
		// this won't apply if there's only one set of atoms (iAtoms), where we would want all-to-all contacts
		if (noOverlap) return list;


		for (int xind=0;xind<cells.length;xind++) {
			for (int yind=0;yind<cells[xind].length;yind++) {
				for (int zind=0;zind<cells[xind][yind].length;zind++) {
					// distances of points within this cell
					GridCell thisCell = cells[xind][yind][zind];
					if (thisCell==null) continue;

					list.addAll(thisCell.getContactsWithinCell());

					// distances of points from this box to all neighbouring boxes: 26 iterations (26 neighbouring boxes)
					for (int x=xind-1;x<=xind+1;x++) {
						for (int y=yind-1;y<=yind+1;y++) {
							for (int z=zind-1;z<=zind+1;z++) {
								if (x==xind && y==yind && z==zind) continue;

								if (x>=0 && x<cells.length && y>=0 && y<cells[x].length && z>=0 && z<cells[x][y].length) {
									if (cells[x][y][z] == null) continue;

									list.addAll(thisCell.getContactsToOtherCell(cells[x][y][z]));
								}
							}
						}
					}
				}
			}
		}

		return list;
	}

	/**
	 * Fast determination of whether any atoms from a given set fall within
	 * the cutoff of iAtoms. If {@link #addAtoms(Atom[], Atom[])} was called
	 * with two sets of atoms, contacts to either set are considered.
	 * @param atoms
	 * @return
	 */
	public boolean hasAnyContact(Point3d[] atoms) {
		return hasAnyContact(Arrays.asList(atoms));
	}
	public boolean hasAnyContact(Collection<Point3d> atoms) {
		for(Point3d atom : atoms) {
			// Calculate Grid cell for the atom
			int xind = xintgrid2xgridindex(getFloor(atom.x));
			int yind = yintgrid2ygridindex(getFloor(atom.y));
			int zind = zintgrid2zgridindex(getFloor(atom.z));

			// Consider 3x3x3 grid of cells around point
			for (int x=xind-1;x<=xind+1;x++) {
				if( x<0 || cells.length<=x) continue;
				for (int y=yind-1;y<=yind+1;y++) {
					if( y<0 || cells[x].length<=y ) continue;
					for (int z=zind-1;z<=zind+1;z++) {
						if( z<0 || cells[x][y].length<=z ) continue;
						
						GridCell cell = cells[x][y][z];
						// Check for contacts in this cell
						if(cell != null && cell.hasContactToAtom(iAtoms, jAtoms, atom, cutoff)) {
							return true;
						}
					}
				}
			}
		}
		return false;
	}

	public double getCutoff() {
		return cutoff;
	}

	/**
	 * Tells whether (after having added atoms to grid) the i and j grids are not overlapping.
	 * Overlap is defined as enclosing bounds of the 2 grids being no more than one cell size apart.
	 * @return true if the 2 grids don't overlap, false if they do
	 */
	public boolean isNoOverlap() {
		return noOverlap;
	}
	
	protected Point3d[] getIAtoms() {
		return iAtoms;
	}
	
	protected Point3d[] getJAtoms() {
		return jAtoms;
	}

}
