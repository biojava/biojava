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
package org.biojava.nbio.structure.asa;

import org.biojava.nbio.structure.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;




/**
 * Class to calculate Accessible Surface Areas based on
 * the rolling ball algorithm by Shrake and Rupley.
 *
 * The code is adapted from a python implementation at http://boscoh.com/protein/asapy
 * (now source is available at https://github.com/boscoh/asa).
 * Thanks to Bosco K. Ho for a great piece of code and for his fantastic blog.
 *
 * See
 * Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms.
 * Lysozyme and Insulin." JMB (1973) 79:351-371.
 * Lee, B., and Richards, F.M. "The interpretation of Protein Structures: Estimation of
 * Static Accessibility" JMB (1971) 55:379-400
 * @author duarte_j
 *
 */
public class AsaCalculator {

	private static final Logger logger = LoggerFactory.getLogger(AsaCalculator.class);

	// Bosco uses as default 960, Shrake and Rupley seem to use in their paper 92 (not sure if this is actually the same parameter)
	public static final int DEFAULT_N_SPHERE_POINTS = 960;
	public static final double DEFAULT_PROBE_SIZE = 1.4;
	public static final int DEFAULT_NTHREADS = 1;



	// Chothia's amino acid atoms vdw radii
	public static final double TRIGONAL_CARBON_VDW = 1.76;
	public static final double TETRAHEDRAL_CARBON_VDW = 1.87;
	public static final double TRIGONAL_NITROGEN_VDW = 1.65;
	public static final double TETRAHEDRAL_NITROGEN_VDW = 1.50;
	public static final double SULFUR_VDW = 1.85;
	public static final double OXIGEN_VDW = 1.40;

	// Chothia's nucleotide atoms vdw radii
	public static final double NUC_CARBON_VDW = 1.80;
	public static final double NUC_NITROGEN_VDW = 1.60;
	public static final double PHOSPHOROUS_VDW = 1.90;





	private class AsaCalcWorker implements Runnable {

		private int i;
		private double[] asas;

		public AsaCalcWorker(int i, double[] asas) {
			this.i = i;
			this.asas = asas;
		}

		@Override
		public void run() {
			asas[i] = calcSingleAsa(i);
		}
	}


	private Point3d[] atomCoords;
	private Atom[] atoms;
	private double[] radii;
	private double probe;
	private int nThreads;
	private Point3d[] spherePoints;
	private double cons;

	/**
	 * Constructs a new AsaCalculator. Subsequently call {@link #calculateAsas()}
	 * or {@link #getGroupAsas()} to calculate the ASAs
	 * Only non-Hydrogen atoms are considered in the calculation.
	 * @param structure
	 * @param probe
	 * @param nSpherePoints
	 * @param nThreads
	 * @param hetAtoms if true HET residues are considered, if false they aren't, equivalent to
	 * NACCESS' -h option
	 * @see StructureTools.getAllNonHAtomArray
	 */
	public AsaCalculator(Structure structure, double probe, int nSpherePoints, int nThreads, boolean hetAtoms) {
		this.atoms = StructureTools.getAllNonHAtomArray(structure, hetAtoms);
		this.atomCoords = Calc.atomsToPoints(atoms);
		this.probe = probe;
		this.nThreads = nThreads;

		// initialising the radii by looking them up through AtomRadii
		radii = new double[atomCoords.length];
		for (int i=0;i<atomCoords.length;i++) {
			radii[i] = getRadius(atoms[i]);
		}

		// initialising the sphere points to sample
		spherePoints = generateSpherePoints(nSpherePoints);

		cons = 4.0 * Math.PI / nSpherePoints;
	}

	/**
	 * Constructs a new AsaCalculator. Subsequently call {@link #calculateAsas()}
	 * or {@link #getGroupAsas()} to calculate the ASAs.
	 * @param atoms an array of atoms not containing Hydrogen atoms
	 * @param probe the probe size
	 * @param nSpherePoints the number of points to be used in generating the spherical
	 * dot-density, the more points the more accurate (and slower) calculation
	 * @param nThreads the number of parallel threads to use for the calculation
	 * @throws IllegalArgumentException if any atom in the array is a Hydrogen atom
	 */
	public AsaCalculator(Atom[] atoms, double probe, int nSpherePoints, int nThreads) {
		this.atoms = atoms;
		this.atomCoords = Calc.atomsToPoints(atoms);
		this.probe = probe;
		this.nThreads = nThreads;

		for (Atom atom:atoms) {
			if (atom.getElement()==Element.H)
				throw new IllegalArgumentException("Can't calculate ASA for an array that contains Hydrogen atoms ");
		}

		// initialising the radii by looking them up through AtomRadii
		radii = new double[atoms.length];
		for (int i=0;i<atoms.length;i++) {
			radii[i] = getRadius(atoms[i]);
		}

		// initialising the sphere points to sample
		spherePoints = generateSpherePoints(nSpherePoints);

		cons = 4.0 * Math.PI / nSpherePoints;
	}

	/**
	 * Constructs a new AsaCalculator. Subsequently call {@link #calcSingleAsa(int)} 
	 * to calculate the atom ASAs. The given radius parameter will be taken as the radius for
	 * all points given. No ASA calculation per group will be possible with this constructor, so
	 * usage of {@link #getGroupAsas()} will result in a NullPointerException.
	 * @param atomCoords
	 * 				the coordinates representing the center of atoms
	 * @param probe 
	 * 				the probe size
	 * @param nSpherePoints 
	 * 				the number of points to be used in generating the spherical
	 * 				dot-density, the more points the more accurate (and slower) calculation
	 * @param nThreads 
	 * 				the number of parallel threads to use for the calculation
	 * @param radius 
	 * 				the radius that will be assign to all given coordinates
	 */
	public AsaCalculator(Point3d[] atomCoords, double probe, int nSpherePoints, int nThreads, double radius) {
		this.atoms = null;
		this.atomCoords = atomCoords;
		this.probe = probe;
		this.nThreads = nThreads;

		// initialising the radii to the given radius for all atoms
		radii = new double[atomCoords.length];
		for (int i=0;i<atomCoords.length;i++) {
			radii[i] = radius;
		}

		// initialising the sphere points to sample
		spherePoints = generateSpherePoints(nSpherePoints);

		cons = 4.0 * Math.PI / nSpherePoints;
	}

	/**
	 * Calculates ASA for all atoms and return them as a GroupAsa
	 * array (one element per residue in structure) containing ASAs per residue
	 * and per atom.
	 * The sorting of Groups in returned array is as specified by {@link org.biojava.nbio.structure.ResidueNumber}
	 * @return
	 */
	public GroupAsa[] getGroupAsas() {

		TreeMap<ResidueNumber, GroupAsa> asas = new TreeMap<ResidueNumber, GroupAsa>();

		double[] asasPerAtom = calculateAsas();

		for (int i=0;i<atomCoords.length;i++) {
			Group g = atoms[i].getGroup();
			if (!asas.containsKey(g.getResidueNumber())) {
				GroupAsa groupAsa = new GroupAsa(g);
				groupAsa.addAtomAsaU(asasPerAtom[i]);
				asas.put(g.getResidueNumber(), groupAsa);
			} else {
				GroupAsa groupAsa = asas.get(g.getResidueNumber());
				groupAsa.addAtomAsaU(asasPerAtom[i]);
			}
		}

		return asas.values().toArray(new GroupAsa[asas.size()]);
	}

	/**
	 * Calculates the Accessible Surface Areas for the atoms given in constructor and with parameters given.
	 * Beware that the parallel implementation is quite memory hungry. It scales well as long as there is
	 * enough memory available.
	 * @return an array with asa values corresponding to each atom of the input array
	 */
	public double[] calculateAsas() {

		double[] asas = new double[atomCoords.length];

		if (nThreads<=1) { // (i.e. it will also be 1 thread if 0 or negative number specified)
			for (int i=0;i<atomCoords.length;i++) {
				asas[i] = calcSingleAsa(i);
			}

		} else {
			// NOTE the multithreaded calculation does not scale up well in some systems,
			// why? I guess some memory/garbage collect problem? I tried increasing Xmx in pc8201 but didn't help

			// Following scaling tests are for 3hbx, calculating ASA of full asym unit (6 chains):

			// SCALING test done in merlinl01 (12 cores, Xeon X5670  @ 2.93GHz, 24GB RAM)
			//1 threads, time:  8.8s -- x1.0
			//2 threads, time:  4.4s -- x2.0
			//3 threads, time:  2.9s -- x3.0
			//4 threads, time:  2.2s -- x3.9
			//5 threads, time:  1.8s -- x4.9
			//6 threads, time:  1.6s -- x5.5
			//7 threads, time:  1.4s -- x6.5
			//8 threads, time:  1.3s -- x6.9

			// SCALING test done in pc8201 (4 cores, Core2 Quad Q9550  @ 2.83GHz, 8GB RAM)
			//1 threads, time: 17.2s -- x1.0
			//2 threads, time:  9.7s -- x1.8
			//3 threads, time:  7.7s -- x2.2
			//4 threads, time:  7.9s -- x2.2

			// SCALING test done in eppic01 (16 cores, Xeon E5-2650 0  @ 2.00GHz, 128GB RAM)
			//1 threads, time: 10.7s -- x1.0
			//2 threads, time:  5.6s -- x1.9
			//3 threads, time:  3.6s -- x3.0
			//4 threads, time:  2.8s -- x3.9
			//5 threads, time:  2.3s -- x4.8
			//6 threads, time:  1.8s -- x6.0
			//7 threads, time:  1.6s -- x6.8
			//8 threads, time:  1.3s -- x8.0
			//9 threads, time:  1.3s -- x8.5
			//10 threads, time:  1.1s -- x10.0
			//11 threads, time:  1.0s -- x10.9
			//12 threads, time:  0.9s -- x11.4



			ExecutorService threadPool = Executors.newFixedThreadPool(nThreads);


			for (int i=0;i<atomCoords.length;i++) {
				threadPool.submit(new AsaCalcWorker(i,asas));
			}

			threadPool.shutdown();

			while (!threadPool.isTerminated());

		}

		return asas;
	}

	/**
	 * Returns list of 3d coordinates of points on a sphere using the
	 * Golden Section Spiral algorithm.
	 * @param nSpherePoints the number of points to be used in generating the spherical dot-density
	 * @return
	 */
	private Point3d[] generateSpherePoints(int nSpherePoints) {
		Point3d[] points = new Point3d[nSpherePoints];
		double inc = Math.PI * (3.0 - Math.sqrt(5.0));
		double offset = 2.0 / nSpherePoints;
		for (int k=0;k<nSpherePoints;k++) {
			double y = k * offset - 1.0 + (offset / 2.0);
			double r = Math.sqrt(1.0 - y*y);
			double phi = k * inc;
			points[k] = new Point3d(Math.cos(phi)*r, y, Math.sin(phi)*r);
		}
		return points;
	}

	/**
	 * Returns list of indices of atoms within probe distance to atom k.
	 * @param k index of atom for which we want neighbor indices
	 */
	private ArrayList<Integer> findNeighborIndices(int k) {
		// looking at a typical protein case, number of neighbours are from ~10 to ~50, with an average of ~30
		// Thus 40 seems to be a good compromise for the starting capacity
		ArrayList<Integer> neighbor_indices = new ArrayList<Integer>(40);

		double radius = radii[k] + probe + probe;

		for (int i=0;i<atomCoords.length;i++) {
			if (i==k) continue;

			double dist = 0;

			dist = atomCoords[i].distance(atomCoords[k]);

			if (dist < radius + radii[i]) {
				neighbor_indices.add(i);
			}

		}

		return neighbor_indices;
	}

	private double calcSingleAsa(int i) {
		Point3d atom_i = atomCoords[i];
		ArrayList<Integer> neighbor_indices = findNeighborIndices(i);
		int n_neighbor = neighbor_indices.size();
		int j_closest_neighbor = 0;
		double radius = probe + radii[i];

		int n_accessible_point = 0;

		for (Point3d point: spherePoints){
			boolean is_accessible = true;
			Point3d test_point = new Point3d(point.x*radius + atom_i.x,
					point.y*radius + atom_i.y,
					point.z*radius + atom_i.z);

			int[] cycled_indices = new int[n_neighbor];
			int arind = 0;
			for (int ind=j_closest_neighbor;ind<n_neighbor;ind++) {
				cycled_indices[arind] = ind;
				arind++;
			}
			for (int ind=0;ind<j_closest_neighbor;ind++){
				cycled_indices[arind] = ind;
				arind++;
			}

			for (int j: cycled_indices) {
				Point3d atom_j = atomCoords[neighbor_indices.get(j)];
				double r = radii[neighbor_indices.get(j)] + probe;
				double diff_sq = test_point.distanceSquared(atom_j);
				if (diff_sq < r*r) {
					j_closest_neighbor = j;
					is_accessible = false;
					break;
				}
			}
			if (is_accessible) {
				n_accessible_point++;
			}
		}
		return cons*n_accessible_point*radius*radius;
	}

	/**
	 * Gets the radius for given amino acid and atom
	 * @param aa
	 * @param atom
	 * @return
	 */
	private static double getRadiusForAmino(AminoAcid amino, Atom atom) {

		if (atom.getElement().equals(Element.H)) return Element.H.getVDWRadius();
		// some unusual entries (e.g. 1tes) contain Deuterium atoms in standard aminoacids
		if (atom.getElement().equals(Element.D)) return Element.D.getVDWRadius();

		String atomCode = atom.getName();
		char aa = amino.getAminoType();

		// here we use the values that Chothia gives in his paper (as NACCESS does)
		if (atom.getElement()==Element.O) {
			return OXIGEN_VDW;
		}
		else if (atom.getElement()==Element.S) {
			return SULFUR_VDW;
		}
		else if (atom.getElement()==Element.N) {
			if (atomCode.equals("NZ")) return TETRAHEDRAL_NITROGEN_VDW; // tetrahedral Nitrogen
			return TRIGONAL_NITROGEN_VDW;								// trigonal Nitrogen
		}
		else if (atom.getElement()==Element.C) { // it must be a carbon
			if (atomCode.equals("C") ||
					atomCode.equals("CE1") || atomCode.equals("CE2") || atomCode.equals("CE3") ||
					atomCode.equals("CH2") ||
					atomCode.equals("CZ") || atomCode.equals("CZ2") || atomCode.equals("CZ3")) {
				return TRIGONAL_CARBON_VDW; 							// trigonal Carbon
			}
			else if (atomCode.equals("CA") || atomCode.equals("CB") ||
					atomCode.equals("CE") ||
					atomCode.equals("CG1") || atomCode.equals("CG2")) {
				return TETRAHEDRAL_CARBON_VDW;							// tetrahedral Carbon
			}
			// the rest of the cases (CD, CD1, CD2, CG) depend on amino acid
			else {
				switch (aa) {
				case 'F':
				case 'W':
				case 'Y':
				case 'H':
				case 'D':
				case 'N':
					return TRIGONAL_CARBON_VDW;

				case 'P':
				case 'K':
				case 'R':
				case 'M':
				case 'I':
				case 'L':
					return TETRAHEDRAL_CARBON_VDW;

				case 'Q':
				case 'E':
					if (atomCode.equals("CD")) return TRIGONAL_CARBON_VDW;
					else if (atomCode.equals("CG")) return TETRAHEDRAL_CARBON_VDW;

				default:
					logger.info("Unexpected carbon atom "+atomCode+" for aminoacid "+aa+", assigning its standard vdw radius");
					return Element.C.getVDWRadius();
				}
			}

			// not any of the expected atoms
		} else {
			// non standard aas, (e.g. MSE, LLP) will always have this problem,
			logger.info("Unexpected atom "+atomCode+" for aminoacid "+aa+ " ("+amino.getPDBName()+"), assigning its standard vdw radius");
			

			return atom.getElement().getVDWRadius();
		}
	}


	/**
	 * Gets the radius for given nucleotide atom
	 * @param atom
	 * @return
	 */
	private static double getRadiusForNucl(NucleotideImpl nuc, Atom atom) {

		if (atom.getElement().equals(Element.H)) return Element.H.getVDWRadius();
		if (atom.getElement().equals(Element.D)) return Element.D.getVDWRadius();

		if (atom.getElement()==Element.C) return NUC_CARBON_VDW;

		if (atom.getElement()==Element.N) return NUC_NITROGEN_VDW;

		if (atom.getElement()==Element.P) return PHOSPHOROUS_VDW;

		if (atom.getElement()==Element.O) return OXIGEN_VDW;

		logger.info("Unexpected atom "+atom.getName()+" for nucleotide "+nuc.getPDBName()+", assigning its standard vdw radius");
		return atom.getElement().getVDWRadius();
	}


	/**
	 * Gets the van der Waals radius of the given atom following the values defined by
	 * Chothia (1976) J.Mol.Biol.105,1-14
	 * NOTE: the vdw values defined by the paper assume no Hydrogens and thus "inflates"
	 * slightly the heavy atoms to account for Hydrogens. Thus this method cannot be used
	 * in a structure that contains Hydrogens!
	 *
	 * If atom is neither part of a nucleotide nor of a standard aminoacid,
	 * the default vdw radius for the element is returned. If atom is of
	 * unknown type (element) the vdw radius of {@link #Element().N} is returned
	 *
	 * @param atom
	 * @return
	 */
	public static double getRadius(Atom atom) {

		if (atom.getElement()==null) {
			logger.warn("Unrecognised atom "+atom.getName()+" with serial "+atom.getPDBserial()+
					", assigning the default vdw radius (Nitrogen vdw radius).");
			return Element.N.getVDWRadius();
		}

		Group res = atom.getGroup();

		if (res==null) {
			logger.warn("Unknown parent residue for atom "+atom.getName()+" with serial "+
					atom.getPDBserial()+", assigning its default vdw radius");
			return atom.getElement().getVDWRadius();
		}

		GroupType type = res.getType();

		if (type == GroupType.AMINOACID) return getRadiusForAmino(((AminoAcid)res), atom);

		if (type == GroupType.NUCLEOTIDE) return getRadiusForNucl((NucleotideImpl)res,atom);


		return atom.getElement().getVDWRadius();
	}

}
