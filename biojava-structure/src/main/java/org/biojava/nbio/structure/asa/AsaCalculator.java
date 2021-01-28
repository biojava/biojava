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
import org.biojava.nbio.structure.contact.Contact;
import org.biojava.nbio.structure.contact.Grid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


/**
 * Class to calculate Accessible Surface Areas based on
 * the rolling ball algorithm by Shrake and Rupley.
 * <p>
 * The code is adapted from a python implementation at http://boscoh.com/protein/asapy
 * (now source is available at https://github.com/boscoh/asa).
 * Thanks to Bosco K. Ho for a great piece of code and for his fantastic blog.
 * <p>
 * A few optimizations come from Eisenhaber et al, J Comp Chemistry 1994
 * (https://onlinelibrary.wiley.com/doi/epdf/10.1002/jcc.540160303)
 * <p>
 * See
 * Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms.
 * Lysozyme and Insulin." JMB (1973) 79:351-371.
 * Lee, B., and Richards, F.M. "The interpretation of Protein Structures: Estimation of
 * Static Accessibility" JMB (1971) 55:379-400
 *
 * @author Jose Duarte
 */
public class AsaCalculator {

	private static final Logger logger = LoggerFactory.getLogger(AsaCalculator.class);

	/**
	 * The default value for number of sphere points to sample.
	 * See this paper for a nice study on the effect of this parameter: https://f1000research.com/articles/5-189/v1
	 */
	public static final int DEFAULT_N_SPHERE_POINTS = 1000;
	public static final double DEFAULT_PROBE_SIZE = 1.4;
	public static final int DEFAULT_NTHREADS = 1;

	private static final boolean DEFAULT_USE_SPATIAL_HASHING = true;



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

		private final int i;
		private final double[] asas;

		private AsaCalcWorker(int i, double[] asas) {
			this.i = i;
			this.asas = asas;
		}

		@Override
		public void run() {
			asas[i] = calcSingleAsa(i);
		}
	}

	static class IndexAndDistance {
		final int index;
		final double dist;
		IndexAndDistance(int index, double dist) {
			this.index = index;
			this.dist = dist;
		}
	}


	private final Point3d[] atomCoords;
	private final Atom[] atoms;
	private final double[] radii;
	private final double probe;
	private final int nThreads;
	private Vector3d[] spherePoints;
	private double cons;
	private IndexAndDistance[][] neighborIndices;

	private boolean useSpatialHashingForNeighbors;

	/**
	 * Constructs a new AsaCalculator. Subsequently call {@link #calculateAsas()}
	 * or {@link #getGroupAsas()} to calculate the ASAs
	 * Only non-Hydrogen atoms are considered in the calculation.
	 * @param structure the structure, all non-H atoms will be used
	 * @param probe the probe size
	 * @param nSpherePoints the number of points to be used in generating the spherical
	 *                         dot-density, the more points the more accurate (and slower) calculation
	 * @param nThreads the number of parallel threads to use for the calculation
	 * @param hetAtoms if true HET residues are considered, if false they aren't, equivalent to
	 * NACCESS' -h option
	 */
	public AsaCalculator(Structure structure, double probe, int nSpherePoints, int nThreads, boolean hetAtoms) {
		this.atoms = StructureTools.getAllNonHAtomArray(structure, hetAtoms);
		this.atomCoords = Calc.atomsToPoints(atoms);
		this.probe = probe;
		this.nThreads = nThreads;

		this.useSpatialHashingForNeighbors = DEFAULT_USE_SPATIAL_HASHING;

		// initialising the radii by looking them up through AtomRadii
		radii = new double[atomCoords.length];
		for (int i=0;i<atomCoords.length;i++) {
			radii[i] = getRadius(atoms[i]);
		}

		initSpherePoints(nSpherePoints);
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

		this.useSpatialHashingForNeighbors = DEFAULT_USE_SPATIAL_HASHING;

		for (Atom atom:atoms) {
			if (atom.getElement()==Element.H)
				throw new IllegalArgumentException("Can't calculate ASA for an array that contains Hydrogen atoms ");
		}

		// initialising the radii by looking them up through AtomRadii
		radii = new double[atoms.length];
		for (int i=0;i<atoms.length;i++) {
			radii[i] = getRadius(atoms[i]);
		}

		initSpherePoints(nSpherePoints);
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

		this.useSpatialHashingForNeighbors = DEFAULT_USE_SPATIAL_HASHING;

		// initialising the radii to the given radius for all atoms
		radii = new double[atomCoords.length];
		for (int i=0;i<atomCoords.length;i++) {
			radii[i] = radius;
		}

		initSpherePoints(nSpherePoints);
	}

	private void initSpherePoints(int nSpherePoints) {

		logger.debug("Will use {} sphere points", nSpherePoints);

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

		TreeMap<ResidueNumber, GroupAsa> asas = new TreeMap<>();

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

		return asas.values().toArray(new GroupAsa[0]);
	}

	/**
	 * Calculates the Accessible Surface Areas for the atoms given in constructor and with parameters given.
	 * Beware that the parallel implementation is quite memory hungry. It scales well as long as there is
	 * enough memory available.
	 * @return an array with asa values corresponding to each atom of the input array
	 */
	public double[] calculateAsas() {

		double[] asas = new double[atomCoords.length];

		long start = System.currentTimeMillis();
		if (useSpatialHashingForNeighbors) {
			logger.debug("Will use spatial hashing to find neighbors");
			neighborIndices = findNeighborIndicesSpatialHashing();
		} else {
			logger.debug("Will not use spatial hashing to find neighbors");
			neighborIndices = findNeighborIndices();
		}
		long end = System.currentTimeMillis();
		logger.debug("Took {} s to find neighbors", (end-start)/1000.0);

		start = System.currentTimeMillis();
		if (nThreads<=1) { // (i.e. it will also be 1 thread if 0 or negative number specified)
			logger.debug("Will use 1 thread for ASA calculation");
			for (int i=0;i<atomCoords.length;i++) {
				asas[i] = calcSingleAsa(i);
			}

		} else {
			logger.debug("Will use {} threads for ASA calculation", nThreads);

			ExecutorService threadPool = Executors.newFixedThreadPool(nThreads);

			for (int i=0;i<atomCoords.length;i++) {
				threadPool.submit(new AsaCalcWorker(i,asas));
			}

			threadPool.shutdown();

			while (!threadPool.isTerminated());

		}
		end = System.currentTimeMillis();
		logger.debug("Took {} s to calculate all {} atoms ASAs (excluding neighbors calculation)", (end-start)/1000.0, atomCoords.length);

		return asas;
	}

	/**
	 * Set the useSpatialHashingForNeighbors flag to use spatial hashing to calculate neighbors (true) or all-to-all
	 * distance calculation (false). Default is {@value DEFAULT_USE_SPATIAL_HASHING}.
	 * Use for testing performance only.
	 * @param useSpatialHashingForNeighbors the flag
	 */
	void setUseSpatialHashingForNeighbors(boolean useSpatialHashingForNeighbors) {
		this.useSpatialHashingForNeighbors = useSpatialHashingForNeighbors;
	}

	/**
	 * Returns list of 3d coordinates of points on a unit sphere using the
	 * Golden Section Spiral algorithm.
	 * @param nSpherePoints the number of points to be used in generating the spherical dot-density
	 * @return the array of points as Vector3d objects
	 */
	private Vector3d[] generateSpherePoints(int nSpherePoints) {
		Vector3d[] points = new Vector3d[nSpherePoints];
		double inc = Math.PI * (3.0 - Math.sqrt(5.0));
		double offset = 2.0 / nSpherePoints;
		for (int k=0;k<nSpherePoints;k++) {
			double y = k * offset - 1.0 + (offset / 2.0);
			double r = Math.sqrt(1.0 - y*y);
			double phi = k * inc;
			points[k] = new Vector3d(Math.cos(phi)*r, y, Math.sin(phi)*r);
		}
		return points;
	}

	/**
	 * Returns the 2-dimensional array with neighbor indices for every atom.
	 * @return 2-dimensional array of size: n_atoms x n_neighbors_per_atom
	 */
	IndexAndDistance[][] findNeighborIndices() {

		// looking at a typical protein case, number of neighbours are from ~10 to ~50, with an average of ~30
		int initialCapacity = 60;

		IndexAndDistance[][] nbsIndices = new IndexAndDistance[atomCoords.length][];

		for (int k=0; k<atomCoords.length; k++) {
			double radius = radii[k] + probe + probe;

			List<IndexAndDistance> thisNbIndices = new ArrayList<>(initialCapacity);

			for (int i = 0; i < atomCoords.length; i++) {
				if (i == k) continue;

				double dist = atomCoords[i].distance(atomCoords[k]);

				if (dist < radius + radii[i]) {
					thisNbIndices.add(new IndexAndDistance(i, dist));
				}
			}

			IndexAndDistance[] indicesArray = thisNbIndices.toArray(new IndexAndDistance[0]);
			nbsIndices[k] = indicesArray;
		}
		return nbsIndices;
	}

	/**
	 * Returns the 2-dimensional array with neighbor indices for every atom,
	 * using spatial hashing to avoid all to all distance calculation.
	 * @return 2-dimensional array of size: n_atoms x n_neighbors_per_atom
	 */
	IndexAndDistance[][] findNeighborIndicesSpatialHashing() {

		// looking at a typical protein case, number of neighbours are from ~10 to ~50, with an average of ~30
		int initialCapacity = 60;

		List<Contact> contactList = calcContacts();
		Map<Integer, List<IndexAndDistance>> indices = new HashMap<>(atomCoords.length);
		for (Contact contact : contactList) {
			// note contacts are stored 1-way only, with j>i
			int i = contact.getI();
			int j = contact.getJ();

			List<IndexAndDistance> iIndices;
			List<IndexAndDistance> jIndices;
			if (!indices.containsKey(i)) {
				iIndices = new ArrayList<>(initialCapacity);
				indices.put(i, iIndices);
			} else {
				iIndices = indices.get(i);
			}
			if (!indices.containsKey(j)) {
				jIndices = new ArrayList<>(initialCapacity);
				indices.put(j, jIndices);
			} else {
				jIndices = indices.get(j);
			}

			double radius = radii[i] + probe + probe;
			double dist = contact.getDistance();
			if (dist < radius + radii[j]) {
				iIndices.add(new IndexAndDistance(j, dist));
				jIndices.add(new IndexAndDistance(i, dist));
			}
		}

		// convert map to array for fast access
		IndexAndDistance[][] nbsIndices = new IndexAndDistance[atomCoords.length][];
		for (Map.Entry<Integer, List<IndexAndDistance>> entry : indices.entrySet()) {
			List<IndexAndDistance> list = entry.getValue();
			IndexAndDistance[] indexAndDistances = list.toArray(new IndexAndDistance[0]);
			nbsIndices[entry.getKey()] = indexAndDistances;
		}

		// important: some atoms might have no neighbors at all: we need to initialise to empty arrays
		for (int i=0; i<nbsIndices.length; i++) {
			if (nbsIndices[i] == null) {
				nbsIndices[i] = new IndexAndDistance[0];
			}
		}

		return nbsIndices;
	}

	Point3d[] getAtomCoords() {
		return atomCoords;
	}

	private List<Contact> calcContacts() {
		if (atomCoords.length == 0)
			return new ArrayList<>();
		double maxRadius = 0;
		OptionalDouble optionalDouble = Arrays.stream(radii).max();
		if (optionalDouble.isPresent())
			maxRadius = optionalDouble.getAsDouble();
		double cutoff = maxRadius + maxRadius + probe + probe;
		logger.debug("Max radius is {}, cutoff is {}", maxRadius, cutoff);
		Grid grid = new Grid(cutoff);
		grid.addCoords(atomCoords);
		return grid.getIndicesContacts();
	}

	private double calcSingleAsa(int i) {
		Point3d atom_i = atomCoords[i];

		int n_neighbor = neighborIndices[i].length;
		IndexAndDistance[] neighbor_indices = neighborIndices[i];
		// Sorting by closest to farthest away neighbors achieves faster runtimes when checking for occluded
		// sphere sample points below. This follows the ideas exposed in
		// Eisenhaber et al, J Comp Chemistry 1994 (https://onlinelibrary.wiley.com/doi/epdf/10.1002/jcc.540160303)
		// This is essential for performance. In my tests this brings down the number of occlusion checks in loop below to
		// an average of n_sphere_points/10 per atom i, producing ~ x4 performance gain overall
		Arrays.sort(neighbor_indices, Comparator.comparingDouble(o -> o.dist));

		double radius_i = probe + radii[i];

		int n_accessible_point = 0;
		// purely for debugging
		int[] numDistsCalced = null;
		if (logger.isDebugEnabled()) numDistsCalced = new int[n_neighbor];

		// now we precalculate anything depending only on i,j in equation 3 in Eisenhaber 1994
		double[] sqRadii = new double[n_neighbor];
		Vector3d[] aj_minus_ais = new Vector3d[n_neighbor];
		for (int nbArrayInd =0; nbArrayInd<n_neighbor; nbArrayInd++) {
			int j = neighbor_indices[nbArrayInd].index;
			double dist = neighbor_indices[nbArrayInd].dist;
			double radius_j = radii[j] + probe;
			// see equation 3 in Eisenhaber 1994
			sqRadii[nbArrayInd] = (dist*dist + radius_i*radius_i - radius_j*radius_j)/(2*radius_i);
			Vector3d aj_minus_ai = new Vector3d(atomCoords[j]);
			aj_minus_ai.sub(atom_i);
			aj_minus_ais[nbArrayInd] = aj_minus_ai;
		}

		for (Vector3d point: spherePoints){
			boolean is_accessible = true;

			// note that the neighbors are sorted by distance, achieving optimal performance in this inner loop
			// See Eisenhaber et al, J Comp Chemistry 1994

			for (int nbArrayInd =0; nbArrayInd<n_neighbor; nbArrayInd++) {

				// see equation 3 in Eisenhaber 1994. This is slightly more efficient than
				// calculating distances to the actual sphere points on atom_i (which would be obtained with:
				// Point3d test_point = new Point3d(point.x*radius + atom_i.x,point.y*radius + atom_i.y,point.z*radius + atom_i.z))
				double dotProd = aj_minus_ais[nbArrayInd].dot(point);

				if (numDistsCalced!=null) numDistsCalced[nbArrayInd]++;

				if (dotProd > sqRadii[nbArrayInd]) {
					is_accessible = false;
					break;
				}
			}
			if (is_accessible) {
				n_accessible_point++;
			}
		}

		// purely for debugging
		if (numDistsCalced!=null) {
			int sum = 0;
			for (int numDistCalcedForJ : numDistsCalced) sum += numDistCalcedForJ;
			logger.debug("Number of sample points distances calculated for neighbors of i={} : average {}, all {}", i, (double) sum / (double) n_neighbor, numDistsCalced);
		}

		return cons*n_accessible_point*radius_i*radius_i;
	}

	/**
	 * Gets the radius for given amino acid and atom
	 * @param amino
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
					logger.info("Unexpected carbon atom {} for aminoacid {}, assigning its standard vdw radius", atomCode, aa);
					return Element.C.getVDWRadius();
				}
			}

			// not any of the expected atoms
		} else {
			// non standard aas, (e.g. MSE, LLP) will always have this problem,
			logger.debug("Unexpected atom {} for aminoacid {} ({}), assigning its standard vdw radius", atomCode, aa, amino.getPDBName());

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
	 * unknown type (element) the vdw radius of {@link Element().N} is returned
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
