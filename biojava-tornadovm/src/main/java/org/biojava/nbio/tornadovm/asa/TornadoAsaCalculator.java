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
package org.biojava.nbio.tornadovm.asa;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.asa.AsaCalculator;
import org.biojava.nbio.structure.asa.GroupAsa;
import org.biojava.nbio.tornadovm.TornadoSupport;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.manchester.tornado.api.GridScheduler;
import uk.ac.manchester.tornado.api.KernelContext;
import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.WorkerGrid;
import uk.ac.manchester.tornado.api.WorkerGrid1D;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;
import uk.ac.manchester.tornado.api.types.arrays.DoubleArray;
import uk.ac.manchester.tornado.api.types.arrays.FloatArray;
import uk.ac.manchester.tornado.api.types.arrays.IntArray;

/**
 * Drop-in replacement for {@link AsaCalculator} that runs the neighbour search and the sphere-point occlusion test
 * on the GPU.
 * <p>
 * Constructors and results match {@link AsaCalculator}: the same atoms, radii, sphere points and neighbour criterion
 * are used, in double precision and with the same operation order. Small structures, or runs without a TornadoVM
 * runtime, are delegated to {@link AsaCalculator}.
 */
public class TornadoAsaCalculator {

	private static final Logger logger = LoggerFactory.getLogger(TornadoAsaCalculator.class);

	/** Below this many atoms the multi-threaded CPU implementation is used. */
	public static final int DEFAULT_MIN_GPU_ATOMS = 1000;

	private final Point3d[] atomCoords;
	private final Atom[] atoms;
	private final double[] radii;
	private final double probe;
	private final int nSpherePoints;
	private final double cons;
	private final int nThreads;

	private int minGpuAtoms = DEFAULT_MIN_GPU_ATOMS;

	private long lastHostNanos;
	private int lastRecomputedAtoms;
	private long lastGpuNanos;

	/** @see AsaCalculator#AsaCalculator(Structure, double, int, int, boolean, int) */
	public TornadoAsaCalculator(Structure structure, double probe, int nSpherePoints, int nThreads, boolean hetAtoms,
			int modelNr) {
		this(StructureTools.getAllNonHAtomArray(structure, hetAtoms, modelNr), probe, nSpherePoints, nThreads);
	}

	/** @see AsaCalculator#AsaCalculator(Structure, double, int, int, boolean) */
	public TornadoAsaCalculator(Structure structure, double probe, int nSpherePoints, int nThreads,
			boolean hetAtoms) {
		this(structure, probe, nSpherePoints, nThreads, hetAtoms, 0);
	}

	/** @see AsaCalculator#AsaCalculator(Atom[], double, int, int) */
	public TornadoAsaCalculator(Atom[] atoms, double probe, int nSpherePoints, int nThreads) {
		for (Atom atom : atoms) {
			// as AsaCalculator
			if (atom.getElement() == Element.H) {
				throw new IllegalArgumentException("Can't calculate ASA for an array that contains Hydrogen atoms ");
			}
		}
		this.nThreads = nThreads;
		this.atoms = atoms;
		this.atomCoords = Calc.atomsToPoints(atoms);
		this.probe = probe;
		this.nSpherePoints = nSpherePoints;
		this.cons = 4.0 * Math.PI / nSpherePoints;
		this.radii = new double[atoms.length];
		for (int i = 0; i < atoms.length; i++) {
			radii[i] = AsaCalculator.getRadius(atoms[i]);
		}
	}

	/** @see AsaCalculator#AsaCalculator(Point3d[], double, int, int, double) */
	public TornadoAsaCalculator(Point3d[] atomCoords, double probe, int nSpherePoints, int nThreads, double radius) {
		this.nThreads = nThreads;
		this.atoms = null;
		this.atomCoords = atomCoords;
		this.probe = probe;
		this.nSpherePoints = nSpherePoints;
		this.cons = 4.0 * Math.PI / nSpherePoints;
		this.radii = new double[atomCoords.length];
		Arrays.fill(radii, radius);
	}

	/** Sets the minimum number of atoms for which the GPU is used. */
	public void setMinGpuAtoms(int minGpuAtoms) {
		this.minGpuAtoms = minGpuAtoms;
	}

	/** @return number of atoms the last single-precision GPU run left to the exact CPU path */
	public int getLastRecomputedAtoms() {
		return lastRecomputedAtoms;
	}

	/** @return host time of the last GPU run spent bucketing atoms into grid cells, in ns */
	public long getLastHostNanos() {
		return lastHostNanos;
	}

	/** @return time of the last GPU run spent in TornadoVM (transfers, kernels, and JIT if cold), in ns */
	public long getLastGpuNanos() {
		return lastGpuNanos;
	}

	/** @see AsaCalculator#getGroupAsas() */
	public GroupAsa[] getGroupAsas() {
		TreeMap<ResidueNumber, GroupAsa> asas = new TreeMap<>();
		double[] asasPerAtom = calculateAsas();
		for (int i = 0; i < atomCoords.length; i++) {
			Group g = atoms[i].getGroup();
			GroupAsa groupAsa = asas.get(g.getResidueNumber());
			if (groupAsa == null) {
				groupAsa = new GroupAsa(g);
				asas.put(g.getResidueNumber(), groupAsa);
			}
			groupAsa.addAtomAsaU(asasPerAtom[i]);
		}
		return asas.values().toArray(new GroupAsa[0]);
	}

	/**
	 * Same result as {@link AsaCalculator#calculateAsas()}: on the GPU for large structures, else with
	 * {@link CpuAsa} (single-threaded if nThreads &lt;= 1, as AsaCalculator).
	 */
	public double[] calculateAsas() {
		if (TornadoSupport.useGpu("asa", atomCoords.length, minGpuAtoms)) {
			try {
				return calculateAsasGpu();
			} catch (RuntimeException | Error e) {
				TornadoSupport.disable("asa", e);
			}
		}
		return calculateAsasCpu();
	}

	/** Runs the CPU path ({@link CpuAsa}). */
	public double[] calculateAsasCpu() {
		return CpuAsa.calculate(atomCoords, radii, probe, nSpherePoints, nThreads > 1);
	}

	/** Runs the GPU path unconditionally (no threshold, no fallback). */
	public double[] calculateAsasGpu() {
		int n = atomCoords.length;
		double[] asas = new double[n];
		if (n == 0) {
			return asas;
		}
		long t0 = System.nanoTime();
		CellGrid cells = new CellGrid();
		long t1 = System.nanoTime();

		IntArray counts = runKernel(cells);
		int ambiguous = 0;
		for (int i = 0; i < n; i++) {
			if (counts.get(i) == -1) {
				logger.info("Atom {} has more than {} neighbours, using the CPU implementation", i,
						AsaKernels.LOCAL_CAPACITY);
				return calculateAsasCpu();
			}
			if (counts.get(i) == -2) {
				ambiguous++;
			}
		}
		long t2 = System.nanoTime();

		for (int i = 0; i < n; i++) {
			double radius = probe + radii[i];
			asas[i] = cons * counts.get(i) * radius * radius;
		}
		if (ambiguous > 0) {
			// single-precision kernel: redo the atoms with a sphere point near an occlusion boundary exactly
			int[] redo = new int[ambiguous];
			for (int i = 0, k = 0; i < n; i++) {
				if (counts.get(i) == -2) {
					redo[k++] = i;
				}
			}
			double[] exact = CpuAsa.calculate(atomCoords, radii, probe, nSpherePoints, nThreads > 1, redo);
			for (int i : redo) {
				asas[i] = exact[i];
			}
			lastRecomputedAtoms = ambiguous;
		} else {
			lastRecomputedAtoms = 0;
		}
		lastHostNanos = t1 - t0;
		lastGpuNanos = t2 - t1;
		return asas;
	}

	private IntArray runKernel(CellGrid cells) {
		boolean fp32 = useFp32();
		Engine engine;
		synchronized (ENGINES) {
			engine = ENGINES.computeIfAbsent(fp32 ? -nSpherePoints : nSpherePoints, key -> new Engine(nSpherePoints, fp32));
		}
		try {
			return engine.run(atomCoords, radii, probe, cells);
		} catch (RuntimeException e) {
			if (!fp32 && isFp64Unsupported(e)) {
				logger.info("The device has no FP64: using the single-precision ASA kernel with exact CPU re-checks");
				fp64Unsupported = true;
				return runKernel(cells);
			}
			throw e;
		}
	}

	/** Set once a device without FP64 has been seen. */
	private static volatile boolean fp64Unsupported;

	/** {@code biojava.tornado.asa.precision}: auto (default: FP64, FP32 if the device has no FP64), fp32 or fp64. */
	private static boolean useFp32() {
		String precision = System.getProperty("biojava.tornado.asa.precision", "auto");
		return "fp32".equalsIgnoreCase(precision) || ("auto".equalsIgnoreCase(precision) && fp64Unsupported);
	}

	private static boolean isFp64Unsupported(Throwable e) {
		for (Throwable t = e; t != null; t = t.getCause()) {
			if (t.getClass().getSimpleName().contains("FP64")) {
				return true;
			}
		}
		return false;
	}

	/** One persistent execution plan per number of sphere points, shared by all calculators. */
	private static final Map<Integer, Engine> ENGINES = new HashMap<>();

	/**
	 * A TornadoVM execution plan whose buffers are sized for a capacity of atoms and grid cells rather than for one
	 * structure, so that it is compiled once and re-executed for every structure that fits. It is rebuilt with a
	 * larger capacity when a structure does not fit.
	 */
	private static final class Engine {
		private final int nPoints;
		private final boolean fp32;
		private final DoubleArray spherePoints;
		private final FloatArray spherePointsF;
		private FloatArray coordsF;
		private FloatArray radiiF;
		private FloatArray gridF;
		private IntArray atomCell;
		private int atomCapacity;
		private int cellCapacity;
		private DoubleArray coords;
		private DoubleArray radii;
		private IntArray cellStart;
		private IntArray cellAtoms;
		private DoubleArray grid;
		private IntArray dims;
		private IntArray counts;
		private WorkerGrid worker;
		private GridScheduler gridScheduler;
		private TornadoExecutionPlan plan;

		Engine(int nPoints, boolean fp32) {
			this.nPoints = nPoints;
			this.fp32 = fp32;
			double[] points = generateSpherePoints(nPoints);
			this.spherePoints = DoubleArray.fromArray(points);
			float[] pointsF = new float[points.length];
			for (int k = 0; k < points.length; k++) {
				pointsF[k] = (float) points[k];
			}
			this.spherePointsF = FloatArray.fromArray(pointsF);
		}

		synchronized IntArray run(Point3d[] atomCoords, double[] atomRadii, double probe, CellGrid cells) {
			int n = atomCoords.length;
			int nCells = cells.nx * cells.ny * cells.nz;
			if (plan == null || n > atomCapacity || nCells + 1 > cellCapacity) {
				build(Math.max(n, atomCapacity), Math.max(nCells + 1, cellCapacity));
			}
			if (fp32) {
				// coordinates relative to the atom's grid cell, computed in double: single precision then keeps
				// ~5e-7 A of accuracy whatever the size of the structure
				for (int i = 0; i < n; i++) {
					int cx = cells.atomCell[3 * i];
					int cy = cells.atomCell[3 * i + 1];
					int cz = cells.atomCell[3 * i + 2];
					coordsF.set(3 * i, (float) (atomCoords[i].x - (cells.minX + cx * cells.cellSize)));
					coordsF.set(3 * i + 1, (float) (atomCoords[i].y - (cells.minY + cy * cells.cellSize)));
					coordsF.set(3 * i + 2, (float) (atomCoords[i].z - (cells.minZ + cz * cells.cellSize)));
					radiiF.set(i, (float) atomRadii[i]);
					atomCell.set(3 * i, cx);
					atomCell.set(3 * i + 1, cy);
					atomCell.set(3 * i + 2, cz);
				}
				gridF.set(0, 0f);
				gridF.set(1, 0f);
				gridF.set(2, 0f);
				gridF.set(3, (float) cells.cellSize);
				gridF.set(4, (float) probe);
			} else {
				for (int i = 0; i < n; i++) {
					coords.set(3 * i, atomCoords[i].x);
					coords.set(3 * i + 1, atomCoords[i].y);
					coords.set(3 * i + 2, atomCoords[i].z);
					radii.set(i, atomRadii[i]);
				}
				grid.set(0, cells.minX);
				grid.set(1, cells.minY);
				grid.set(2, cells.minZ);
				grid.set(3, cells.cellSize);
				grid.set(4, probe);
			}
			for (int i = 0; i < n; i++) {
				cellAtoms.set(i, cells.cellAtoms[i]);
			}
			for (int c = 0; c <= nCells; c++) {
				cellStart.set(c, cells.cellStart[c]);
			}
			dims.set(0, n);
			dims.set(1, cells.nx);
			dims.set(2, cells.ny);
			dims.set(3, cells.nz);
			worker.setGlobalWork((long) n * AsaKernels.GROUP_SIZE, 1, 1);

			long t0 = System.nanoTime();
			plan.withGridScheduler(gridScheduler).execute();
			if (Boolean.getBoolean("biojava.tornado.trace")) {
				System.err.printf("asa: n=%d capacity=%d execute %.1f ms%n", n, atomCapacity,
						(System.nanoTime() - t0) / 1e6);
			}
			IntArray result = new IntArray(n);
			for (int i = 0; i < n; i++) {
				result.set(i, counts.get(i));
			}
			return result;
		}

		private void build(int minAtoms, int minCells) {
			if (plan != null) {
				try {
					plan.close();
				} catch (Exception e) {
					logger.warn("Could not release the previous ASA execution plan", e);
				}
			}
			atomCapacity = Math.max(1 << 14, Integer.highestOneBit(minAtoms - 1) << 1);
			cellCapacity = Math.max(1 << 16, Integer.highestOneBit(minCells - 1) << 1);
			cellAtoms = new IntArray(atomCapacity);
			cellStart = new IntArray(cellCapacity);
			dims = new IntArray(4);
			counts = new IntArray(atomCapacity);
			KernelContext ctx = new KernelContext();

			TaskGraph taskGraph;
			if (fp32) {
				coordsF = new FloatArray(3 * atomCapacity);
				radiiF = new FloatArray(atomCapacity);
				atomCell = new IntArray(3 * atomCapacity);
				gridF = new FloatArray(5);
				taskGraph = new TaskGraph("asa")
						.transferToDevice(DataTransferMode.FIRST_EXECUTION, spherePointsF)
						.transferToDevice(DataTransferMode.EVERY_EXECUTION, coordsF, radiiF, cellStart, cellAtoms,
								atomCell, gridF, dims)
						.task("fused", AsaKernels::fusedAsaFp32, ctx, coordsF, radiiF, cellStart, cellAtoms,
								atomCell, gridF, dims, spherePointsF, counts, nPoints)
						.transferToHost(DataTransferMode.EVERY_EXECUTION, counts);
			} else {
				coords = new DoubleArray(3 * atomCapacity);
				radii = new DoubleArray(atomCapacity);
				grid = new DoubleArray(5);
				taskGraph = new TaskGraph("asa")
						.transferToDevice(DataTransferMode.FIRST_EXECUTION, spherePoints)
						.transferToDevice(DataTransferMode.EVERY_EXECUTION, coords, radii, cellStart, cellAtoms, grid,
								dims)
						.task("fused", AsaKernels::fusedAsa, ctx, coords, radii, cellStart, cellAtoms, grid, dims,
								spherePoints, counts, nPoints)
						.transferToHost(DataTransferMode.EVERY_EXECUTION, counts);
			}
			worker = new WorkerGrid1D(AsaKernels.GROUP_SIZE);
			worker.setLocalWork(AsaKernels.GROUP_SIZE, 1, 1);
			gridScheduler = new GridScheduler("asa.fused", worker);
			plan = new TornadoExecutionPlan(taskGraph.snapshot());
		}
	}

	/**
	 * A uniform grid of cubic cells whose side is at least the largest possible neighbour distance, so that all
	 * neighbours of an atom are in its own or the 26 adjacent cells. Atoms are bucketed by cell with a counting sort.
	 */
	private final class CellGrid {
		final double minX, minY, minZ, cellSize;
		final int nx, ny, nz;
		final int[] cellStart;
		final int[] cellAtoms;
		final int[] atomCell;

		CellGrid() {
			int n = atomCoords.length;
			double maxRadius = Arrays.stream(radii).max().orElse(0);
			double loX = Double.MAX_VALUE, loY = Double.MAX_VALUE, loZ = Double.MAX_VALUE;
			double hiX = -Double.MAX_VALUE, hiY = -Double.MAX_VALUE, hiZ = -Double.MAX_VALUE;
			for (Point3d p : atomCoords) {
				loX = Math.min(loX, p.x);
				loY = Math.min(loY, p.y);
				loZ = Math.min(loZ, p.z);
				hiX = Math.max(hiX, p.x);
				hiY = Math.max(hiY, p.y);
				hiZ = Math.max(hiZ, p.z);
			}
			// same cutoff as AsaCalculator.calcContacts, never smaller than any neighbour distance
			double side = Math.max(maxRadius + maxRadius + probe + probe, 1e-3);
			// cap the number of cells for sparse structures
			while (cellCount(hiX - loX, hiY - loY, hiZ - loZ, side) > Math.max(8L * n, 1L << 20)) {
				side *= 1.5;
			}
			minX = loX;
			minY = loY;
			minZ = loZ;
			cellSize = side;
			nx = (int) ((hiX - loX) / side) + 1;
			ny = (int) ((hiY - loY) / side) + 1;
			nz = (int) ((hiZ - loZ) / side) + 1;

			int nCells = nx * ny * nz;
			int[] cellOf = new int[n];
			int[] start = new int[nCells + 1];
			atomCell = new int[3 * n];
			for (int i = 0; i < n; i++) {
				Point3d p = atomCoords[i];
				int cx = (int) ((p.x - minX) / cellSize);
				int cy = (int) ((p.y - minY) / cellSize);
				int cz = (int) ((p.z - minZ) / cellSize);
				atomCell[3 * i] = cx;
				atomCell[3 * i + 1] = cy;
				atomCell[3 * i + 2] = cz;
				cellOf[i] = (cx * ny + cy) * nz + cz;
				start[cellOf[i] + 1]++;
			}
			for (int c = 0; c < nCells; c++) {
				start[c + 1] += start[c];
			}
			int[] fill = Arrays.copyOf(start, nCells);
			int[] sorted = new int[n];
			for (int i = 0; i < n; i++) {
				sorted[fill[cellOf[i]]++] = i;
			}
			cellStart = start;
			cellAtoms = sorted;
		}

		private long cellCount(double dx, double dy, double dz, double side) {
			return ((long) (dx / side) + 1) * ((long) (dy / side) + 1) * ((long) (dz / side) + 1);
		}
	}

	/** Golden Section Spiral points, as {@code AsaCalculator.generateSpherePoints}. */
	static double[] generateSpherePoints(int nSpherePoints) {
		double[] points = new double[3 * nSpherePoints];
		double inc = Math.PI * (3.0 - Math.sqrt(5.0));
		double offset = 2.0 / nSpherePoints;
		for (int k = 0; k < nSpherePoints; k++) {
			double y = k * offset - 1.0 + (offset / 2.0);
			double r = Math.sqrt(1.0 - y * y);
			double phi = k * inc;
			points[3 * k] = Math.cos(phi) * r;
			points[3 * k + 1] = y;
			points[3 * k + 2] = Math.sin(phi) * r;
		}
		return points;
	}
}
