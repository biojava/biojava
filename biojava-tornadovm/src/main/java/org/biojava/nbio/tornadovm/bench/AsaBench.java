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
package org.biojava.nbio.tornadovm.bench;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.asa.AsaCalculator;
import org.biojava.nbio.tornadovm.asa.TornadoAsaCalculator;

/**
 * ASA parity check and timing: BioJava (1 thread, all cores) vs TornadoVM.
 * Usage: AsaBench PDBID [PDBID...]; property bench.reps (default 5).
 */
public class AsaBench {

	public static void main(String[] args) throws Exception {
		int reps = Integer.getInteger("bench.reps", 5);
		int cores = Runtime.getRuntime().availableProcessors();
		System.out.printf("%-6s %8s | %10s %10s %10s | %10s | %9s %9s %9s | %s%n", "pdb", "atoms", "bj1 ms",
				"bjN ms", "leanN ms", "gpu ms", "vs bj1", "vs bjN", "vs leanN", "parity");
		for (String id : args) {
			Structure s = Structures.load(id);
			Atom[] atoms = StructureTools.getAllNonHAtomArray(s, false, 0);

			double[] ref = null;
			double cpu1 = Double.MAX_VALUE, cpuN = Double.MAX_VALUE, gpu = Double.MAX_VALUE;
			for (int r = 0; r < reps; r++) {
				long t = System.nanoTime();
				ref = new AsaCalculator(atoms, AsaCalculator.DEFAULT_PROBE_SIZE, AsaCalculator.DEFAULT_N_SPHERE_POINTS, 1).calculateAsas();
				cpu1 = Math.min(cpu1, (System.nanoTime() - t) / 1e6);
			}
			for (int r = 0; r < reps; r++) {
				long t = System.nanoTime();
				new AsaCalculator(atoms, AsaCalculator.DEFAULT_PROBE_SIZE, AsaCalculator.DEFAULT_N_SPHERE_POINTS, cores).calculateAsas();
				cpuN = Math.min(cpuN, (System.nanoTime() - t) / 1e6);
			}
			double[] radii = new double[atoms.length];
			for (int i = 0; i < atoms.length; i++) radii[i] = AsaCalculator.getRadius(atoms[i]);
			javax.vecmath.Point3d[] pts = org.biojava.nbio.structure.Calc.atomsToPoints(atoms);
			String dump = System.getProperty("bench.dump");
			if (dump != null) {
				// atoms and BioJava radii, for external tools (e.g. FreeSASA's C API)
				try (java.io.PrintWriter w = new java.io.PrintWriter(new java.io.File(dump, id + ".xyzr"))) {
					w.println(atoms.length);
					for (int i = 0; i < atoms.length; i++) w.printf(java.util.Locale.ROOT, "%.3f %.3f %.3f %.4f%n", pts[i].x, pts[i].y, pts[i].z, radii[i]);
				}
			}
			double lean = Double.MAX_VALUE;
			double[] leanAsas = null;
			for (int r = 0; r < reps + 1; r++) {
				long t = System.nanoTime();
				leanAsas = org.biojava.nbio.tornadovm.asa.CpuAsa.calculate(pts, radii, AsaCalculator.DEFAULT_PROBE_SIZE, AsaCalculator.DEFAULT_N_SPHERE_POINTS, true);
				if (r > 0) lean = Math.min(lean, (System.nanoTime() - t) / 1e6);
			}
			if (!java.util.Arrays.equals(ref, leanAsas)) System.out.println("  lean CPU differs from BioJava: " + parity(ref, leanAsas));
			double[] got = null;
			double nb = 0, tv = 0;
			for (int r = 0; r < reps + 1; r++) {
				TornadoAsaCalculator calc = new TornadoAsaCalculator(atoms, AsaCalculator.DEFAULT_PROBE_SIZE, AsaCalculator.DEFAULT_N_SPHERE_POINTS, cores);
				long t = System.nanoTime();
				got = calc.calculateAsasGpu();
				double ms = (System.nanoTime() - t) / 1e6;
				if (r == 0) {
					System.out.printf("  [%s cold GPU run incl. JIT: %.1f ms, atoms re-checked on CPU: %d]%n", id, ms, calc.getLastRecomputedAtoms());
				} else if (ms < gpu) {
					gpu = ms;
					nb = calc.getLastHostNanos() / 1e6;
					tv = calc.getLastGpuNanos() / 1e6;
				}
			}
			System.out.printf("%-6s %8d | %10.1f %10.1f %10.1f | %10.1f | %8.1fx %8.1fx %8.1fx | %s%n", id, atoms.length,
					cpu1, cpuN, lean, gpu, cpu1 / gpu, cpuN / gpu, lean / gpu, parity(ref, got));
		}
	}

	static String parity(double[] ref, double[] got) {
		double totRef = 0, totGot = 0, maxAbs = 0;
		int diff = 0;
		for (int i = 0; i < ref.length; i++) {
			totRef += ref[i];
			totGot += got[i];
			double d = Math.abs(ref[i] - got[i]);
			if (d > 1e-9) diff++;
			maxAbs = Math.max(maxAbs, d);
		}
		return String.format("total %.2f vs %.2f (rel %.1e), atoms differing %d, max |d| %.3f A^2", totRef, totGot,
				Math.abs(totRef - totGot) / totRef, diff, maxAbs);
	}
}
