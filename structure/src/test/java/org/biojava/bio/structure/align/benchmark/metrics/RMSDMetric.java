package org.biojava.bio.structure.align.benchmark.metrics;


import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.benchmark.MultipleAlignment;
import org.biojava.bio.structure.align.model.AFPChain;


public abstract class RMSDMetric {
	public static class Reference extends Metric {

		@Override
		public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2) {
			try {

				List<Atom[]> structures = new ArrayList<Atom[]>(2);
				structures.add(ca1);
				structures.add(ca2);
				int[][] optAln = reference.getAlignmentMatrix(structures);
				
				// Create new arrays for the subset of atoms in the alignment.
				Atom[] ca1aligned = new Atom[reference.size()];
				Atom[] ca2aligned = new Atom[reference.size()];
				int pos=0;
				for(int i=0;i<optAln[0].length;i++) {
					ca1aligned[pos] = ca1[optAln[0][pos]];
					ca2aligned[pos] = (Atom) ca2[optAln[1][pos]].clone();
					pos++;

				}

				return SVDSuperimposer.getRMS(ca1aligned, ca2aligned);
			} catch (StructureException e) {
				e.printStackTrace();
				return Double.NaN;
			}
		}

		@Override
		public String getName() {
			return "Ref_RMSD";
		}

	}

	public static class Alignment extends Metric {

		@Override
		public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2) {
			// Create new arrays for the subset of atoms in the alignment.
			Atom[] ca1aligned = new Atom[align.getOptLength()];
			Atom[] ca2aligned = new Atom[align.getOptLength()];
			int pos=0;
			int[] blockLens = align.getOptLen();
			int[][][] optAln = align.getOptAln();
			//for(int block=0;block< align.getBlockNum();block++) {
			assert(align.getBlockNum() == optAln.length);
			for(int block=0;block< optAln.length;block++) {
				//for(int i=0;i<blockLens[block];i++) {
				assert(blockLens[block] == optAln[block][0].length);
				for(int i=0;i<optAln[block][0].length;i++) {
					ca1aligned[pos] = ca1[optAln[block][0][i]];
					ca2aligned[pos] = (Atom) ca2[optAln[block][1][i]].clone();
					pos++;
				}
			}

			try {
				return SVDSuperimposer.getRMS(ca1aligned, ca2aligned);
			} catch (StructureException e) {
				e.printStackTrace();
				return Double.NaN;
			}
		}

		@Override
		public String getName() {
			return "Aln_RMSD";
		}
		
		@Override
		public String format(double result) {
			return String.format("%.3f", result);
		}
	}
}
